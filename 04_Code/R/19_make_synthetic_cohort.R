## =============================================================================
## 19_make_synthetic_cohort.R
##
## Generates a fully synthetic cohort with the EXACT column schema of the
## confidential real data, so that the published analysis pipeline runs
## end-to-end without access to patient records.
##
## WHY THIS EXISTS. The Data and Code Availability statement promises a synthetic
## dataset "of the same structure" that "allows the full analysis pipeline to be
## executed end-to-end". The pre-existing simulated_*.csv files do not satisfy
## that: their schema differs from real_*.csv and they omit columns the pipeline
## consumes -- log_bcrabl and log_abl in the longitudinal file (used by
## 03_build_stan_data.R), and n_visits, has_complete_covariates,
## baseline/last/min_log_mrd and most ph_* columns in the patient file (used by
## 04_prepare_joint_model_data.R). Running the pipeline on them fails.
##
## The files written here carry every column the real files carry, in the same
## order, with the same types, so `prefix = "synthetic"` is a drop-in for
## `prefix = "real"` throughout.
##
## GENERATIVE MODEL. Parameters are the posterior means of the fitted joint model
## (Table 5): a quadratic trajectory in log(1+t) with independent random
## intercept and slope, Gaussian observation error, left-censoring at the assay
## floor, and events from the complementary log-log interval hazard. Visit
## spacing reproduces the cohort's marginal features (median gap ~6 months,
## 87 patients, ~5-6 visits each). NO PATIENT DATA IS READ OR USED.
##
## The synthetic data are for verifying that the code runs and produces output of
## the right shape. They are not the study data and will not reproduce the
## reported estimates; the DOI'd archive records the real estimates in the
## results files.
##
## RUN from the project root:
##   Rscript 04_Code\R\19_make_synthetic_cohort.R
## OUTPUT: 03_Data/Processed/synthetic_{longitudinal,interval_survival,patient_level}_analysis.csv
## =============================================================================

set.seed(20260801)
root <- normalizePath(".", winslash = "/")
outdir <- file.path(root, "03_Data", "Processed")
if (!dir.exists(outdir)) stop("Run from the project root; ", outdir, " not found.")

## ---- parameters: posterior means of the fitted model (Table 5) --------------
TH <- list(beta0 = -2.107, beta_time = -3.561, beta_time2 = 0.502,
           beta_bm = 0.607, sigma_y = 1.813, tau_b0 = 0.968, tau_b1 = 2.203,
           gamma0 = -4.036, gamma_time = -0.228, gamma_gap = -1.831,
           alpha_mrd = -1.275)
FLOOR <- -5.0          # assay floor on the log10 scale
CD    <- -4.5          # DMR threshold
N     <- 87L           # cohort size

long <- list(); iv <- list(); pat <- list()

for (i in seq_len(N)) {
  b0 <- rnorm(1, 0, TH$tau_b0); b1 <- rnorm(1, 0, TH$tau_b1)

  ## visit times in MONTHS: lognormal gaps with median ~6, follow-up 18-135 mo
  horizon <- runif(1, 18, 135)
  t <- cumsum(rlnorm(40, log(6), 0.5))
  t <- c(round(runif(1, 0.5, 3), 3), t)
  t <- sort(unique(round(t[t <= horizon], 3)))
  if (length(t) < 2) t <- c(1, 7)
  ell <- log1p(t / 12)                       # model works on years

  ## ---- longitudinal ----
  bm <- rbinom(length(t), 1, 0.9)
  mu <- TH$beta0 + (TH$beta_time + b1) * ell + TH$beta_time2 * ell^2 +
        b0 + TH$beta_bm * bm
  y_lat <- rnorm(length(t), mu, TH$sigma_y)
  y     <- pmax(y_lat, FLOOR)                # reported at the floor
  ratio <- 10^y
  abl   <- round(runif(length(t), 2e4, 1e5))
  bcr   <- round(ratio * abl)
  gap   <- c(t[1], diff(t))

  long[[i]] <- data.frame(
    patient_id = sprintf("S%04d", i), patient_num = i,
    visit_index = seq_along(t), t_months = t, gap_months = gap,
    sample_type = ifelse(bm == 1, "BM", "PB"), sample_bm = bm,
    bcr_abl_copy = bcr, abl_copy = abl, ratio = ratio, IS = ratio,
    log_bcrabl = ifelse(bcr > 0, log10(bcr), 0),
    log_abl = log10(abl), log_IS = y, log_mrd = y,
    cmr = as.integer(y <= FLOOR + 1e-9), dmr = as.integer(y <= CD),
    dmr_from_log_mrd = as.integer(y <= CD), stringsAsFactors = FALSE)

  ## ---- at-risk intervals and the event, from the interval hazard ----
  ts <- head(t, -1); te <- tail(t, -1)
  midl <- log1p(0.5 * (ts + te) / 12); gapl <- log1p((te - ts) / 12)
  m_mid <- TH$beta0 + (TH$beta_time + b1) * midl + TH$beta_time2 * midl^2 + b0
  eta <- TH$gamma0 + TH$gamma_time * midl + TH$gamma_gap * gapl + TH$alpha_mrd * m_mid
  p <- 1 - exp(-exp(eta) * pmax((te - ts) / 12, 1e-6))
  ev <- rbinom(length(p), 1, p); first <- which(ev == 1)[1]
  keep <- if (is.na(first)) seq_along(ts) else seq_len(first)
  e <- rep(0L, length(keep)); if (!is.na(first)) e[first] <- 1L

  iv[[i]] <- data.frame(
    patient_id = sprintf("S%04d", i), patient_num = i,
    visit_index = keep, t_start = ts[keep], t_end = te[keep],
    gap_months = te[keep] - ts[keep], log_gap = log1p(te[keep] - ts[keep]),
    event_interval = e, stringsAsFactors = FALSE)

  ## ---- patient level ----
  delta <- as.integer(!is.na(first))
  tt <- if (delta == 1) te[first] else max(t)
  ph <- round(pmax(0, 100 * exp(-runif(1, 0.3, 1.2) * seq_len(15))), 1)
  pat[[i]] <- data.frame(
    patient_id = sprintf("S%04d", i), patient_num = i, dmr_event = delta,
    time_to_dmr_or_censor = tt, followup_months = max(t),
    n_visits = length(t), n_intervals = length(keep),
    baseline_log_mrd = y[1], last_log_mrd = y[length(y)], min_log_mrd = min(y),
    ever_cmr = as.integer(any(y <= FLOOR + 1e-9)),
    age = round(rnorm(1, 45, 14)), duration_years = round(runif(1, 0.5, 15), 1),
    ph_baseline_pct = 100,
    ph_3m_pct = ph[1], ph_6m_pct = ph[2], ph_9m_pct = ph[3], ph_12m_pct = ph[4],
    ph_18m_pct = ph[5], ph_2y_pct = ph[6], ph_3y_pct = ph[7], ph_4y_pct = ph[8],
    ph_5y_pct = ph[9], ph_6y_pct = ph[10], ph_7y_pct = ph[11], ph_8y_pct = ph[12],
    ph_9y_pct = ph[13], ph_10y_pct = ph[14], ph_11y_pct = ph[15], ph_12y_pct = 0,
    sex_male = rbinom(1, 1, 0.55), has_complete_covariates = 1L,
    stringsAsFactors = FALSE)
}

L <- do.call(rbind, long); I <- do.call(rbind, iv); P <- do.call(rbind, pat)

## ---- verify the schema matches the real files exactly -----------------------
## The real files are absent from the public archive by design, so this check is
## skipped there and runs only in the authors' working copy.
check <- function(sim, real_file, label) {
  rf <- file.path(outdir, real_file)
  if (!file.exists(rf)) {
    message(sprintf("  %-22s %d columns (real file absent; schema not verified)",
                    label, length(names(sim))))
    return(sim)
  }
  want <- names(utils::read.csv(rf, nrows = 1, check.names = FALSE))
  miss <- setdiff(want, names(sim)); extra <- setdiff(names(sim), want)
  message(sprintf("  %-22s columns %d/%d%s%s", label, length(names(sim)), length(want),
                  if (length(miss)) paste0(" MISSING:", paste(miss, collapse = ",")) else "",
                  if (length(extra)) paste0(" EXTRA:", paste(extra, collapse = ",")) else ""))
  if (length(miss)) stop("Schema mismatch in ", label, "; the pipeline would fail.")
  sim[, want, drop = FALSE]
}
message("Schema check against the real files:")
L <- check(L, "real_longitudinal_analysis.csv",      "longitudinal")
I <- check(I, "real_interval_survival_analysis.csv", "interval")
P <- check(P, "real_patient_level_analysis.csv",     "patient")

write.csv(L, file.path(outdir, "synthetic_longitudinal_analysis.csv"),      row.names = FALSE)
write.csv(I, file.path(outdir, "synthetic_interval_survival_analysis.csv"), row.names = FALSE)
write.csv(P, file.path(outdir, "synthetic_patient_level_analysis.csv"),     row.names = FALSE)

message(sprintf("\nWrote synthetic cohort: %d patients, %d observations, %d intervals, DMR rate %.2f",
                nrow(P), nrow(L), nrow(I), mean(P$dmr_event)))
message("Floor rate ", round(mean(L$log_mrd <= -5 + 1e-9), 3),
        " | median gap ", round(median(L$gap_months[L$gap_months > 0]), 2), " months")
message("\nUse prefix = \"synthetic\" in 03_build_stan_data.R / 04_prepare_joint_model_data.R.")
