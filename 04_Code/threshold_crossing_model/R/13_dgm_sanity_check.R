## =============================================================================
## 13_dgm_sanity_check.R  —  does the simulated mechanism look like the cohort?
##
## An ADEMP simulation study is only persuasive if its data-generating mechanism
## resembles the application. This script prints the handful of marginal
## summaries a reviewer will check, side by side with the analysed cohort:
##   visits per patient, visit gap, follow-up, assay-floor rate, DMR rate.
##
## It exists because the pre-2026-07-29 mechanism used Exp(rate 4) visit times to
## t = 5, i.e. a mean gap of 3 months, against an observed median gap of 6
## months -- every simulated dataset was roughly twice as informative as the real
## one. Run this after changing any DGM switch, and BEFORE launching a long run.
##
## RUN (from the R/ folder):
##   Rscript 13_dgm_sanity_check.R
##   $env:DGM_VERSION='v1_legacy'; Rscript 13_dgm_sanity_check.R   # compare
## =============================================================================

.chk_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
source(file.path(.chk_dir, "dgm_common.R"))
dgm_describe()

N_CHECK <- as.integer(Sys.getenv("N_CHECK", "300"))
R_CHECK <- as.integer(Sys.getenv("R_CHECK", "20"))
set.seed(20260729)

CENTRAL <- list(trajectory = "quadratic", floor_kind = "fixed",
                err_family = "gaussian", re_family = "gaussian", informative = FALSE)

## ---- summaries of one dataset ----------------------------------------------
describe_ds <- function(long, intervals = NULL, floor_col = "floor_ind") {
  v   <- as.numeric(table(long$pid))
  gap <- unlist(lapply(split(long$t, long$pid), function(x) diff(sort(x))))
  fu  <- vapply(split(long$t, long$pid), max, numeric(1))
  fl  <- long[[floor_col]]
  dmr <- if (!is.null(intervals))
           mean(vapply(split(intervals$event, intervals$pid), function(e) as.numeric(any(e == 1)), numeric(1)))
         else NA_real_
  c(visits_med = median(v), visits_mean = mean(v),
    gap_med_mo = 12 * median(gap), gap_p_gt6mo = mean(gap > 0.5),
    fu_med_y = median(fu), fu_max_y = max(fu),
    floor_rate = mean(fl), dmr_rate = dmr,
    y_med = median(long$y))
}

## ---- simulated ------------------------------------------------------------
sim <- do.call(rbind, lapply(seq_len(R_CHECK), function(r) {
  ds <- simulate_dataset(N_CHECK, CENTRAL)
  describe_ds(ds$long, ds$intervals, "floor_ind")
}))
sim_mean <- colMeans(sim)

## ---- REAL cohort summaries -------------------------------------------------
## Reference marginals come from the actual CML cohort in 03_Data/Processed,
## NOT from outputs/primary_cohort.rds. The latter is synthetic --
## 02_fit_models.R:67 builds it as simulate_cohort(n = 200, seed = 101) -- so
## comparing against it only shows that this generator reproduces another
## generator, which is not evidence that the study resembles the application.
##
## Only marginal SUMMARIES are read here (visit counts, gaps, follow-up, floor
## and DMR rates). No patient-level data enters the simulation, and nothing in
## this folder is fitted to the cohort; the sandbox stays simulation-only.
## Real files are in months; this folder works in years.
REAL_DIR <- Sys.getenv("REAL_DATA_DIR",
                       file.path(DGM_DIR, "..", "..", "..", "03_Data", "Processed"))
read_real <- function() {
  fl <- file.path(REAL_DIR, "real_longitudinal_analysis.csv")
  fp <- file.path(REAL_DIR, "real_patient_level_analysis.csv")
  if (!file.exists(fl) || !file.exists(fp)) return(NULL)
  L <- utils::read.csv(fl, stringsAsFactors = FALSE)
  P <- utils::read.csv(fp, stringsAsFactors = FALSE)
  gaps <- L$gap_months[is.finite(L$gap_months) & L$gap_months > 0]
  fu   <- P$followup_months / 12
  c(visits_med  = stats::median(P$n_visits),
    visits_mean = mean(P$n_visits),
    gap_med_mo  = stats::median(gaps),
    gap_p_gt6mo = mean(gaps > 6),
    fu_med_y    = stats::median(fu),
    fu_max_y    = max(fu),
    floor_rate  = mean(L$log_mrd <= -5 + 1e-9, na.rm = TRUE),
    dmr_rate    = mean(P$dmr_event, na.rm = TRUE),
    y_med       = stats::median(L$log_mrd, na.rm = TRUE))
}
obs <- read_real()
if (is.null(obs))
  message("NOTE: real cohort summaries not found under ", REAL_DIR,
          " -- set REAL_DATA_DIR. Printing simulated summaries only.")

## ---- report ---------------------------------------------------------------
lab <- c(visits_med="visits/patient (median)", visits_mean="visits/patient (mean)",
         gap_med_mo="visit gap, months (median)", gap_p_gt6mo="P(gap > 6 months)",
         fu_med_y="follow-up, years (median)", fu_max_y="follow-up, years (max)",
         floor_rate="assay-floor rate", dmr_rate="documented DMR rate",
         y_med="observed y (median)")
cat(sprintf("\n%-28s %12s %12s %10s\n", "quantity", "simulated", "cohort", "ratio"))
cat(strrep("-", 66), "\n")
for (k in names(lab)) {
  s <- sim_mean[[k]]; o <- if (is.null(obs)) NA_real_ else obs[[k]]
  rr <- if (is.na(o) || o == 0) NA_real_ else s / o
  cat(sprintf("%-28s %12.3f %12s %10s\n", lab[[k]], s,
              if (is.na(o)) "-" else sprintf("%.3f", o),
              if (is.na(rr)) "-" else sprintf("%.2f", rr)))
}
cat(strrep("-", 66), "\n")

## Leave an artefact: a console-only check cannot be audited later, and this is
## the check that validates the largest change to the mechanism.
chk <- data.frame(dgm_tag = dgm_tag(), quantity = names(lab), label = unname(lab),
                  simulated = round(as.numeric(sim_mean[names(lab)]), 4),
                  cohort = if (is.null(obs)) NA_real_ else round(as.numeric(obs[names(lab)]), 4),
                  row.names = NULL)
chk$ratio <- round(chk$simulated / chk$cohort, 3)
chk$n_check <- N_CHECK; chk$r_check <- R_CHECK
chk$checked_at <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
write.csv(chk, file.path(DGM_OUT, paste0("dgm_sanity_check_", dgm_tag(), ".csv")),
          row.names = FALSE)
cat("Wrote outputs/dgm_sanity_check_", dgm_tag(), ".csv\n", sep = "")
cat("Simulated column averages", R_CHECK, "datasets of n =", N_CHECK, "\n")
cat("A ratio far from 1.00 for the visit or follow-up rows means the study is\n")
cat("not evaluating the estimator under the conditions of the application.\n\n")

if (!is.null(obs)) {
  bad <- names(lab)[c(1,3,5)][abs(log(sim_mean[c(1,3,5)] / obs[c(1,3,5)])) > log(1.5)]
  if (length(bad))
    warning("Simulated and observed differ by >1.5x on: ",
            paste(lab[bad], collapse = "; "), call. = FALSE, immediate. = TRUE)
}
