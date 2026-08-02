## =====================================================================
## 16_fit_unified_model.R  —  ROUTE B fitting ladder
##
## Fits glw_unified_joint.stan, in which the three models of the paper are
## data-flag options of ONE program sharing one latent trajectory and one
## set of random effects:
##
##   M1  link=hazard,   no spline            -> joint longitudinal-interval
##   M2  link=crossing, no spline            -> latent threshold-crossing
##   M3  link=crossing, +durability +relapse -> multi-state
##   M3s M3 + spline                          -> multi-state, flexible traj.
##
## Fitting them as a LADDER (simple -> complex, each initialised from the
## previous) is deliberate: the full model asks sigma_thr, sigma_zeta,
## tau_b0, tau_b1 and sigma_y to separate simultaneously on 87 patients
## with ~48% floor-censored observations and few relapses. The ladder makes
## any identifiability failure visible at the step that causes it.
##
## Units: all times are converted to YEARS (source CSVs are in months).
## =====================================================================

suppressPackageStartupMessages({library(cmdstanr); library(posterior)})

find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("03_Data", "04_Code")))) &&
         dirname(d) != d) d <- dirname(d)
  d
}
root <- find_root()
proc <- file.path(root, "03_Data/Processed")
out  <- file.path(root, "08_Model"); dir.create(out, showWarnings = FALSE)

C_D        <- -4.5      # DMR threshold
FLOOR      <- -5.0      # set to -4.5 for the censoring-point sensitivity
KNOT_YEARS <- 2
DUR_WINDOW <- 2         # durability window w, years

## ---- read source data ------------------------------------------------
long <- read.csv(file.path(proc, "real_longitudinal_analysis.csv"))
intv <- read.csv(file.path(proc, "real_interval_survival_analysis.csv"))
pat  <- read.csv(file.path(proc, "real_patient_level_analysis.csv"))

long <- long[order(long$patient_num, long$t_months), ]
N_pat <- length(unique(long$patient_num))
stopifnot(all(sort(unique(long$patient_num)) == seq_len(N_pat)))

## ---- onset intervals per patient (crossing form) ---------------------
## event patients: the interval whose end carries the first documented DMR
## censored patients: (0, last t_end]
ons <- do.call(rbind, lapply(split(intv, intv$patient_num), function(d) {
  d <- d[order(d$t_start), ]
  ev <- which(d$event_interval == 1)
  if (length(ev)) {
    data.frame(id = d$patient_num[1], L = d$t_start[ev[1]],
               R = d$t_end[ev[1]], event = 1L)
  } else {
    data.frame(id = d$patient_num[1], L = 0,
               R = max(d$t_end), event = 0L)
  }
}))
stopifnot(nrow(ons) == N_pat)

## ---- durability + relapse rows (multi-state), from observed data -----
## durability: among patients with documented onset, did every measurement
## in (T_D, T_D + w] stay at/below c_D?
## relapse (confirmed): two consecutive post-onset measurements above c_D.
onset_time <- setNames(ons$R, ons$id)
dur <- rel <- NULL
for (i in ons$id[ons$event == 1]) {
  td <- onset_time[[as.character(i)]]
  li <- long[long$patient_num == i & long$t_months > td, ]
  li <- li[order(li$t_months), ]
  # durability window
  w  <- li[li$t_months <= td + DUR_WINDOW * 12, ]
  if (nrow(w))
    dur <- rbind(dur, data.frame(id = i, lo = td/12,
                                 hi = (td + DUR_WINDOW*12)/12,
                                 ok = as.integer(all(w$log_mrd <= C_D))))
  # confirmed relapse: first pair of consecutive values above c_D
  above <- li$log_mrd > C_D
  k <- if (length(above) >= 2) which(above[-length(above)] & above[-1])[1] else NA
  if (!is.na(k)) {
    rel <- rbind(rel, data.frame(id = i, L = td/12,
                                 R = li$t_months[k + 1]/12, event = 1L))
  } else if (nrow(li)) {
    rel <- rbind(rel, data.frame(id = i, L = td/12,
                                 R = max(li$t_months)/12, event = 0L))
  }
}
message(sprintf("durability rows: %d (durable %d) | relapse rows: %d (events %d)",
                nrow(dur), sum(dur$ok), nrow(rel), sum(rel$event)))

## ---- assemble the unified data list ----------------------------------
mk_data <- function(link, spline, with_multistate) {
  list(
    link_type = link, use_spline = spline,
    knot_years = KNOT_YEARS, c_D = C_D, floor_value = FLOOR,
    N_obs = nrow(long), N_pat = N_pat,
    id_obs = as.integer(long$patient_num),
    y = as.numeric(long$log_mrd),
    is_floor = as.integer(long$log_mrd <= FLOOR),
    t_obs = as.numeric(long$t_months) / 12,
    sample_bm = as.numeric(long$sample_bm),
    N_int = nrow(intv), id_int = as.integer(intv$patient_num),
    t_start = intv$t_start/12, t_end = intv$t_end/12,
    gap = intv$gap_months/12,
    event_interval = as.integer(intv$event_interval),
    N_ons = nrow(ons), id_ons = as.integer(ons$id),
    ons_L = ons$L/12, ons_R = ons$R/12, ons_event = as.integer(ons$event),
    N_dur = if (with_multistate) nrow(dur) else 0L,
    id_dur = if (with_multistate) as.integer(dur$id) else integer(0),
    dur_lo = if (with_multistate) dur$lo else numeric(0),
    dur_hi = if (with_multistate) dur$hi else numeric(0),
    dur_ok = if (with_multistate) as.integer(dur$ok) else integer(0),
    N_rel = if (with_multistate) nrow(rel) else 0L,
    id_rel = if (with_multistate) as.integer(rel$id) else integer(0),
    rel_L = if (with_multistate) rel$L else numeric(0),
    rel_R = if (with_multistate) rel$R else numeric(0),
    rel_event = if (with_multistate) as.integer(rel$event) else integer(0)
  )
}

mod <- cmdstan_model(file.path(root, "04_Code/Stan/glw_unified_joint.stan"))

## ---- explicit initial values -----------------------------------------
## Stan's default init draws sigma_y as low as ~0.14, and with 239
## left-censored observations normal_lcdf(floor | mu, sigma_y) then
## underflows to -Inf ("Log probability evaluates to log(0)"). Starting
## near the known posterior avoids that entirely.
mk_init <- function(link, spline) function() {
  z <- list(
    beta0      = rnorm(1, -2.0, 0.2),
    beta_time  = rnorm(1, -3.5, 0.2),
    beta_time2 = abs(rnorm(1, 0.5, 0.1)),
    beta_bm    = rnorm(1,  0.6, 0.2),
    sigma_y    = abs(rnorm(1, 1.8, 0.1)),
    tau_b      = abs(rnorm(2, c(1.0, 2.2), 0.1)),
    z_b        = matrix(rnorm(2 * N_pat, 0, 0.1), nrow = 2)
  )
  if (link == 0) {
    z$gamma0     <- as.array(rnorm(1, -4.0, 0.3))
    z$gamma_time <- as.array(rnorm(1, -0.2, 0.1))
    z$gamma_gap  <- as.array(rnorm(1, -1.8, 0.2))
    z$alpha_mrd  <- as.array(rnorm(1, -1.3, 0.1))
  } else {
    z$sigma_thr  <- as.array(abs(rnorm(1, 0.3, 0.05)))
  }
  if (spline == 1) {
    z$sigma_zeta <- as.array(abs(rnorm(1, 0.5, 0.1)))
    z$z_zeta     <- rnorm(N_pat, 0, 0.1)
  }
  z
}

fit_one <- function(tag, link, spline, ms, adapt = 0.99, init = NULL) {
  message("\n=== fitting ", tag, " ===")
  if (is.null(init)) init <- mk_init(link, spline)
  f <- mod$sample(data = mk_data(link, spline, ms),
                  chains = 4, parallel_chains = 4,
                  iter_warmup = 1500, iter_sampling = 1500,
                  adapt_delta = adapt, max_treedepth = 12,
                  seed = 20260726, refresh = 500, init = init)
  ## NOTE: select against the fit's OWN variable list. An earlier version
  ## filtered through a "[1]"-only helper, which silently dropped tau_b[2]
  ## (the random-slope SD) from every summary.
  want <- c("beta0","beta_time","beta_time2","beta_bm","sigma_y",
            "tau_b[1]","tau_b[2]","gamma0[1]","gamma_time[1]","gamma_gap[1]",
            "alpha_mrd[1]","sigma_thr[1]","sigma_zeta[1]",
            "nadir_years","nadir_depth")
  s <- f$summary(intersect(want, f$summary()$variable))
  write.csv(s, file.path(out, paste0("unified_", tag, "_summary.csv")),
            row.names = FALSE)
  print(f$diagnostic_summary())
  print(s[, c("variable","mean","q5","q95","rhat","ess_bulk")])
  f
}

## ---- the ladder ------------------------------------------------------
m1  <- fit_one("M1_hazard",            link = 0, spline = 0, ms = FALSE)
m2  <- fit_one("M2_crossing",          link = 1, spline = 0, ms = FALSE)
m3  <- fit_one("M3_multistate",        link = 1, spline = 0, ms = TRUE, adapt = 0.995)
m3s <- fit_one("M3s_multistate_spline",link = 1, spline = 1, ms = TRUE, adapt = 0.995)

## ---- LOO on the longitudinal block -----------------------------------
## Valid for comparing spline vs no spline (same link, same data block).
## NOT like-for-like across link_type, whose event blocks differ in
## dimension (275 intervals vs 87 patients) -- do not rank M1 vs M2 by LOO.
if (requireNamespace("loo", quietly = TRUE)) {
  l3  <- loo::loo(m3$draws("log_lik"))
  l3s <- loo::loo(m3s$draws("log_lik"))
  print(loo::loo_compare(list(M3 = l3, M3s = l3s)))
  capture.output(loo::loo_compare(list(M3 = l3, M3s = l3s)),
                 file = file.path(out, "unified_loo_spline_vs_none.txt"))
}

message("\nDone. Summaries in 08_Model/unified_*_summary.csv")
message("Check: does sigma_thr stay stable from M2 -> M3 -> M3s? ",
        "If it collapses when the spline enters, the spline is still ",
        "absorbing trajectory misspecification.")
