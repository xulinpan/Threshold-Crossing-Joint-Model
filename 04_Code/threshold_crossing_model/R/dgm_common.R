## =============================================================================
## dgm_common.R  —  SINGLE SOURCE OF TRUTH for every simulation study
##
## WHY THIS FILE EXISTS. The data-generating mechanism used to be duplicated in
## 06_simulation_redesign.R, 08_simulation_multistate.R, 10_sbc.R and
## 11_comparators.R. The four copies had already drifted apart (different TRUTH
## values, different visit processes, different trajectory departures), so the
## simulation study, the SBC check and the comparator study were no longer
## evaluating the same estimand. Everything generative now lives here; the
## driver scripts source this file and add only their own queue/reporting logic.
##
## THREE SWITCHES control the mechanism. All are read from the environment so a
## run is fully described by its env vars, and all are recorded in the output.
##
##   DGM_VERSION      v1_legacy | v2_cohort   (default v2_cohort)
##       v1_legacy  reproduces the mechanism behind the checkpoints written
##                  before 2026-07-29: TRUTH_LEGACY, and visits ~ Exp(rate 4) to
##                  t = 5, i.e. a mean gap of 3 months.
##       v2_cohort  the corrected mechanism. Visit gaps are lognormal with
##                  MEDIAN 0.5 years and follow-up is U(1.5, 7) years, matching
##                  the real cohort (median observed gap 6.0 months). v1's
##                  3-month spacing made every simulated dataset roughly twice
##                  as informative as the application, which flattered the
##                  operating characteristics and biased gamma2 in particular
##                  (gamma2 is the coefficient on log(1 + gap)).
##
##   TRUTH_SOURCE     legacy | posterior     (default follows DGM_VERSION)
##       legacy     the hard-coded values previously used in 06/11.
##       posterior  read from outputs/posterior_summary_interval.csv.
##
##       !! READ THIS BEFORE JUSTIFYING TRUTH IN THE MANUSCRIPT !!
##       "posterior" does NOT mean "fitted to the CML cohort". This folder is a
##       simulation-only methods sandbox: 02_fit_models.R:67 builds its primary
##       cohort as simulate_cohort(n = 200, seed = 101), so
##       posterior_summary_interval.csv is a SIMULATION-RECOVERY fit -- a
##       posterior for data generated from TRUTH in 00_setup.R. Using it as
##       TRUTH is therefore circular, and it must not be described in the paper
##       as "the posterior means of the fitted CML model" (00_setup.R:46 makes
##       that claim; it is wrong). The real cohort lives in
##       03_Data/Processed/real_*.csv and is analysed by 04_Code/R/, not here.
##
##       Both options are legitimate CHOICES of a parameter vector for a
##       simulation study. Neither is an empirical calibration. State the vector
##       and its provenance explicitly in the manuscript; do not claim it was
##       estimated from patient data. 13_dgm_sanity_check.R compares the
##       resulting marginals against the real cohort summaries, which is the
##       defensible way to argue the mechanism resembles the application.
##
##   EVENT_MECHANISM  hazard | threshold     (default hazard)
##       hazard     events drawn from the SAME complementary-log-log interval
##                  hazard that interval_hazard_joint.stan fits. The central
##                  cell is then correctly specified and gamma0..2, alpha have
##                  genuine true values, so parameter recovery is meaningful.
##       threshold  documented DMR is the deterministic event {first visit with
##                  observed y <= cD}, which is how DMR is ascertained in the
##                  real data. Under this mechanism gamma0..2 and alpha have NO
##                  generative counterpart and truth_for_scenario() returns NA
##                  for them, so no coverage is reported against a value that
##                  does not exist.
##
##       These two are mutually exclusive descriptions of the same endpoint and
##       the manuscript must commit to one. Keeping both here makes the choice
##       explicit and lets "threshold" be reported as a misspecification cell
##       for the hazard model rather than being silently assumed away.
##
## CHECKPOINT SAFETY. dgm_tag() returns a short string identifying the
## mechanism. Driver scripts MUST include it in the checkpoint directory so fits
## generated under different mechanisms can never be pooled by the aggregator.
## =============================================================================

## ---- 0. paths ---------------------------------------------------------------
dgm_script_dir <- function() {
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) return(dirname(normalizePath(sub("^--file=", "", m[1]))))
  of <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
  if (!is.null(of)) return(dirname(normalizePath(of)))
  getwd()
}
DGM_DIR  <- dgm_script_dir()
DGM_ROOT <- normalizePath(file.path(DGM_DIR, ".."), mustWork = FALSE)
DGM_OUT  <- file.path(DGM_ROOT, "outputs")

## ---- 1. switches ------------------------------------------------------------
DGM_VERSION <- Sys.getenv("DGM_VERSION", "v2_cohort")
if (!DGM_VERSION %in% c("v1_legacy", "v2_cohort"))
  stop("DGM_VERSION must be 'v1_legacy' or 'v2_cohort', got: ", DGM_VERSION)

EVENT_MECHANISM <- Sys.getenv("EVENT_MECHANISM", "hazard")
if (!EVENT_MECHANISM %in% c("hazard", "threshold"))
  stop("EVENT_MECHANISM must be 'hazard' or 'threshold', got: ", EVENT_MECHANISM)

TRUTH_SOURCE <- Sys.getenv("TRUTH_SOURCE",
                           if (DGM_VERSION == "v1_legacy") "legacy" else "posterior")
if (!TRUTH_SOURCE %in% c("legacy", "posterior"))
  stop("TRUTH_SOURCE must be 'legacy' or 'posterior', got: ", TRUTH_SOURCE)

## ---- 2. constants -----------------------------------------------------------
FLOOR <- -5.0     # assay floor c_F (reported value when an assay reads below it)
C_D   <- -4.5     # DMR threshold c_D

## ---- 3. TRUTH ---------------------------------------------------------------
PARAMS <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1",
            "gamma0","gamma1","gamma2","alpha")

TRUTH_LEGACY <- list(beta0=-2.07, beta1=-3.56, beta2=0.50, beta_bm=0.61,
                     sigma_y=1.81, tau0=0.97, tau1=2.20,
                     gamma0=-4.09, gamma1=-0.23, gamma2=-1.83, alpha=-1.27)

## Posterior means of the fitted interval model on the real cohort. Read from
## disk rather than hard-coded so TRUTH cannot silently go stale again.
truth_from_posterior <- function(
    file = file.path(DGM_OUT, "posterior_summary_interval.csv")) {
  if (!file.exists(file))
    stop("TRUTH_SOURCE='posterior' needs ", file,
         "\n  Run 03_numeric_results.R first, or set TRUTH_SOURCE=legacy.")
  s <- utils::read.csv(file, stringsAsFactors = FALSE)
  miss <- setdiff(PARAMS, s$variable)
  if (length(miss)) stop("posterior summary is missing: ", paste(miss, collapse=", "))
  out <- as.list(s$mean[match(PARAMS, s$variable)]); names(out) <- PARAMS
  out
}

TRUTH <- if (TRUTH_SOURCE == "legacy") TRUTH_LEGACY else truth_from_posterior()

## ---- 3b. prior sets ---------------------------------------------------------
## The priors used to be hard-coded in interval_hazard_joint.stan. The v2_cohort
## study showed undercoverage concentrated in the parameters whose values sit in
## the prior tails (gamma2 0.75-0.83, beta1 0.73, beta2 0.79, alpha 0.88), and
## 15_prior_sensitivity.R showed that doubling the scales restores essentially
## all of it. Priors are therefore a study FACTOR now, passed as data to
## interval_hazard_joint_pp.stan.
##
## "current" reproduces interval_hazard_joint.stan exactly. The others scale the
## regression and hazard SDs only; locations are untouched, so widening cannot be
## accused of nudging the prior toward the answer.
PRIOR_SET <- Sys.getenv("PRIOR_SET", "current")

.prior_base <- list(
  pr_beta0 = c(-2.5, 2), pr_beta1 = c(-1, 1), pr_beta2 = c(0, 0.5),
  pr_beta_bm = c(0, 1), pr_sigma_y = 1, pr_tau0 = 1, pr_tau1 = 1,
  pr_gamma0 = c(-2, 2), pr_gamma1 = c(0, 1), pr_gamma2 = c(0, 1),
  pr_alpha = c(-0.5, 0.75))
.prior_widen <- function(p, k) {
  for (nm in c("pr_beta1","pr_beta2","pr_gamma0","pr_gamma1","pr_gamma2","pr_alpha"))
    p[[nm]][2] <- p[[nm]][2] * k
  p
}
PRIOR_SETS <- list(
  current = .prior_base,
  wide2   = .prior_widen(.prior_base, 2),
  wide4   = .prior_widen(.prior_base, 4),
  vague   = { p <- .prior_widen(.prior_base, 10)
              p$pr_sigma_y <- 0.5; p$pr_tau0 <- 0.5; p$pr_tau1 <- 0.5; p })
if (!PRIOR_SET %in% names(PRIOR_SETS))
  stop("PRIOR_SET must be one of: ", paste(names(PRIOR_SETS), collapse = ", "),
       " -- got: ", PRIOR_SET)

## Stan data entries for the chosen prior set. gq = 0 suppresses the
## per-observation predictive quantities, which otherwise add ~4,650 columns per
## draw (~43 MB per chain) that a simulation study never reads.
stan_priors <- function(set = PRIOR_SET, gq = 0L, prior_only = 0L)
  c(PRIOR_SETS[[set]], list(gq = as.integer(gq), prior_only = as.integer(prior_only)))

## The prior set joins the tag so runs cannot pool. "current" is left off so the
## tag reproduces the pre-2026-07-30 value and completed checkpoints stay valid.
dgm_tag <- function() {
  base <- sprintf("%s_%s_%s", DGM_VERSION, TRUTH_SOURCE, EVENT_MECHANISM)
  if (identical(PRIOR_SET, "current")) base else paste0(base, "_", PRIOR_SET)
}

dgm_describe <- function() {
  message(sprintf("DGM: %s | truth=%s | events=%s", DGM_VERSION, TRUTH_SOURCE, EVENT_MECHANISM))
  message(sprintf("  visits: %s", if (DGM_VERSION == "v1_legacy")
    "Exp(rate 4) to t=5  [mean gap 3.0 months -- LEGACY, ~2x denser than cohort]"
    else "lognormal median 0.5 y, follow-up U(1.5, 7) y  [matches cohort]"))
  message(sprintf("  truth : beta1=%.3f sigma_y=%.3f tau1=%.3f alpha=%.3f",
                  TRUTH$beta1, TRUTH$sigma_y, TRUTH$tau1, TRUTH$alpha))
  invisible(NULL)
}

## ---- 4. random-effect / floor families --------------------------------------
## t3 draws are rescaled to the SAME standard deviation as the Gaussian case, so
## the departure is tail shape only and tau0/tau1/sigma_y keep their true values.
rre <- function(family, n, sd) {
  if (family == "gaussian") rnorm(n, 0, sd)
  else if (family == "t3")  sd * rt(n, 3) / sqrt(3)
  else stop("unknown re/err family: ", family)
}
floor_draw <- function(kind, n) {
  if (kind == "fixed") rep(FLOOR, n)
  else if (kind == "heterogeneous") FLOOR + rnorm(n, 0, 0.3)
  else stop("unknown floor kind: ", kind)
}

## ---- 5. latent trajectory ---------------------------------------------------
## ell = log(1 + t); b1 is a random slope on ell. "quadratic" is the form that
## interval_hazard_joint.stan assumes; the others are departures.
traj_fun <- function(kind, ell, b0, b1, th = TRUTH) {
  base <- th$beta0 + b0
  switch(kind,
    quadratic   = base + (th$beta1 + b1) * ell + th$beta2 * ell^2,
    monotone    = base + (th$beta1 + b1) * ell - 0.5 * ell,
    exponential = base + (th$beta1 + b1) * (1 - exp(-2 * ell)),
    stop("unknown trajectory: ", kind))
}

## ---- 6. visit process -------------------------------------------------------
## Returns visit times for ONE patient. Under informative monitoring the visit
## intensity depends on the patient's own latent trajectory (NHPP by thinning),
## which is what "informative" has to mean. The previous implementation in
## 06_simulation_redesign.R merely duplicated ~30% of visit times independently
## of the trajectory -- it created zero-length intervals and tested nothing.
VISIT_INF <- list(d0 = -0.7, d1 = 0.35, d2 = 0.2, wsd = 0.5)

simulate_visits_one <- function(b0, b1, trajectory = "quadratic",
                                informative = FALSE, th = TRUTH) {
  if (DGM_VERSION == "v1_legacy") { gap_med <- 0.25; horizon <- 5.0 }
  else                           { gap_med <- 0.50; horizon <- runif(1, 1.5, 7.0) }

  if (!informative) {
    if (DGM_VERSION == "v1_legacy") {
      t <- 0; v <- numeric(0)
      while (t < horizon) { t <- t + rexp(1, 1 / gap_med); if (t < horizon) v <- c(v, t) }
    } else {
      t <- 0; v <- numeric(0)
      repeat {
        t <- t + rlnorm(1, log(gap_med), 0.5)
        if (t > horizon) break
        v <- c(v, t)
      }
    }
  } else {
    ## NHPP thinning: intensity rises as the latent value approaches the
    ## threshold, so sicker/faster-responding patients are seen more often.
    p <- VISIT_INF; w <- rnorm(1, 0, p$wsd)
    lam_of <- function(tt) {
      m <- traj_fun(trajectory, log1p(tt), b0, b1, th)
      exp(p$d0 + p$d1 * (m + 4.5) + p$d2 * log1p(tt) + w)
    }
    grid    <- seq(1e-3, horizon, length.out = 50)
    lam_max <- max(lam_of(grid)) * 1.2 + 1e-6
    t <- 0; v <- numeric(0)
    repeat {
      t <- t + rexp(1, lam_max)
      if (t > horizon) break
      if (runif(1) < lam_of(t) / lam_max) v <- c(v, t)
    }
  }
  ## Guarantee an early baseline visit and at least two distinct times, so every
  ## patient contributes >= 1 at-risk interval of strictly positive length.
  ## v1_legacy pinned baseline at exactly 0.25 -- kept, so this branch really does
  ## reproduce the pre-2026-07-29 mechanism; v2 jitters it, which is what a real
  ## enrolment window looks like.
  baseline <- if (DGM_VERSION == "v1_legacy") 0.25 else runif(1, 0.02, 0.30)
  v <- sort(unique(round(c(baseline, v), 3)))
  if (length(v) < 2) v <- c(v, round(v[1] + gap_med, 3))
  v
}

## ---- 7. one dataset ---------------------------------------------------------
## scen is a one-row data.frame / list with: trajectory, floor_kind, err_family,
## re_family, informative (logical).
simulate_dataset <- function(n, scen, th = TRUTH) {
  informative <- isTRUE(as.logical(scen$informative))
  b0 <- rre(scen$re_family, n, th$tau0)
  b1 <- rre(scen$re_family, n, th$tau1)
  floors <- floor_draw(scen$floor_kind, n)

  long <- vector("list", n); intervals <- vector("list", n)
  for (i in seq_len(n)) {
    vt  <- simulate_visits_one(b0[i], b1[i], scen$trajectory, informative, th)
    ell <- log1p(vt)

    ## ---- longitudinal ----
    m  <- traj_fun(scen$trajectory, ell, b0[i], b1[i], th)
    bm <- rbinom(length(vt), 1, 0.9)            # 90% bone marrow / 10% peripheral
    y  <- m + th$beta_bm * bm + rre(scen$err_family, length(m), th$sigma_y)
    fl <- as.integer(y <= floors[i])
    y[fl == 1] <- FLOOR                          # lab REPORTS the nominal floor
    long[[i]] <- data.frame(pid = i, t = vt, ell = ell, y = y,
                            floor_ind = fl, bm = bm)

    ## ---- at-risk intervals between consecutive visits ----
    ts <- head(vt, -1); te <- tail(vt, -1)
    dl <- pmax(te - ts, 1e-6)

    if (EVENT_MECHANISM == "hazard") {
      midl  <- log1p(0.5 * (ts + te))
      gapl  <- log1p(te - ts)
      m_mid <- traj_fun(scen$trajectory, midl, b0[i], b1[i], th)
      logh  <- th$gamma0 + th$gamma1 * midl + th$gamma2 * gapl + th$alpha * m_mid
      ev    <- rbinom(length(dl), 1, 1 - exp(-exp(logh) * dl))
      first <- which(ev == 1)[1]
    } else {
      ## deterministic ascertainment: first visit whose REPORTED y is <= c_D.
      ## Visit k closes interval k, so a crossing at visit k is an event in
      ## interval k (k indexes te).
      dv    <- which(y[-1] <= C_D)
      first <- if (length(dv)) dv[1] else NA_integer_
    }
    keep <- if (is.na(first)) seq_along(ts) else seq_len(first)   # exit at first DMR
    e <- rep(0L, length(keep)); if (!is.na(first)) e[first] <- 1L
    intervals[[i]] <- data.frame(pid = i, t_start = ts[keep], t_end = te[keep],
                                 event = e)
  }
  list(long = do.call(rbind, long), intervals = do.call(rbind, intervals))
}

## ---- 8. Stan data -----------------------------------------------------------
make_stan_data <- function(ds) {
  L <- ds$long; I <- ds$intervals
  list(N = length(unique(L$pid)), Nobs = nrow(L), pid = L$pid, ell = L$ell, y = L$y,
       floor_ind = L$floor_ind, bm = L$bm, cF = FLOOR,
       Nint = nrow(I), pid_int = I$pid,
       midlog    = log1p(0.5 * (I$t_start + I$t_end)),
       gaplog    = log1p(I$t_end - I$t_start),
       delta_len = pmax(I$t_end - I$t_start, 1e-6),
       event     = I$event)
}

## ---- 9. scenario-specific truth --------------------------------------------
## CRITICAL. Under a departure, a working-model parameter may have no
## generative counterpart at all, in which case its "bias" and "coverage" are
## artefacts of comparing to a value that does not exist. Returning NA_real_
## makes the reporting layer skip it instead of publishing, e.g., "coverage
## 0.000 for beta2 under a monotone trajectory" -- which only says that a
## quadratic coefficient was absent from the truth.
##
## Parameters that survive a departure keep their generative value; those that
## are only PSEUDO-true (a KL projection of the working model, not a parameter
## of the DGM) are set to NA and listed in attr(, "undefined").
truth_for_scenario <- function(scen, th = TRUTH) {
  tr <- unlist(th[PARAMS]); names(tr) <- PARAMS
  undef <- character(0)

  ## the deterministic-ascertainment mechanism has no hazard parameters at all
  if (EVENT_MECHANISM == "threshold")
    undef <- c(undef, "gamma0", "gamma1", "gamma2", "alpha")

  if (identical(scen$trajectory, "monotone")) {
    ## generative fixed slope is (beta1 - 0.5); there is no quadratic term, and
    ## alpha multiplies a latent m of the wrong functional form
    tr["beta1"] <- th$beta1 - 0.5
    undef <- c(undef, "beta2", "alpha")
  } else if (identical(scen$trajectory, "exponential")) {
    ## b1 and beta1 scale a bounded function of ell, not ell, so neither the
    ## slope nor the random-slope SD is on a comparable scale
    undef <- c(undef, "beta1", "beta2", "tau1", "alpha")
  }

  undef <- unique(undef)
  tr[undef] <- NA_real_
  attr(tr, "undefined") <- undef
  tr
}

## ---- 10. per-parameter performance measures --------------------------------
## Skips parameters whose truth is NA (see truth_for_scenario).
metrics_one <- function(draws, param, truth) {
  if (is.na(truth))
    return(data.frame(param = param, post_mean = NA_real_, truth = NA_real_,
                      bias = NA_real_, covered = NA_integer_, ci_width = NA_real_,
                      truth_defined = FALSE))
  x  <- as.numeric(draws[[param]])
  qi <- quantile(x, c(0.025, 0.975))
  data.frame(param = param, post_mean = mean(x), truth = truth,
             bias = mean(x) - truth,
             covered = as.integer(truth >= qi[1] & truth <= qi[2]),
             ci_width = as.numeric(diff(qi)), truth_defined = TRUE)
}

## ---- 11. HMC fit + diagnostics ---------------------------------------------
## comp_fail used to be a single flag combining R-hat > 1.01, E-BFMI < 0.2 and
## any divergence. At 4 x 500 draws that lumped genuinely broken fits (R-hat up
## to 1.78 was present in the n=87 cell) together with entirely benign ones
## (R-hat 1.02 and no divergences, which flagged 100% of the n=300 cell). The
## two are now separate: fit_bad = exclude and report; fit_marginal = report only.
fit_hmc <- function(mod, standata, cfg) {
  fit <- mod$sample(data = standata, chains = cfg$chains,
                    parallel_chains = cfg$parallel_chains,
                    iter_warmup = cfg$iter_warmup, iter_sampling = cfg$iter_sampling,
                    adapt_delta = cfg$adapt_delta, max_treedepth = cfg$max_treedepth,
                    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  dr <- fit$draws(format = "draws_df")
  dg <- fit$diagnostic_summary(quiet = TRUE)
  ndiv  <- sum(dg$num_divergent)
  ebfmi <- suppressWarnings(min(dg$ebfmi))
  if (!is.finite(ebfmi)) ebfmi <- NA_real_
  ## Diagnostics are computed over the REPORTED parameters only. Taken over all
  ## quantities the model declares -- 11 parameters plus z0, z1, b0, b1, i.e.
  ## 4N+11, which is 611 at n=150 -- the maximum R-hat exceeds 1.01 almost surely
  ## by chance, and min ESS is driven by whichever latent effect happens to move
  ## least. The v2_cohort run showed the consequence: a "marginal" rate of
  ## 70-99% per cell and a reported min ESS of 6, neither of which says anything
  ## about the estimates being reported. The all-quantity versions are kept
  ## alongside, suffixed _all, so nothing is hidden.
  sm    <- posterior::summarise_draws(dr, "rhat", "ess_bulk", "ess_tail")
  keep  <- sm$variable %in% PARAMS
  mx <- function(x) suppressWarnings(max(x, na.rm = TRUE))
  mn <- function(x) suppressWarnings(min(x, na.rm = TRUE))
  rhat  <- mx(sm$rhat[keep]);     rhat_all <- mx(sm$rhat)
  ess   <- mn(sm$ess_bulk[keep]); ess_all  <- mn(sm$ess_bulk)
  ## Coverage is a property of the 2.5%/97.5% posterior quantiles, so TAIL ESS is
  ## the relevant precision measure for this study, not bulk ESS.
  esst  <- mn(sm$ess_tail[keep])
  list(draws = dr, ndiv = ndiv, ebfmi = ebfmi, rhat_max = rhat, ess_min = ess,
       ess_tail_min = esst, rhat_max_all = rhat_all, ess_min_all = ess_all,
       fit_bad      = as.integer(ndiv > 0 | (!is.na(ebfmi) & ebfmi < 0.2) | rhat > 1.05),
       fit_marginal = as.integer(rhat > 1.01 & rhat <= 1.05))
}

if (identical(environment(), globalenv()) && sys.nframe() == 0) dgm_describe()
