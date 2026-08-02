## 02_fit_models.R ---------------------------------------------------------
## Compile and fit the Stan models (deliverable 1: Bayesian model in Stan).
## cmdstanr is preferred; falls back to rstan if cmdstanr is unavailable.
## Fits both models to a single simulated "primary" cohort and saves draws.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd()) })
source(file.path(this_dir, "00_setup.R"))
source(file.path(this_dir, "01_simulate_data.R"))

## ---- sensible initial values ---------------------------------------------
## Random inits make the threshold-crossing / interval-hazard likelihoods
## evaluate to log(0) = -Inf ("Rejecting initial value"), so Stan cannot start.
## These generic, data-agnostic inits put every patient on a plausible declining
## trajectory with a non-degenerate crossing width, giving a finite initial
## log density. z0/z1 start at 0 (population trajectory); sigma_thr starts away
## from 0 to avoid Phi-difference underflow.
init_threshold <- function(N) function() list(
  beta0 = -2.0, beta1 = -2.5, beta2 = 0.3, beta_bm = 0.0,
  sigma_y = 1.5, tau0 = 0.7, tau1 = 1.2, sigma_thr = 0.4,
  z0 = rep(0, N), z1 = rep(0, N))
init_interval <- function(N) function() list(
  beta0 = -2.0, beta1 = -2.5, beta2 = 0.3, beta_bm = 0.0,
  sigma_y = 1.5, tau0 = 0.7, tau1 = 1.2,
  gamma0 = -2.0, gamma1 = 0.0, gamma2 = 0.0, alpha = -0.3,
  z0 = rep(0, N), z1 = rep(0, N))

## ---- backend abstraction -------------------------------------------------
fit_stan <- function(stan_file, data, init = NULL,
                     iter_warmup = 1000, iter_sampling = 1000,
                     chains = 4, adapt_delta = 0.95, max_treedepth = 12,
                     seed = 20260709) {
  if (requireNamespace("cmdstanr", quietly = TRUE)) {
    mod <- cmdstanr::cmdstan_model(stan_file)
    fit <- mod$sample(data = data, chains = chains, parallel_chains = chains,
                      iter_warmup = iter_warmup, iter_sampling = iter_sampling,
                      adapt_delta = adapt_delta, max_treedepth = max_treedepth,
                      seed = seed, refresh = 200,
                      init = if (is.null(init)) 0.5 else init)
    ## Materialize draws + diagnostics NOW so the saved object does not depend
    ## on cmdstanr's temporary CSV files (which vanish between R sessions).
    ds <- fit$diagnostic_summary()
    diag <- data.frame(num_divergent = sum(ds$num_divergent),
                       num_max_treedepth = sum(ds$num_max_treedepth),
                       ebfmi_min = min(ds$ebfmi))
    return(list(backend = "cmdstanr",
                draws = posterior::as_draws_df(fit$draws()), diag = diag))
  } else if (requireNamespace("rstan", quietly = TRUE)) {
    rstan::rstan_options(auto_write = TRUE)
    sm  <- rstan::stan_model(stan_file)
    fit <- rstan::sampling(sm, data = data, chains = chains,
                           warmup = iter_warmup, iter = iter_warmup + iter_sampling,
                           control = list(adapt_delta = adapt_delta,
                                          max_treedepth = max_treedepth),
                           seed = seed, refresh = 200,
                           init = if (is.null(init)) "random" else init)
    sp <- rstan::get_sampler_params(fit, inc_warmup = FALSE)
    diag <- data.frame(num_divergent = sum(sapply(sp, function(x) sum(x[, "divergent__"]))),
                       num_max_treedepth = NA_real_, ebfmi_min = NA_real_)
    return(list(backend = "rstan",
                draws = posterior::as_draws_df(fit), diag = diag))
  }
  stop("Neither cmdstanr nor rstan is installed. Install one to fit the model.")
}

## ---- primary simulated cohort (n = 200) ----------------------------------
coh <- simulate_cohort(n = 200, informative = FALSE, seed = 101)
saveRDS(coh, file.path(DIR_OUT, "primary_cohort.rds"))

## ---- fit threshold-crossing model ----------------------------------------
## NOTE. threshold_crossing_joint.stan is the SUPERSEDED grid formulation: it
## approximates the running minimum by a soft minimum over a G-point grid with
## sharpness kappa. The manuscript reports the CLOSED-FORM formulation instead
## (multistate/stan/multistate_spline.stan), which evaluates the running minimum
## at the interval endpoints and the clamped vertex. No number in the paper comes
## from the fit below.
##
## The diagnostics it writes to outputs/hmc_diagnostics.csv show 219 divergent
## transitions. That is expected and is the point: at kappa = 50 the soft minimum
## is effectively a hard minimum, so the log density is non-differentiable where
## the arg-minimum switches between grid points and HMC cannot integrate across
## it. 17_calibration_study.R quantifies this -- divergences in 100% of simulated
## fits, mean 128 of 2000 transitions, roughly 50x the cost of the interval
## model -- and the manuscript cites it as the motivation for the closed form.
## The fit is retained only so that comparison remains reproducible.
data_thr <- build_stan_threshold(coh)
res_thr  <- fit_stan(file.path(DIR_STAN, "threshold_crossing_joint.stan"), data_thr,
                     init = init_threshold(data_thr$N))
saveRDS(res_thr, file.path(DIR_OUT, "fit_threshold.rds"))

## ---- fit interval-hazard comparator --------------------------------------
data_int <- build_stan_interval(coh)
res_int  <- fit_stan(file.path(DIR_STAN, "interval_hazard_joint.stan"), data_int,
                     init = init_interval(data_int$N))
saveRDS(res_int, file.path(DIR_OUT, "fit_interval.rds"))

message("Fitting complete (backend: ", res_thr$backend, "). Draws saved to ", DIR_OUT)
