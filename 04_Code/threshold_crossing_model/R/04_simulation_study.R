## =============================================================================
## DEPRECATED 2026-07-29 — DO NOT REPORT RESULTS FROM THIS SCRIPT.
##
## This study estimates by cmdstanr $optimize(), i.e. the mode of the JOINT
## posterior over the hyperparameters AND the latent random effects z0/z1. The
## joint mode is not a consistent estimator of variance components, and the
## output shows exactly that pathology: the bias in tau1 GROWS with sample size
## (+9.01 at n=87, +11.70 at n=150, +16.11 at n=300 against a truth of 2.203).
## sigma_y is biased low by ~0.3 at every n and mean predicted risk is 0.67
## against an observed rate of 0.79. None of outputs/simulation_summary.csv or
## simulation_raw.csv is reportable.
##
## USE 06_simulation_redesign.R INSTEAD (full HMC, correctly specified central
## cell, quantile-based coverage). This file is kept only because fit_one() is
## still referenced elsewhere; it should not be re-run for the manuscript.
## =============================================================================

## 04_simulation_study.R ---------------------------------------------------
## Deliverable 3: simulate the model and evaluate operating characteristics.
## Repeated simulation across scenarios (assay-floor handling, sample size,
## informative monitoring). Uses cmdstanr $optimize (penalized MLE / posterior
## mode) per replicate for speed; set USE_MCMC = TRUE for full coverage.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd()) })
source(file.path(this_dir, "00_setup.R"))
source(file.path(this_dir, "01_simulate_data.R"))

N_REP   <- 100                       # replicates per scenario (reduce for a quick run)
USE_MCMC <- FALSE                    # TRUE -> short MCMC (slow) enabling coverage
stopifnot(requireNamespace("cmdstanr", quietly = TRUE))
mod_int <- cmdstanr::cmdstan_model(file.path(DIR_STAN, "interval_hazard_joint.stan"))

## predicted P(documented DMR by end) for interval model from a parameter vector
predict_interval <- function(coh, par, b0, b1) {
  sapply(seq_len(coh$n), function(i) {
    ti <- coh$long$t[coh$long$pid == i]; yi <- coh$long$y[coh$long$pid == i]
    edges <- c(0, ti); dv <- which(yi <= C_D)
    kmax <- if (length(dv) > 0) dv[1] else length(ti)
    surv <- 1
    for (k in seq_len(kmax)) {
      L <- edges[k]; R <- edges[k + 1]; mid <- log1p((L + R) / 2)
      m <- par["beta0"] + par["beta1"] * mid + par["beta2"] * mid^2 + b0[i] + b1[i] * mid
      logh <- par["gamma0"] + par["gamma1"] * mid + par["gamma2"] * log1p(R - L) + par["alpha"] * m
      surv <- surv * exp(-exp(logh) * (R - L))
    }
    1 - surv
  })
}

## generic init to avoid log(0) rejections at random starting values
init_int <- function(N) function() list(
  beta0 = -2.0, beta1 = -2.5, beta2 = 0.3, beta_bm = 0.0,
  sigma_y = 1.5, tau0 = 0.7, tau1 = 1.2,
  gamma0 = -2.0, gamma1 = 0.0, gamma2 = 0.0, alpha = -0.3,
  z0 = rep(0, N), z1 = rep(0, N))

fit_one <- function(coh, exact = FALSE) {
  sd <- build_stan_interval(coh)
  if (exact) sd$floor_ind <- rep(0L, sd$Nobs)      # treat floor values as exact
  ini <- init_int(sd$N)
  if (USE_MCMC) {
    f <- mod_int$sample(data = sd, chains = 2, iter_warmup = 400, iter_sampling = 400,
                        adapt_delta = 0.95, refresh = 0, show_messages = FALSE, init = ini)
    dd <- posterior::as_draws_df(f$draws())
    est <- sapply(c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1",
                    "gamma0","gamma1","gamma2","alpha"), function(v) mean(dd[[v]]))
    se  <- sapply(names(est), function(v) sd(dd[[v]]))
    b0  <- colMeans(dd[grep("^b0\\[", names(dd))]); b1 <- colMeans(dd[grep("^b1\\[", names(dd))])
  } else {
    f <- mod_int$optimize(data = sd, jacobian = TRUE, refresh = 0, seed = 1, init = ini)
    s <- f$summary()
    getv <- function(v) s$estimate[s$variable == v]
    est <- sapply(c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1",
                    "gamma0","gamma1","gamma2","alpha"), getv)
    se  <- rep(NA_real_, length(est)); names(se) <- names(est)
    b0  <- sapply(seq_len(coh$n), function(i) getv(sprintf("b0[%d]", i)))
    b1  <- sapply(seq_len(coh$n), function(i) getv(sprintf("b1[%d]", i)))
  }
  pred <- predict_interval(coh, est, b0, b1); obs <- coh$pat$delta
  list(est = est, se = se, brier = mean((pred - obs)^2),
       mean_pred = mean(pred), obs_rate = mean(obs))
}

scen <- list(
  list(name = "ni_n87",  n = 87,  inf = FALSE),
  list(name = "ni_n150", n = 150, inf = FALSE),
  list(name = "ni_n300", n = 300, inf = FALSE),
  list(name = "inf_n150",n = 150, inf = TRUE)
)
ests <- c("floor", "exact")

results <- list()
for (sc in scen) {
  for (est_kind in ests) {
    if (sc$name == "inf_n150" && est_kind == "exact") next
    for (r in seq_len(N_REP)) {
      coh <- simulate_cohort(sc$n, informative = sc$inf, seed = 1000 * match(sc$name, sapply(scen, `[[`, "name")) + r)
      fit <- tryCatch(fit_one(coh, exact = (est_kind == "exact")), error = function(e) NULL)
      if (is.null(fit)) next
      results[[length(results) + 1]] <- data.frame(
        scen = sc$name, n = sc$n, informative = sc$inf, estimator = est_kind, rep = r,
        beta1 = fit$est["beta1"], sigma_y = fit$est["sigma_y"], tau1 = fit$est["tau1"],
        alpha = fit$est["alpha"], brier = fit$brier,
        mean_pred = fit$mean_pred, obs_rate = fit$obs_rate,
        se_beta1 = fit$se["beta1"], se_tau1 = fit$se["tau1"])
    }
    message(sprintf("done %s / %s", sc$name, est_kind))
  }
}
res <- do.call(rbind, results)
saveRDS(res, file.path(DIR_OUT, "simulation_raw.rds"))
write.csv(res, file.path(DIR_OUT, "simulation_raw.csv"), row.names = FALSE)

## ---- summary: bias / RMSE by scenario x estimator ----
agg <- do.call(rbind, lapply(split(res, list(res$scen, res$estimator), drop = TRUE), function(d) {
  data.frame(scen = d$scen[1], n = d$n[1], estimator = d$estimator[1], nrep = nrow(d),
             bias_beta1 = mean(d$beta1) - TRUTH$beta1, rmse_beta1 = sqrt(mean((d$beta1 - TRUTH$beta1)^2)),
             bias_tau1  = mean(d$tau1)  - TRUTH$tau1,  rmse_tau1  = sqrt(mean((d$tau1  - TRUTH$tau1)^2)),
             bias_sigma = mean(d$sigma_y) - TRUTH$sigma_y,
             brier = mean(d$brier), mean_pred = mean(d$mean_pred), obs_rate = mean(d$obs_rate))
}))
write.csv(agg, file.path(DIR_OUT, "simulation_summary.csv"), row.names = FALSE)
cat("\n==== SIMULATION SUMMARY ====\n"); print(agg, digits = 3)
message("Simulation study complete -> ", DIR_OUT)
