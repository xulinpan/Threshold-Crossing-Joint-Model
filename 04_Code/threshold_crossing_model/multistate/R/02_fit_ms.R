## 02_fit_ms.R -------------------------------------------------------------
## Fit the bidirectional multi-state threshold-crossing model in Stan and
## report parameter recovery plus onset / durability / relapse calibration.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else getwd() })
source(file.path(this_dir, "01_simulate_ms.R"))
suppressMessages({ library(cmdstanr); library(posterior) })

OUT <- file.path(this_dir, "..", "outputs"); dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
STAN <- file.path(this_dir, "..", "stan", "multistate_threshold.stan")

## ---- simulate a demonstration cohort ----
coh <- simulate_ms(n = 500, seed = 101)
sd  <- build_stan_ms(coh)

## ---- sensible inits (random starts make crossing terms evaluate to log 0) ----
init_ms <- function() list(beta0 = -2.0, beta1 = -3.0, beta2 = 1.0, beta_bm = 0.0,
                           sigma_y = 1.0, tau0 = 0.6, tau1 = 0.8, sigma_thr = 0.3,
                           z0 = rep(0, sd$N), z1 = rep(0, sd$N))

mod <- cmdstan_model(STAN)
fit <- mod$sample(data = sd, chains = 4, parallel_chains = 4,
                  iter_warmup = 1000, iter_sampling = 1000,
                  adapt_delta = 0.99, max_treedepth = 12,
                  init = init_ms, seed = 20260709, refresh = 200)

dr <- as_draws_df(fit$draws())
saveRDS(list(draws = dr, coh = coh), file.path(OUT, "fit_multistate.rds"))

## ---- parameter recovery ----
keys <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_thr")
sm <- summarise_draws(dr, "mean","sd",
        ~quantile2(.x, probs = c(0.05, 0.95)), "rhat","ess_bulk")
sm <- subset(sm, variable %in% keys)
sm$truth <- unlist(TRUTH_MS[sm$variable])
write.csv(sm, file.path(OUT, "ms_posterior_summary.csv"), row.names = FALSE)

## ---- diagnostics ----
ds <- fit$diagnostic_summary()
diag <- data.frame(divergences = sum(ds$num_divergent),
                   max_treedepth = sum(ds$num_max_treedepth),
                   ebfmi_min = min(ds$ebfmi),
                   max_rhat = max(sm$rhat), min_ess = min(sm$ess_bulk))
write.csv(diag, file.path(OUT, "ms_diagnostics.csv"), row.names = FALSE)

## ---- calibration of the three transition probabilities ----
drdf <- as.data.frame(dr)
grab <- function(pre) colMeans(drdf[grep(sprintf("^%s\\[", pre), names(drdf))])
calib <- function(pred, obs, lab) {
  b <- glm(obs ~ qlogis(pmin(pmax(pred, 1e-4), 1 - 1e-4)), family = binomial())
  data.frame(target = lab, brier = mean((pred - obs)^2),
             cal_intercept = coef(b)[1], cal_slope = coef(b)[2],
             mean_pred = mean(pred), obs_rate = mean(obs))
}
P <- coh$pat
cal <- rbind(
  calib(grab("p_onset_byC"),   P$onset,                          "onset by C"),
  calib(grab("p_sustained_W"), as.integer(P$onset == 1 & P$sustained == 1), "sustained DMR (W)"),
  calib(grab("p_relapse_byC"), as.integer(P$onset == 1 & P$relapse == 1),   "relapse by C"))
write.csv(cal, file.path(OUT, "ms_transition_calibration.csv"), row.names = FALSE)

cat("\n==== PARAMETER RECOVERY ====\n"); print(sm[, c("variable","mean","q5","q95","rhat","truth")], digits = 3)
cat("\n==== DIAGNOSTICS ====\n"); print(diag, digits = 3)
cat("\n==== TRANSITION CALIBRATION ====\n"); print(cal, digits = 3)
message("Multi-state fit complete -> ", OUT)
