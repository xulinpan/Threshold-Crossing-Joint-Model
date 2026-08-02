## 03_numeric_results.R ----------------------------------------------------
## Deliverable 2: numeric results from the fitted Bayesian models.
## Posterior summaries, HMC diagnostics, assay-floor PPC, patient-level DMR
## calibration + Brier score, and dynamic-prediction probabilities.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd()) })
source(file.path(this_dir, "00_setup.R"))
source(file.path(this_dir, "01_simulate_data.R"))

coh     <- readRDS(file.path(DIR_OUT, "primary_cohort.rds"))
res_thr <- readRDS(file.path(DIR_OUT, "fit_threshold.rds"))
res_int <- readRDS(file.path(DIR_OUT, "fit_interval.rds"))

## ---- draws + diagnostics are materialized by 02_fit_models.R -------------
## (Saved objects hold a posterior::draws_df, not cmdstanr temp-CSV pointers.)
get_draws    <- function(res) res$draws
diag_summary <- function(res) res$diag

## ---- 1. posterior summaries ----------------------------------------------
key_thr <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_thr")
key_int <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1",
             "gamma0","gamma1","gamma2","alpha")

## Note: "q5"/"q95" are output column names, not summary functions; produce
## them with quantile2() (whose default probs are 0.05 and 0.95).
sumfun <- function(dr) posterior::summarise_draws(
  dr, "mean", "sd",
  ~posterior::quantile2(.x, probs = c(0.05, 0.95)),
  "rhat", "ess_bulk", "ess_tail")
sm_thr <- sumfun(get_draws(res_thr))
sm_int <- sumfun(get_draws(res_int))
sm_thr <- subset(sm_thr, variable %in% key_thr)
sm_int <- subset(sm_int, variable %in% key_int)
sm_thr$truth <- unlist(TRUTH[sm_thr$variable])
sm_int$truth <- unlist(TRUTH[sm_int$variable])
write.csv(sm_thr, file.path(DIR_OUT, "posterior_summary_threshold.csv"), row.names = FALSE)
write.csv(sm_int, file.path(DIR_OUT, "posterior_summary_interval.csv"), row.names = FALSE)

## ---- 2. diagnostics ------------------------------------------------------
diag <- rbind(cbind(model = "threshold", diag_summary(res_thr)),
              cbind(model = "interval",  diag_summary(res_int)))
diag$max_rhat  <- c(max(sm_thr$rhat, na.rm = TRUE), max(sm_int$rhat, na.rm = TRUE))
diag$min_ess   <- c(min(sm_thr$ess_bulk, na.rm = TRUE), min(sm_int$ess_bulk, na.rm = TRUE))
write.csv(diag, file.path(DIR_OUT, "hmc_diagnostics.csv"), row.names = FALSE)

## ---- 3. assay-floor posterior predictive check ---------------------------
## Plain data.frame for column indexing (avoids draws_df metadata warnings).
dr <- as.data.frame(get_draws(res_thr))
pfloor_cols <- grep("^p_floor\\[", names(dr), value = TRUE)
post_floor  <- mean(colMeans(dr[pfloor_cols]))
obs_floor   <- mean(coh$long$floor)
## non-floor predictive coverage (90/95) from y_rep
yrep_cols <- grep("^y_rep\\[", names(dr), value = TRUE)
nonfloor  <- which(coh$long$floor == 0)
cover <- sapply(c(0.90, 0.95), function(lev) {
  a <- (1 - lev) / 2
  lo <- apply(dr[yrep_cols[nonfloor]], 2, quantile, a)
  hi <- apply(dr[yrep_cols[nonfloor]], 2, quantile, 1 - a)
  mean(coh$long$y[nonfloor] >= lo & coh$long$y[nonfloor] <= hi)
})
ppc <- data.frame(obs_floor_rate = obs_floor, post_mean_floor_prob = post_floor,
                  cover90 = cover[1], cover95 = cover[2])
write.csv(ppc, file.path(DIR_OUT, "ppc_summary.csv"), row.names = FALSE)

## ---- 4. patient-level DMR calibration + Brier ----------------------------
pdmr_cols <- grep("^pDMR_end\\[", names(dr), value = TRUE)
pred <- colMeans(dr[pdmr_cols])
obs  <- coh$pat$delta
brier <- mean((pred - obs)^2)
cal <- glm(obs ~ qlogis(pmin(pmax(pred, 1e-4), 1 - 1e-4)), family = binomial())
calib <- data.frame(model = "threshold", brier = brier,
                    cal_intercept = coef(cal)[1], cal_slope = coef(cal)[2],
                    mean_pred = mean(pred), obs_rate = mean(obs))
write.csv(calib, file.path(DIR_OUT, "calibration_summary.csv"), row.names = FALSE)
saveRDS(list(pred = pred, obs = obs), file.path(DIR_OUT, "dmr_pred_obs.rds"))

## ---- 5. dynamic prediction pi_i(s, H) ------------------------------------
## running minimum via soft-min on a grid, using posterior-mean parameters
b0m <- colMeans(dr[grep("^b0\\[", names(dr))])
b1m <- colMeans(dr[grep("^b1\\[", names(dr))])
pm  <- sapply(key_thr, function(v) mean(dr[[v]]))
Mrun <- function(t, u0, u1, G = GRID, kappa = KAPPA) {
  gg <- ell(seq(0, t, length.out = G))
  mv <- pm["beta0"] + pm["beta1"] * gg + pm["beta2"] * gg^2 + u0 + u1 * gg
  -log(sum(exp(-kappa * mv))) / kappa
}
land <- c(0.5, 1.0, 1.5, 2.0); hor <- c(0.5, 1.0, 2.0)
rows <- list()
for (s in land) for (H in hor) {
  atrisk <- which(coh$pat$true_onset > s)     # patients DMR-free at landmark (oracle risk set)
  if (length(atrisk) < 5) next
  pi_i <- sapply(atrisk, function(i) {
    Ms  <- Mrun(s, b0m[i], b1m[i]); MsH <- Mrun(s + H, b0m[i], b1m[i])
    1 - pnorm((MsH - C_D) / pm["sigma_thr"]) / max(pnorm((Ms - C_D) / pm["sigma_thr"]), 1e-6)
  })
  realized <- as.integer(coh$pat$true_onset[atrisk] <= s + H)
  rows[[length(rows) + 1]] <- data.frame(
    landmark = s, horizon = H, n_at_risk = length(atrisk),
    mean_pred = mean(pi_i), observed = mean(realized),
    brier = mean((pi_i - realized)^2))
}
dynpred <- do.call(rbind, rows)
write.csv(dynpred, file.path(DIR_OUT, "dynamic_prediction.csv"), row.names = FALSE)

cat("\n==== POSTERIOR SUMMARY (threshold-crossing) ====\n"); print(sm_thr, digits = 3)
cat("\n==== DIAGNOSTICS ====\n"); print(diag, digits = 3)
cat("\n==== PPC ====\n"); print(ppc, digits = 3)
cat("\n==== CALIBRATION / BRIER ====\n"); print(calib, digits = 3)
cat("\n==== DYNAMIC PREDICTION ====\n"); print(dynpred, digits = 3)
message("Numeric results written to ", DIR_OUT)
