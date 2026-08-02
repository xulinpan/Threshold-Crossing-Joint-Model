## 05_figures.R ------------------------------------------------------------
## Deliverable 4: high-quality figures (res = 600 dpi) using RColorBrewer.
## Reads outputs from 01-04; each panel is skipped gracefully if its input
## is missing so the script can be run after a partial pipeline.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd()) })
source(file.path(this_dir, "00_setup.R"))
source(file.path(this_dir, "01_simulate_data.R"))

fig <- function(name) file.path(DIR_FIG, name)

## ============================ FIGURE 1 ====================================
## Simulated latent log-MRD trajectories with population trend and thresholds.
coh <- if (file.exists(file.path(DIR_OUT, "primary_cohort.rds")))
  readRDS(file.path(DIR_OUT, "primary_cohort.rds")) else simulate_cohort(120, seed = 3)

open_png(fig("figR_01_trajectories.png"), width = 7, height = 5)
op <- par(mar = c(4.2, 4.2, 1.5, 1))
ids <- unique(coh$long$pid)
xr <- range(coh$long$t); yr <- range(coh$long$y) + c(-0.2, 0.2)
plot(NA, xlim = xr, ylim = yr, xlab = "Years since treatment start",
     ylab = "log-MRD", las = 1)
grays <- adjustcolor(PAL_BLUE[5], 0.35)
for (id in ids) {
  d <- coh$long[coh$long$pid == id, ]
  lines(d$t, d$y, col = grays, lwd = 0.7)
}
## smoothed population trend
lo <- loess(y ~ t, data = coh$long, span = 0.5)
tt <- seq(xr[1], xr[2], length.out = 100)
lines(tt, predict(lo, newdata = data.frame(t = tt)), col = COL$trend, lwd = 3)
abline(h = C_D, col = COL$dmr, lwd = 2, lty = 2)
abline(h = C_F, col = COL$cmr, lwd = 2, lty = 3)
legend("topright", bty = "n", cex = 0.85,
       legend = c("patient trajectory", "smoothed trend",
                  "DMR threshold (-4.5)", "CMR / floor (-5.0)"),
       col = c(grays, COL$trend, COL$dmr, COL$cmr),
       lwd = c(0.7, 3, 2, 2), lty = c(1, 1, 2, 3))
par(op); dev.off()

## ============================ FIGURE 2 ====================================
## Assay-floor posterior predictive check (needs fitted threshold model).
ppc_file <- file.path(DIR_OUT, "ppc_summary.csv")
if (file.exists(ppc_file) && file.exists(file.path(DIR_OUT, "fit_threshold.rds"))) {
  res_thr <- readRDS(file.path(DIR_OUT, "fit_threshold.rds"))
  dr <- as.data.frame(res_thr$draws)   # plain df for column indexing (no draws_df warnings)
  mu_cols <- grep("^mu_obs\\[", names(dr), value = TRUE)
  mu_hat  <- colMeans(dr[mu_cols])
  pf      <- colMeans(dr[grep("^p_floor\\[", names(dr), value = TRUE)])
  open_png(fig("figR_02_ppc_floor.png"), width = 8, height = 4)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2, 1))
  ## (a) observed vs posterior-mean fitted
  plot(mu_hat, coh$long$y, pch = 19, cex = 0.4,
       col = ifelse(coh$long$floor == 1, COL$cens, COL$obs),
       xlab = "Posterior mean fitted", ylab = "Observed log-MRD", main = "(a) Fitted vs observed", las = 1)
  abline(0, 1, col = "grey40", lty = 2)
  legend("topleft", bty = "n", pch = 19, col = c(COL$obs, COL$cens),
         legend = c("non-floor", "floor"), cex = 0.85)
  ## (b) floor probability vs time
  ord <- order(coh$long$t)
  plot(coh$long$t[ord], pf[ord], type = "n", ylim = c(0, 1),
       xlab = "Years since treatment start", ylab = "P(log-MRD <= floor)",
       main = "(b) Assay-floor probability", las = 1)
  points(coh$long$t, pf, pch = 19, cex = 0.4,
         col = ifelse(coh$long$floor == 1, COL$cens, COL$obs))
  par(op); dev.off()
}

## ============================ FIGURE 3 ====================================
## Patient-level DMR calibration (needs 03 output).
po_file <- file.path(DIR_OUT, "dmr_pred_obs.rds")
if (file.exists(po_file)) {
  po <- readRDS(po_file); pred <- po$pred; obs <- po$obs
  brks <- quantile(pred, seq(0, 1, length.out = 6)); brks[1] <- -Inf; brks[length(brks)] <- Inf
  bin <- cut(pred, brks); mp <- tapply(pred, bin, mean); op_ <- tapply(obs, bin, mean)
  nb <- tapply(pred, bin, length)
  open_png(fig("figR_03_dmr_calibration.png"), width = 5.2, height = 5)
  op <- par(mar = c(4.2, 4.2, 1.5, 1))
  plot(mp, op_, xlim = 0:1, ylim = 0:1, pch = 19, cex = 1 + 2 * nb / max(nb),
       col = COL$thresh, xlab = "Mean predicted DMR probability",
       ylab = "Observed DMR rate", las = 1)
  abline(0, 1, col = "grey40", lty = 2)
  lines(mp, op_, col = COL$thresh, lwd = 2)
  par(op); dev.off()
}

## ============================ FIGURE 4 ====================================
## Dynamic prediction: predicted vs observed by landmark/horizon.
dp_file <- file.path(DIR_OUT, "dynamic_prediction.csv")
if (file.exists(dp_file)) {
  dp <- read.csv(dp_file); hs <- sort(unique(dp$horizon))
  cols <- setNames(PAL_DARK[seq_along(hs)], hs)
  open_png(fig("figR_04_dynamic_prediction.png"), width = 6, height = 5)
  op <- par(mar = c(4.2, 4.2, 1.5, 1))
  plot(NA, xlim = range(dp$landmark), ylim = c(0, 1), xlab = "Landmark time (years)",
       ylab = "P(DMR within horizon)", las = 1)
  for (H in hs) {
    d <- dp[dp$horizon == H, ]
    lines(d$landmark, d$mean_pred, col = cols[as.character(H)], lwd = 2)
    points(d$landmark, d$mean_pred, col = cols[as.character(H)], pch = 19)
    points(d$landmark, d$observed, col = cols[as.character(H)], pch = 1, cex = 1.3)
  }
  legend("topright", bty = "n", lwd = 2, col = cols,
         legend = paste0("H = ", hs, " yr"), title = "predicted (line/solid)")
  legend("bottomright", bty = "n", pch = 1, legend = "observed (open)", cex = 0.85)
  par(op); dev.off()
}

## ============================ FIGURE 5 ====================================
## Simulation: floor vs exact bias for tau1 and beta1 across n.
sim_file <- file.path(DIR_OUT, "simulation_raw.rds")
if (file.exists(sim_file)) {
  res <- readRDS(sim_file)
  ns <- c(87, 150, 300); scn <- c("ni_n87", "ni_n150", "ni_n300")
  open_png(fig("figR_05_sim_floor_bias.png"), width = 8.5, height = 4)
  op <- par(mfrow = c(1, 2), mar = c(4.2, 4.4, 2, 1))
  for (v in c("tau1", "beta1")) {
    truev <- TRUTH[[v]]
    plot(NA, xlim = c(70, 320), ylim = if (v == "tau1") c(0.5, 2.6) else c(-4, -2),
         xlab = "sample size n", ylab = bquote(hat(.(as.name(v)))),
         main = if (v == "tau1") expression(tau[b1]) else expression(beta[time]), las = 1)
    for (est in c("floor", "exact")) {
      mus <- ses <- numeric(length(ns))
      for (j in seq_along(scn)) {
        d <- res[res$scen == scn[j] & res$estimator == est, v]
        mus[j] <- mean(d); ses[j] <- sd(d) / sqrt(length(d))
      }
      cc <- if (est == "floor") COL$floor else COL$exact
      lines(ns, mus, col = cc, lwd = 2, type = "b", pch = if (est == "floor") 19 else 17)
      arrows(ns, mus - 1.96 * ses, ns, mus + 1.96 * ses, angle = 90, code = 3, length = 0.03, col = cc)
    }
    abline(h = truev, lty = 2, col = "grey30")
    if (v == "tau1") legend("bottomright", bty = "n", col = c(COL$floor, COL$exact),
                            pch = c(19, 17), lwd = 2, legend = c("floor left-censored", "floor as exact"))
  }
  par(op); dev.off()

  ## ========================== FIGURE 6 ==================================
  ## Simulation calibration-in-the-large (mean predicted vs observed).
  agg <- aggregate(cbind(mean_pred, obs_rate) ~ scen + estimator, data = res, mean)
  open_png(fig("figR_06_sim_calibration.png"), width = 5.5, height = 5)
  op <- par(mar = c(4.2, 4.2, 1.5, 1))
  plot(agg$mean_pred, agg$obs_rate, xlim = c(0.3, 0.95), ylim = c(0.3, 0.95),
       pch = ifelse(agg$estimator == "floor", 19, 17),
       col = ifelse(agg$estimator == "floor", COL$floor, COL$exact),
       xlab = "Mean predicted DMR probability", ylab = "Observed DMR rate", las = 1, cex = 1.4)
  abline(0, 1, lty = 2, col = "grey40")
  text(agg$mean_pred, agg$obs_rate, labels = agg$scen, pos = 3, cex = 0.6)
  legend("topleft", bty = "n", pch = c(19, 17), col = c(COL$floor, COL$exact),
         legend = c("floor left-censored", "floor as exact"))
  par(op); dev.off()
}

message("Figures (600 dpi) written to ", DIR_FIG)
