## =====================================================================
## 14_fit_serial_correlation.R
##
## Audit item A2/B1: refit the joint model with a within-patient OU/AR(1)
## serial residual so that sigma_y is pure measurement error. Builds the
## serial-correlation bookkeeping from the real stan data, reorders all
## observation-level vectors by (patient, time), and fits
## glw_joint_interval_dmr_serial.stan.
##
## STATUS: NEW MODEL — validate on simulated data first. Compare the
## recovered sigma_y (expect it to drop toward assay CV ~0.2-0.3 on log10),
## tau_b1, and alpha_mrd against the primary c_F=-5.0 fit.
##
## Requires: cmdstanr with a working toolchain.
## =====================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("03_Data", "04_Code")))) &&
         dirname(d) != d) d <- dirname(d)
  d
}
root <- find_root()

stan_file <- file.path(root, "04_Code/Stan/glw_joint_interval_dmr_serial.stan")
rds_file  <- file.path(root, "03_Data/Processed/stan_data_real_joint_interval_dmr.rds")
out_dir   <- file.path(root, "08_Model"); dir.create(out_dir, showWarnings = FALSE)

sd0 <- readRDS(rds_file)
req <- c("N_obs","id_obs","y","is_floor","t_obs","sample_bm","floor_value")
stopifnot(all(req %in% names(sd0)))

## ---- reorder observation-level vectors by (patient, time) ------------
ord <- order(sd0$id_obs, sd0$t_obs)
obs_vecs <- c("id_obs","y","is_floor","t_obs","sample_bm")
sd <- sd0
for (v in obs_vecs) sd[[v]] <- sd0[[v]][ord]

## ---- build first_in_patient, prev_obs, dt_prev -----------------------
N <- sd$N_obs
first_in_patient <- integer(N)
prev_obs <- integer(N)
dt_prev  <- numeric(N)
for (n in seq_len(N)) {
  if (n == 1 || sd$id_obs[n] != sd$id_obs[n - 1]) {
    first_in_patient[n] <- 1L
    prev_obs[n] <- 0L
    dt_prev[n]  <- 0
  } else {
    first_in_patient[n] <- 0L
    prev_obs[n] <- n - 1L            # valid because sorted -> prev index < n
    dt_prev[n]  <- sd$t_obs[n] - sd$t_obs[n - 1]
  }
}
stopifnot(all(prev_obs < seq_len(N)))         # required by the Stan model
sd$first_in_patient <- first_in_patient
sd$prev_obs <- prev_obs
sd$dt_prev  <- dt_prev

## ---- fit -------------------------------------------------------------
mod <- cmdstan_model(stan_file)
fit <- mod$sample(
  data = sd, chains = 4, parallel_chains = 4,
  iter_warmup = 1500, iter_sampling = 1500,
  adapt_delta = 0.995, max_treedepth = 12,
  seed = 20260726, refresh = 200
)

pars <- c("beta0","beta_time","beta_time2","beta_bm","sigma_y",
          "sigma_u","ell","sigma_total","gamma0","gamma_time",
          "gamma_gap","alpha_mrd","tau_b[1]","tau_b[2]")
summ <- fit$summary(pars)
write.csv(summ, file.path(out_dir, "serial_correlation_summary.csv"),
          row.names = FALSE)
saveRDS(fit, file.path(out_dir, "stan_joint_interval_dmr_serial_fit.rds"))
print(fit$diagnostic_summary())

cat("\n=== serial-correlation model: key parameters ===\n")
print(summ[summ$variable %in%
  c("sigma_y","sigma_u","ell","tau_b[2]","alpha_mrd","beta_time"),
  c("variable","mean","q5","q95","rhat")])
cat("\nImplied assay CV (log10 -> natural): exp(log(10)*sigma_y_mean) not needed;\n",
    "report sigma_y on log10 scale and compare with lab analytical CV.\n")
message("\nDone. If sigma_y drops toward ~0.2-0.3 and tau_b1/alpha are stable,",
        " this supports the serial-correlation refinement (Discussion).")
