## =====================================================================
## 13_sensitivity_censoring_point.R
##
## Audit item A1/B2: the primary analysis left-censors at c_F = -5.0, but
## because no retained observation falls in (-5.0, -4.5], deep values were
## almost certainly collapsed at MR4.5, in which case the correct censoring
## point is c_F = -4.5. This script re-fits the primary (independent) joint
## model at c_F = -4.5 and compares it with the c_F = -5.0 fit.
##
## Why only `floor_value` changes: the set of floor observations (is_floor)
## is identical under both censoring points, because every floor observation
## sits at -5.0 and there are no observations in (-5.0, -4.5]. The left-
## censored contribution Phi((c_F - mu)/sigma_y) is the only term affected.
##
## Requires: cmdstanr with a working toolchain (the repo already ships
## compiled .exe models, so cmdstan is configured).
## Run from the project root or anywhere; paths are resolved relative to it.
## =====================================================================

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
})

## ---- locate project root (folder containing 03_Data, 04_Code) --------
find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("03_Data", "04_Code")))) &&
         dirname(d) != d) d <- dirname(d)
  d
}
root <- find_root()
message("Project root: ", root)

stan_file <- file.path(root, "04_Code/Stan/glw_joint_interval_dmr_independent.stan")
rds_file  <- file.path(root, "03_Data/Processed/stan_data_real_joint_interval_dmr.rds")
out_dir   <- file.path(root, "08_Model")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- load real stan data ---------------------------------------------
stan_data <- readRDS(rds_file)
stopifnot(!is.null(stan_data$floor_value), !is.null(stan_data$is_floor))
message("Original floor_value = ", stan_data$floor_value,
        " ; n floor obs = ", sum(stan_data$is_floor),
        " / ", stan_data$N_obs)

## sanity: confirm no observations lie in (-5.0, -4.5] among NON-floor obs,
## which would invalidate the "is_floor unchanged" argument.
nonfloor_in_gap <- with(stan_data,
  sum(is_floor == 0 & y > -5.0 & y <= -4.5))
message("Non-floor observations in (-5.0, -4.5]: ", nonfloor_in_gap,
        "  (expected 0)")

## ---- refit at c_F = -4.5 ---------------------------------------------
stan_data_45 <- stan_data
stan_data_45$floor_value <- -4.5

mod <- cmdstan_model(stan_file)

fit_45 <- mod$sample(
  data            = stan_data_45,
  chains          = 4,
  parallel_chains = 4,
  iter_warmup     = 1000,
  iter_sampling   = 1000,
  adapt_delta     = 0.99,
  max_treedepth   = 12,
  seed            = 20260726,
  refresh         = 200
)

## ---- save summaries and diagnostics ----------------------------------
pars <- c("beta0","beta_time","beta_time2","beta_bm","sigma_y",
          "gamma0","gamma_time","gamma_gap","alpha_mrd","tau_b[1]","tau_b[2]")

summ_45 <- fit_45$summary(pars)
write.csv(summ_45,
  file.path(out_dir, "sensitivity_cF_minus4p5_summary.csv"),
  row.names = FALSE)

diag_45 <- fit_45$diagnostic_summary()
saveRDS(fit_45, file.path(out_dir, "stan_joint_interval_dmr_cF_minus4p5_fit.rds"))
writeLines(
  c(sprintf("num_divergent: %s", paste(diag_45$num_divergent, collapse=",")),
    sprintf("num_max_treedepth: %s", paste(diag_45$num_max_treedepth, collapse=",")),
    sprintf("ebfmi: %s", paste(round(diag_45$ebfmi,3), collapse=","))),
  file.path(out_dir, "sensitivity_cF_minus4p5_diagnostics.txt"))

## ---- side-by-side comparison with the c_F=-5.0 primary fit ------------
primary <- tryCatch(
  read.csv(file.path(out_dir,
    "stan_joint_interval_dmr_independent_summary.csv")),
  error = function(e) NULL)

cat("\n=== c_F = -4.5 posterior means (95% CrI) ===\n")
print(summ_45[, c("variable","mean","q5","q95")])
if (!is.null(primary)) {
  cat("\nCompare with the c_F = -5.0 primary summary in\n",
      file.path(out_dir,
        "stan_joint_interval_dmr_independent_summary.csv"), "\n")
}

message("\nDone. Wire sensitivity_cF_minus4p5_summary.csv into the ",
        "Prior/sensitivity discussion of the manuscript (Section 4).")
