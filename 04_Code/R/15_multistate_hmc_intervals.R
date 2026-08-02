## =====================================================================
## 15_multistate_hmc_intervals.R
##
## Audit item B4: the multi-state (threshold-crossing) table in the
## manuscript reports marginal-mode estimates with ASYMPTOTIC intervals
## that are far narrower than Hamiltonian Monte Carlo intervals on the same
## data. A full HMC threshold fit already exists under
##   04_Code/threshold_crossing_model/outputs/fit_threshold.rds
## This script extracts proper 95% HMC credible intervals to replace the
## asymptotic ones, and audits convergence.
##
## IMPORTANT CAVEAT: the existing threshold fit recorded 219 divergent
## transitions (04_Code/threshold_crossing_model/outputs/hmc_diagnostics.csv).
## Divergences bias the posterior geometry; the intervals below should be
## treated as provisional and the model RE-RUN at higher adapt_delta
## (>= 0.995) with reparameterisation before the intervals are quoted.
## =====================================================================

suppressPackageStartupMessages({
  library(posterior)
})

find_root <- function(start = getwd()) {
  d <- normalizePath(start, winslash = "/", mustWork = FALSE)
  while (!all(dir.exists(file.path(d, c("03_Data", "04_Code")))) &&
         dirname(d) != d) d <- dirname(d)
  d
}
root <- find_root()
tc   <- file.path(root, "04_Code/threshold_crossing_model/outputs")
out  <- file.path(root, "08_Model"); dir.create(out, showWarnings = FALSE)

## ---- convergence audit -----------------------------------------------
diag_path <- file.path(tc, "hmc_diagnostics.csv")
if (file.exists(diag_path)) {
  diag <- read.csv(diag_path)
  cat("=== HMC diagnostics (existing threshold fit) ===\n"); print(diag)
  bad <- diag[diag$model == "threshold", , drop = FALSE]
  if (nrow(bad) && bad$num_divergent[1] > 0)
    warning(sprintf(
      "Threshold fit has %d divergent transitions; re-run at higher adapt_delta before quoting intervals.",
      bad$num_divergent[1]))
}

## ---- extract 95% HMC intervals from the fitted draws -----------------
fit_path <- file.path(tc, "fit_threshold.rds")
if (!file.exists(fit_path)) stop("fit_threshold.rds not found: run the ",
  "threshold_crossing_model pipeline (04_Code/threshold_crossing_model/R/run_all.R) first.")

## The pipeline (02_fit_models.R) saves a LIST: res$draws is a
## posterior::draws_df and res$diag holds the diagnostics. Handle that,
## and also tolerate a bare draws object or a cmdstanr fit.
res <- readRDS(fit_path)
draws <- if (is.list(res) && !is.null(res[["draws"]]) &&
             !is.function(res[["draws"]])) {
  posterior::as_draws_df(res[["draws"]])
} else if (inherits(res, "draws")) {
  posterior::as_draws_df(res)
} else if (!is.null(res$draws) && is.function(res$draws)) {
  posterior::as_draws_df(res$draws())
} else {
  stop("Unrecognised fit_threshold.rds structure; expected res$draws.")
}

pars <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_thr")
pars <- intersect(pars, posterior::variables(draws))

## Match the pipeline's summarise pattern (03_numeric_results.R): explicit
## measures + quantile2 for the interval columns (q2.5/q97.5 for 95%).
## Pass measure names as strings so the output columns are reliably named
## (mean, sd, q2.5, q97.5, rhat, ess_bulk, ess_tail).
tab <- posterior::summarise_draws(
  posterior::subset_draws(draws, variable = pars),
  "mean", "sd",
  ~posterior::quantile2(.x, probs = c(0.025, 0.975)),
  "rhat", "ess_bulk", "ess_tail")

write.csv(tab, file.path(out, "multistate_hmc_intervals_95.csv"),
          row.names = FALSE)

cat("\n=== Multi-state HMC 95% credible intervals ===\n")
show_cols <- intersect(c("variable","mean","q2.5","q97.5","rhat","ess_bulk"),
                       names(tab))
print(as.data.frame(tab)[, show_cols])

## ---- flag any parameter with a still-missing / degenerate interval ----
if ("truth" %in% names(read.csv(file.path(tc, "posterior_summary_threshold.csv"))))
  message("\nNOTE: posterior_summary_threshold.csv carries a 'truth' column, ",
          "indicating a SIMULATION-recovery run. Confirm fit_threshold.rds is ",
          "the REAL-cohort fit (primary_cohort.rds) before wiring into Table 9; ",
          "if not, refit on the real cohort.")

message("\nDone. Replace the asymptotic intervals in the multi-state table ",
        "with multistate_hmc_intervals_95.csv, and re-run at adapt_delta>=0.995 ",
        "to clear the 219 divergences before final quoting.")
