## run_all.R ---------------------------------------------------------------
## Orchestrates the threshold-crossing pipeline. From this folder:
##     Rscript run_all.R
## Requires cmdstanr (preferred) or rstan, plus posterior and RColorBrewer.
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else getwd() })
step <- function(f) { message("\n=== Running ", f, " ==="); source(file.path(this_dir, f)) }

step("00_setup.R")
step("02_fit_models.R")        # deliverable 1: Stan Bayesian fit (threshold + interval)
step("03_numeric_results.R")   # deliverable 2: posterior summaries, diagnostics, calibration
step("04_simulation_study.R")  # deliverable 3: simulation study across scenarios
step("05_figures.R")           # deliverable 4: 600-dpi RColorBrewer figures

message("\nAll steps complete. Outputs in ./outputs (csv), ./outputs/figures (600 dpi PNG).")
