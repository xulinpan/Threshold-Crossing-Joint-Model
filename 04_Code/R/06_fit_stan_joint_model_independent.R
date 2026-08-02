options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !nzchar(x)) y else x
}

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) "") %||%
  "04_Code/R/06_fit_stan_joint_model_independent.R"
script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "06_fit_stan_joint_model.R"))

fit_stan_joint_model_independent <- function(
  chains = as_env_int("GLW_STAN_CHAINS", 4),
  parallel_chains = as_env_int("GLW_STAN_PARALLEL_CHAINS", min(chains, 4)),
  iter_warmup = as_env_int("GLW_STAN_WARMUP", 1000),
  iter_sampling = as_env_int("GLW_STAN_SAMPLING", 1000),
  seed = as_env_int("GLW_STAN_SEED", 20260709),
  model_prefix = Sys.getenv(
    "GLW_MODEL_PREFIX",
    unset = "stan_joint_interval_dmr_independent"
  )
) {
  project_root <- find_project_root()
  setwd(project_root)

  require_cmdstanr()

  source(file.path("04_Code", "R", "04_prepare_joint_model_data.R"))
  stan_data <- build_joint_interval_data("real")

  stan_file <- file.path(project_root, "04_Code", "Stan", "glw_joint_interval_dmr_independent_studentt.stan")
  model_dir <- file.path(project_root, "08_Model")
  dir.create(model_dir, showWarnings = FALSE, recursive = TRUE)

  mod <- cmdstanr::cmdstan_model(stan_file)
  fit <- mod$sample(
    data = stan_data,
    seed = seed,
    chains = chains,
    parallel_chains = parallel_chains,
    iter_warmup = iter_warmup,
    iter_sampling = iter_sampling,
    refresh = 100
  )

  fit_rds <- file.path(model_dir, paste0(model_prefix, "_fit.rds"))
  summary_csv <- file.path(model_dir, paste0(model_prefix, "_summary.csv"))
  draws_dir <- file.path(model_dir, paste0(model_prefix, "_draws"))
  dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

  saveRDS(fit, fit_rds)
  write.csv(fit$summary(), summary_csv, row.names = FALSE)
  fit$save_output_files(dir = draws_dir)

  cat("Saved renewed Stan fit:", fit_rds, "\n")
  cat("Saved renewed Stan summary:", summary_csv, "\n")
  cat("Saved renewed CmdStan CSV output files:", draws_dir, "\n")
  invisible(fit)
}

if (sys.nframe() == 0) {
  fit_stan_joint_model_independent()
}
