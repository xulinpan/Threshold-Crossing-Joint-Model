options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x) || !nzchar(x)) y else x
}

parents_of <- function(path, max_depth = 5) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  out <- path
  for (i in seq_len(max_depth)) {
    path <- dirname(path)
    out <- c(out, path)
  }
  unique(out)
}

rstudio_active_dir <- function() {
  if (!requireNamespace("rstudioapi", quietly = TRUE) || !rstudioapi::isAvailable()) {
    return(character(0))
  }
  path <- tryCatch(rstudioapi::getActiveDocumentContext()$path, error = function(e) "")
  if (!nzchar(path)) character(0) else dirname(path)
}

find_project_root <- function() {
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||%
    "04_Code/R/06_fit_stan_joint_model.R"
  seed_paths <- c(
    Sys.getenv("GLW_PROJECT_ROOT", unset = ""),
    getwd(),
    dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE)),
    rstudio_active_dir(),
    "D:/research2026/paper01_glw"
  )
  candidates <- unique(unlist(lapply(seed_paths[nzchar(seed_paths)], parents_of), use.names = FALSE))
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(candidate, "03_Data", "Processed")) &&
        dir.exists(file.path(candidate, "04_Code", "Stan"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop(
    "Could not locate project root. In RStudio, run:\n",
    "setwd('D:/research2026/paper01_glw')\n",
    "or:\n",
    "Sys.setenv(GLW_PROJECT_ROOT = 'D:/research2026/paper01_glw')"
  )
}

require_cmdstanr <- function() {
  project_root <- find_project_root()
  local_lib <- file.path(project_root, "04_Code", "R", "library")
  if (dir.exists(local_lib)) {
    .libPaths(unique(c(local_lib, .libPaths())))
  }

  if (!requireNamespace("cmdstanr", quietly = TRUE)) {
    stop(
      "cmdstanr is not installed. Install it first, then rerun this script.\n",
      "Suggested R commands:\n",
      "source('04_Code/R/install_cmdstan_workspace.R')\n",
      call. = FALSE
    )
  }

  local_cmdstan_root <- file.path(project_root, "04_Code", "Stan", "cmdstan")
  if (dir.exists(local_cmdstan_root)) {
    candidates <- c(
      local_cmdstan_root,
      list.dirs(local_cmdstan_root, recursive = FALSE, full.names = TRUE)
    )
    candidates <- candidates[dir.exists(file.path(candidates, "bin"))]
    if (length(candidates) > 0) {
      candidates <- candidates[order(file.info(candidates)$mtime, decreasing = TRUE)]
      try(cmdstanr::set_cmdstan_path(candidates[[1]]), silent = TRUE)
    }
  }

  cmdstan_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
  if (!nzchar(cmdstan_path) || !dir.exists(cmdstan_path)) {
    stop(
      "cmdstanr is installed, but CmdStan is not configured. Run:\n",
      "source('04_Code/R/install_cmdstan_workspace.R')\n",
      "or set the existing CmdStan path with cmdstanr::set_cmdstan_path().",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

as_env_int <- function(name, default) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) return(default)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed <= 0) default else parsed
}

fit_stan_joint_model <- function(
  chains = as_env_int("GLW_STAN_CHAINS", 4),
  parallel_chains = as_env_int("GLW_STAN_PARALLEL_CHAINS", min(chains, 4)),
  iter_warmup = as_env_int("GLW_STAN_WARMUP", 1000),
  iter_sampling = as_env_int("GLW_STAN_SAMPLING", 1000),
  seed = as_env_int("GLW_STAN_SEED", 20260709)
) {
  project_root <- find_project_root()
  setwd(project_root)

  require_cmdstanr()

  source(file.path("04_Code", "R", "04_prepare_joint_model_data.R"))
  stan_data <- build_joint_interval_data("real")

  stan_file <- file.path(project_root, "04_Code", "Stan", "glw_joint_interval_dmr.stan")
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

  fit_rds <- file.path(model_dir, "stan_joint_interval_dmr_fit.rds")
  summary_csv <- file.path(model_dir, "stan_joint_interval_dmr_summary.csv")
  draws_dir <- file.path(model_dir, "stan_joint_interval_dmr_draws")
  dir.create(draws_dir, showWarnings = FALSE, recursive = TRUE)

  saveRDS(fit, fit_rds)
  write.csv(fit$summary(), summary_csv, row.names = FALSE)
  fit$save_output_files(dir = draws_dir)

  cat("Saved Stan fit:", fit_rds, "\n")
  cat("Saved Stan summary:", summary_csv, "\n")
  cat("Saved CmdStan CSV output files:", draws_dir, "\n")
  invisible(fit)
}

if (sys.nframe() == 0) {
  fit_stan_joint_model()
}
