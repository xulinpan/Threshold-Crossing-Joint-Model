options(stringsAsFactors = FALSE)

project_root <- normalizePath("D:/research2026/paper01_glw", winslash = "/", mustWork = TRUE)
local_lib <- file.path(project_root, "04_Code", "R", "library")
cmdstan_dir <- file.path(project_root, "04_Code", "Stan", "cmdstan")

dir.create(local_lib, showWarnings = FALSE, recursive = TRUE)
dir.create(cmdstan_dir, showWarnings = FALSE, recursive = TRUE)
.libPaths(unique(c(local_lib, .libPaths())))

if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages(
    "cmdstanr",
    repos = c("https://mc-stan.org/r-packages/", "https://cloud.r-project.org"),
    lib = local_lib
  )
}

library(cmdstanr)

has_cmdstan <- FALSE
current_path <- tryCatch(cmdstanr::cmdstan_path(), error = function(e) "")
if (nzchar(current_path) && dir.exists(current_path)) {
  has_cmdstan <- TRUE
}

if (!has_cmdstan) {
  cmdstanr::install_cmdstan(dir = cmdstan_dir, cores = max(1, parallel::detectCores() - 1))
}

cat("cmdstanr library:", find.package("cmdstanr"), "\n")
cat("CmdStan path:", cmdstanr::cmdstan_path(), "\n")
