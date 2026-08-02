options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL) %||% "04_Code/R/run_all.R"
script_dir <- dirname(normalizePath(this_file, winslash = "/", mustWork = FALSE))
root_candidates <- c(
  getwd(),
  file.path(script_dir, "..", ".."),
  file.path(getwd(), ".."),
  file.path(getwd(), "..", "..")
)
project_root <- NULL
for (candidate in root_candidates) {
  candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
  if (dir.exists(file.path(candidate, "03_Data", "Raw"))) {
    project_root <- normalizePath(candidate, winslash = "/", mustWork = TRUE)
    break
  }
}
if (is.null(project_root)) stop("Could not locate project root.")
setwd(project_root)

source(file.path("04_Code", "R", "01_clean_generate_model_data.R"))
write_outputs()

source(file.path("04_Code", "R", "02_eda_figures_tables_latex.R"))

source(file.path("04_Code", "R", "03_build_stan_data.R"))
stan_data_real <- build_stan_data(PROCESSED_DIR, "real", covariate_adjusted = FALSE)
saveRDS(stan_data_real, file.path(PROCESSED_DIR, "stan_data_real_molecular_only.rds"))

source(file.path("04_Code", "R", "04_prepare_joint_model_data.R"))
build_joint_interval_data("real")

source(file.path("04_Code", "R", "05_fit_benchmark_models.R"))
fit_benchmark_models()

compile_latex <- Sys.which("pdflatex")
if (nzchar(compile_latex)) {
  old_wd <- getwd()
  setwd(file.path(project_root, "02_LaTeX"))
  on.exit(setwd(old_wd), add = TRUE)
  system2(compile_latex, c("-interaction=nonstopmode", "-halt-on-error", "glw_data_analysis_report.tex"))
  system2(compile_latex, c("-interaction=nonstopmode", "-halt-on-error", "glw_data_analysis_report.tex"))
  setwd(old_wd)
}

source(file.path("04_Code", "R", "99_verify_outputs.R"))

cat("Pipeline complete.\n")
cat("Processed data:", file.path(project_root, "03_Data", "Processed"), "\n")
cat("Figures:", file.path(project_root, "05_Figures"), "\n")
cat("Tables:", file.path(project_root, "06_Tables"), "\n")
cat("Model outputs:", file.path(project_root, "08_Model"), "\n")
cat("LaTeX report:", file.path(project_root, "02_LaTeX", "glw_data_analysis_report.tex"), "\n")
