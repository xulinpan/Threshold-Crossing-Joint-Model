options(stringsAsFactors = FALSE)

source_local <- function(file_name) {
  candidates <- c(
    file.path("04_Code", "R", file_name),
    file_name,
    file.path(dirname(normalizePath(sys.frame(1)$ofile %||% "", winslash = "/", mustWork = FALSE)), file_name),
    file.path("D:/research2026/paper01_glw/04_Code/R", file_name)
  )
  for (candidate in candidates) {
    if (file.exists(candidate)) {
      source(candidate)
      return(invisible(TRUE))
    }
  }
  stop("Could not locate ", file_name, ". In RStudio, run setwd('D:/research2026/paper01_glw') first.")
}

`%||%` <- function(x, y) {
  if (is.null(x) || !nzchar(x)) y else x
}

source_local("01_clean_generate_model_data.R")

safe_scale <- function(x) {
  x <- as.numeric(x)
  center <- mean(x, na.rm = TRUE)
  spread <- stats::sd(x, na.rm = TRUE)
  if (is.na(spread) || spread == 0) spread <- 1
  (x - center) / spread
}

build_stan_data <- function(folder = PROCESSED_DIR, prefix = "real", covariate_adjusted = FALSE) {
  long <- read.csv(file.path(folder, paste0(prefix, "_longitudinal_analysis.csv")), stringsAsFactors = FALSE)
  interval <- read.csv(file.path(folder, paste0(prefix, "_interval_survival_analysis.csv")), stringsAsFactors = FALSE)
  patient <- read.csv(file.path(folder, paste0(prefix, "_patient_level_analysis.csv")), stringsAsFactors = FALSE)

  patient_nums <- sort(unique(long$patient_num))
  patient <- patient[match(patient_nums, patient$patient_num), , drop = FALSE]

  if (covariate_adjusted) {
    covariates <- patient[, c("age", "duration_years", "ph_baseline_pct", "sex_male"), drop = FALSE]
    complete <- stats::complete.cases(covariates)
    if (!all(complete)) {
      stop("Covariate-adjusted Stan data requires complete patient covariates. Use real_patient_level_complete_covariates.csv or set covariate_adjusted = FALSE.")
    }
    z <- cbind(
      age = safe_scale(covariates$age),
      duration_years = safe_scale(covariates$duration_years),
      ph_baseline_pct = safe_scale(covariates$ph_baseline_pct),
      sex_male = as.numeric(covariates$sex_male)
    )
  } else {
    z <- matrix(0, nrow = length(patient_nums), ncol = 0)
  }

  list(
    N_obs = nrow(long),
    N_pat = length(patient_nums),
    N_int = nrow(interval),
    id_obs = as.integer(long$patient_num),
    id_int = as.integer(interval$patient_num),
    y = safe_scale(long$log_mrd),
    t_obs = safe_scale(long$t_months),
    log_bcrabl = safe_scale(long$log_bcrabl),
    log_abl = safe_scale(long$log_abl),
    sample_bm = as.numeric(long$sample_bm),
    t_start = safe_scale(interval$t_start),
    t_end = safe_scale(interval$t_end),
    log_gap = safe_scale(interval$log_gap),
    event_interval = as.integer(interval$event_interval),
    K = ncol(z),
    Z = z,
    is_event_patient = as.integer(patient$dmr_event)
  )
}

if (sys.nframe() == 0) {
  stan_data_real <- build_stan_data(PROCESSED_DIR, "real", covariate_adjusted = FALSE)
  saveRDS(stan_data_real, file.path(PROCESSED_DIR, "stan_data_real_molecular_only.rds"))
  cat("Saved Stan data list to:", file.path(PROCESSED_DIR, "stan_data_real_molecular_only.rds"), "\n")
  cat("N_obs:", stan_data_real$N_obs, "N_pat:", stan_data_real$N_pat, "N_int:", stan_data_real$N_int, "\n")
}
