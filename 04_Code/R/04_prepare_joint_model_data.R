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
  script_file <- tryCatch(sys.frame(1)$ofile, error = function(e) "") %||% ""
  seed_paths <- c(
    Sys.getenv("GLW_PROJECT_ROOT", unset = ""),
    getwd(),
    if (nzchar(script_file)) dirname(script_file) else character(0),
    rstudio_active_dir(),
    "D:/research2026/paper01_glw"
  )
  candidates <- unique(unlist(lapply(seed_paths[nzchar(seed_paths)], parents_of), use.names = FALSE))
  for (candidate in candidates) {
    candidate <- normalizePath(candidate, winslash = "/", mustWork = FALSE)
    if (dir.exists(file.path(candidate, "03_Data", "Processed")) &&
        dir.exists(file.path(candidate, "04_Code", "R"))) {
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

ROOT_DIR <- find_project_root()
PROCESSED_DIR <- file.path(ROOT_DIR, "03_Data", "Processed")
MODEL_DIR <- file.path(ROOT_DIR, "08_Model")

dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

read_processed <- function(name) {
  read.csv(file.path(PROCESSED_DIR, name), stringsAsFactors = FALSE)
}

stop_if_missing <- function(data, required, data_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(data_name, " is missing required columns: ", paste(missing, collapse = ", "))
  }
}

build_joint_interval_data <- function(prefix = "real") {
  long <- read_processed(paste0(prefix, "_longitudinal_analysis.csv"))
  interval <- read_processed(paste0(prefix, "_interval_survival_analysis.csv"))
  patient <- read_processed(paste0(prefix, "_patient_level_analysis.csv"))

  stop_if_missing(
    long,
    c("patient_id", "patient_num", "visit_index", "t_months", "gap_months",
      "sample_bm", "log_mrd"),
    "longitudinal data"
  )
  stop_if_missing(
    interval,
    c("patient_id", "patient_num", "visit_index", "t_start", "t_end",
      "gap_months", "event_interval"),
    "interval data"
  )
  stop_if_missing(
    patient,
    c("patient_id", "patient_num", "dmr_event", "time_to_dmr_or_censor",
      "followup_months", "n_visits", "baseline_log_mrd", "age",
      "duration_years", "sex_male", "has_complete_covariates"),
    "patient-level data"
  )

  long <- long[order(long$patient_num, long$visit_index), , drop = FALSE]
  interval <- interval[order(interval$patient_num, interval$visit_index), , drop = FALSE]
  patient <- patient[order(patient$patient_num), , drop = FALSE]

  patient_nums <- sort(unique(long$patient_num))
  patient_lookup <- setNames(seq_along(patient_nums), patient_nums)

  if (!all(unique(interval$patient_num) %in% patient_nums)) {
    stop("Interval data contain patients absent from longitudinal data.")
  }

  long$model_patient_num <- as.integer(patient_lookup[as.character(long$patient_num)])
  interval$model_patient_num <- as.integer(patient_lookup[as.character(interval$patient_num)])

  event_counts <- aggregate(event_interval ~ patient_num, interval, sum)
  if (any(event_counts$event_interval > 1)) {
    stop("At least one patient has more than one event interval.")
  }

  floor_value <- -5

  stan_data <- list(
    N_obs = nrow(long),
    N_pat = length(patient_nums),
    N_int = nrow(interval),
    id_obs = long$model_patient_num,
    id_int = interval$model_patient_num,
    y = as.numeric(long$log_mrd),
    is_floor = as.integer(long$log_mrd <= floor_value),
    floor_value = floor_value,
    t_obs = pmax(as.numeric(long$t_months) / 12, 0),
    t_start = pmax(as.numeric(interval$t_start) / 12, 0),
    t_end = pmax(as.numeric(interval$t_end) / 12, 0),
    gap = pmax(as.numeric(interval$gap_months) / 12, 1e-6),
    sample_bm = as.numeric(long$sample_bm),
    event_interval = as.integer(interval$event_interval)
  )

  out_rds <- file.path(PROCESSED_DIR, paste0("stan_data_", prefix, "_joint_interval_dmr.rds"))
  saveRDS(stan_data, out_rds)

  summary_rows <- data.frame(
    metric = c(
      "patients",
      "longitudinal_observations",
      "at_risk_intervals",
      "dmr_event_patients",
      "censored_patients",
      "floor_observations_log_mrd_le_minus_5",
      "median_followup_months",
      "median_visits_per_patient",
      "median_visit_gap_months",
      "visit_gaps_over_6_months_pct",
      "visit_gaps_over_12_months_pct",
      "complete_core_covariates_pct"
    ),
    value = c(
      length(patient_nums),
      nrow(long),
      nrow(interval),
      sum(patient$dmr_event == 1),
      sum(patient$dmr_event == 0),
      sum(long$log_mrd <= floor_value, na.rm = TRUE),
      round(stats::median(patient$followup_months, na.rm = TRUE), 1),
      round(stats::median(patient$n_visits, na.rm = TRUE), 1),
      round(stats::median(long$gap_months, na.rm = TRUE), 1),
      round(mean(long$gap_months > 6, na.rm = TRUE) * 100, 1),
      round(mean(long$gap_months > 12, na.rm = TRUE) * 100, 1),
      round(mean(patient$has_complete_covariates == 1, na.rm = TRUE) * 100, 1)
    )
  )

  write.csv(
    summary_rows,
    file.path(MODEL_DIR, paste0("joint_interval_model_data_summary_", prefix, ".csv")),
    row.names = FALSE
  )

  invisible(stan_data)
}

if (sys.nframe() == 0) {
  stan_data <- build_joint_interval_data("real")
  cat("Saved joint interval-DMR Stan data.\n")
  cat("N_obs:", stan_data$N_obs, "N_pat:", stan_data$N_pat, "N_int:", stan_data$N_int, "\n")
}

