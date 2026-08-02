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
  this_file <- tryCatch(sys.frame(1)$ofile, error = function(e) "") %||%
    "04_Code/R/11_prepare_visit_process_data.R"
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

ROOT_DIR <- find_project_root()
PROCESSED_DIR <- file.path(ROOT_DIR, "03_Data", "Processed")
MODEL_DIR <- file.path(ROOT_DIR, "08_Model")
dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

read_processed <- function(name) {
  read.csv(file.path(PROCESSED_DIR, name), stringsAsFactors = FALSE)
}

make_midpoint_quadrature <- function(patient_rows, nodes_per_patient = 32) {
  pieces <- lapply(seq_len(nrow(patient_rows)), function(i) {
    c_years <- pmax(as.numeric(patient_rows$time_to_dmr_or_censor[i]) / 12, 1e-6)
    width <- c_years / nodes_per_patient
    data.frame(
      patient_id = patient_rows$patient_id[i],
      patient_num = patient_rows$patient_num[i],
      t_quad = (seq_len(nodes_per_patient) - 0.5) * width,
      w_quad = rep(width, nodes_per_patient),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, pieces)
}

build_visit_process_data <- function(nodes_per_patient = 32) {
  stan_data <- readRDS(file.path(PROCESSED_DIR, "stan_data_real_joint_interval_dmr.rds"))
  long <- read_processed("real_longitudinal_analysis.csv")
  patient <- read_processed("real_patient_level_analysis.csv")

  long <- long[order(long$patient_num, long$visit_index), , drop = FALSE]
  patient <- patient[order(patient$patient_num), , drop = FALSE]

  patient_nums <- sort(unique(long$patient_num))
  model_lookup <- setNames(seq_along(patient_nums), patient_nums)
  patient$model_patient_num <- as.integer(model_lookup[as.character(patient$patient_num)])
  long$model_patient_num <- as.integer(model_lookup[as.character(long$patient_num)])

  outcome_time <- patient[, c("patient_num", "time_to_dmr_or_censor"), drop = FALSE]
  long <- merge(long, outcome_time, by = "patient_num", all.x = TRUE, sort = FALSE)
  long$t_years <- as.numeric(long$t_months) / 12
  long$censor_years <- as.numeric(long$time_to_dmr_or_censor) / 12

  visit_rows <- long[
    is.finite(long$t_years) &
      long$t_years > 0 &
      long$t_years <= long$censor_years + 1e-8,
    ,
    drop = FALSE
  ]
  visit_rows <- visit_rows[order(visit_rows$model_patient_num, visit_rows$t_years), , drop = FALSE]

  quad <- make_midpoint_quadrature(patient, nodes_per_patient = nodes_per_patient)
  quad$model_patient_num <- as.integer(model_lookup[as.character(quad$patient_num)])
  quad <- quad[order(quad$model_patient_num, quad$t_quad), , drop = FALSE]

  visit_data <- stan_data
  visit_data$N_visit <- nrow(visit_rows)
  visit_data$N_quad <- nrow(quad)
  visit_data$id_visit <- as.integer(visit_rows$model_patient_num)
  visit_data$t_visit <- as.numeric(visit_rows$t_years)
  visit_data$id_quad <- as.integer(quad$model_patient_num)
  visit_data$t_quad <- as.numeric(quad$t_quad)
  visit_data$w_quad <- as.numeric(quad$w_quad)

  out_rds <- file.path(PROCESSED_DIR, "stan_data_real_joint_interval_dmr_visit_process.rds")
  saveRDS(visit_data, out_rds)

  visit_counts <- aggregate(t_years ~ patient_num, visit_rows, length)
  names(visit_counts)[2] <- "n_visit_process"
  patient_counts <- merge(
    patient[, c("patient_id", "patient_num", "time_to_dmr_or_censor", "followup_months", "dmr_event")],
    visit_counts,
    by = "patient_num",
    all.x = TRUE
  )
  patient_counts$n_visit_process[is.na(patient_counts$n_visit_process)] <- 0

  summary_rows <- data.frame(
    metric = c(
      "patients",
      "visit_process_observed_visits",
      "quadrature_nodes",
      "quadrature_nodes_per_patient",
      "median_visit_process_visits_per_patient",
      "median_time_window_months",
      "dmr_event_patients",
      "censored_patients",
      "risk_window_definition"
    ),
    value = c(
      length(patient_nums),
      nrow(visit_rows),
      nrow(quad),
      nodes_per_patient,
      stats::median(patient_counts$n_visit_process),
      round(stats::median(patient_counts$time_to_dmr_or_censor), 2),
      sum(patient_counts$dmr_event == 1),
      sum(patient_counts$dmr_event == 0),
      "from treatment initiation to first documented DMR or censoring; excludes t=0 baseline visits"
    ),
    stringsAsFactors = FALSE
  )

  write.csv(summary_rows, file.path(MODEL_DIR, "visit_process_data_summary.csv"), row.names = FALSE)
  write.csv(patient_counts, file.path(MODEL_DIR, "visit_process_patient_counts.csv"), row.names = FALSE)

  invisible(visit_data)
}

if (sys.nframe() == 0) {
  visit_data <- build_visit_process_data()
  cat("Saved visit-process Stan data.\n")
  cat("N_visit:", visit_data$N_visit, "N_quad:", visit_data$N_quad, "\n")
}
