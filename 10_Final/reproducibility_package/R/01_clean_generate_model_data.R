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
    if (dir.exists(file.path(candidate, "03_Data", "Raw")) &&
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

PROJECT_ROOT <- find_project_root()
RAW_DIR <- file.path(PROJECT_ROOT, "03_Data", "Raw")
PROCESSED_DIR <- file.path(PROJECT_ROOT, "03_Data", "Processed")
PRIVATE_DIR <- file.path(PROCESSED_DIR, "private")
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(PRIVATE_DIR, recursive = TRUE, showWarnings = FALSE)

raw_longitudinal_file <- file.path(RAW_DIR, "glw.csv")
raw_patient_file <- file.path(RAW_DIR, "PH+\u67d3\u8272\u4f53-Table 1.csv")

trim_text <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "na", "null", "NULL")] <- NA_character_
  x
}

parse_number <- function(x) {
  x <- trim_text(x)
  x <- gsub(",", "", x)
  x <- gsub("%", "", x)
  suppressWarnings(as.numeric(x))
}

parse_percent <- function(x) {
  parse_number(x)
}

parse_lab_date <- function(x) {
  x <- trim_text(x)
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  idx <- !is.na(x)
  out[idx] <- as.Date(strptime(x[idx], format = "%m/%d/%y"))
  bad <- idx & is.na(out)
  if (any(bad)) {
    out[bad] <- suppressWarnings(as.Date(x[bad], tryFormats = c("%Y-%m-%d", "%Y/%m/%d", "%Y.%m.%d")))
  }
  out
}

parse_treatment_start <- function(x) {
  x <- trim_text(x)
  out <- as.Date(rep(NA_real_, length(x)), origin = "1970-01-01")
  for (i in seq_along(x)) {
    value <- x[i]
    if (is.na(value)) next
    value <- gsub("\\s+", "", value)
    m <- regexec("^([0-9]{1,2})\u6708([0-9]{1,2})\u65e5([0-9]{4})$", value)
    r <- regmatches(value, m)[[1]]
    if (length(r) == 4) {
      out[i] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(r[4]), as.integer(r[2]), as.integer(r[3])))
      next
    }
    m <- regexec("^([0-9]{4})\\s*([0-9]{1,2})\u6708([0-9]{1,2})\u65e5.*$", value)
    r <- regmatches(value, m)[[1]]
    if (length(r) == 4) {
      out[i] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(r[2]), as.integer(r[3]), as.integer(r[4])))
      next
    }
    m <- regexec("^([0-9]{4})\\s*([0-9]{1,2})\u6708.*$", value)
    r <- regmatches(value, m)[[1]]
    if (length(r) == 3) {
      out[i] <- as.Date(sprintf("%04d-%02d-01", as.integer(r[2]), as.integer(r[3])))
      next
    }
    m <- regexec("^([0-9]{4})[./-]([0-9]{1,2})[./-]([0-9]{1,2})$", value)
    r <- regmatches(value, m)[[1]]
    if (length(r) == 4) {
      out[i] <- as.Date(sprintf("%04d-%02d-%02d", as.integer(r[2]), as.integer(r[3]), as.integer(r[4])))
      next
    }
  }
  out
}

parse_visit_label_months <- function(x) {
  x <- trim_text(x)
  out <- rep(NA_real_, length(x))
  for (i in seq_along(x)) {
    value <- x[i]
    if (is.na(value)) next
    value <- gsub("\\s+", "", value)
    if (grepl("\u670d\u836f\u524d", value)) {
      out[i] <- 0
    } else if (grepl("([0-9.]+)\u4e2a\u6708", value)) {
      m <- regexec("([0-9.]+)\u4e2a\u6708", value)
      out[i] <- parse_number(regmatches(value, m)[[1]][2])
    } else if (grepl("([0-9.]+)\u6708", value)) {
      m <- regexec("([0-9.]+)\u6708", value)
      out[i] <- parse_number(regmatches(value, m)[[1]][2])
    } else if (grepl("([0-9.]+)\u5e74", value)) {
      m <- regexec("([0-9.]+)\u5e74", value)
      out[i] <- 12 * parse_number(regmatches(value, m)[[1]][2])
    } else {
      num <- parse_number(value)
      if (!is.na(num)) out[i] <- 12 * num
    }
  }
  out
}

first_non_missing <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA)
  x[[1]]
}

read_longitudinal_raw <- function() {
  raw <- read.csv(
    raw_longitudinal_file,
    fileEncoding = "UTF-8-BOM",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  raw <- raw[, seq_len(min(11, ncol(raw))), drop = FALSE]
  names(raw) <- c(
    "treatment_start_raw", "duration_raw", "patient_name", "lab_date_raw",
    "sample_type", "bcr_abl_copy_raw", "abl_copy_raw", "ratio_raw",
    "is_raw", "log_mrd_raw", "visit_label_raw"
  )[seq_len(ncol(raw))]
  raw$row_id <- seq_len(nrow(raw))
  raw
}

read_patient_raw <- function() {
  raw <- read.csv(
    raw_patient_file,
    fileEncoding = "UTF-8-BOM",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("", "NA")
  )
  raw <- raw[, seq_len(min(22, ncol(raw))), drop = FALSE]
  names(raw) <- c(
    "source_seq", "patient_name", "sex", "age", "duration_years",
    "ph_baseline_pct", "ph_3m_pct", "ph_6m_pct", "ph_9m_pct",
    "ph_12m_pct", "ph_18m_pct", "ph_2y_pct", "ph_3y_pct",
    "ph_4y_pct", "ph_5y_pct", "ph_6y_pct", "ph_7y_pct",
    "ph_8y_pct", "ph_9y_pct", "ph_10y_pct", "ph_11y_pct",
    "ph_12y_pct"
  )[seq_len(ncol(raw))]
  raw$row_id <- seq_len(nrow(raw))
  raw
}

make_patient_key <- function(long_raw, patient_raw) {
  names_long <- trim_text(long_raw$patient_name)
  names_patient <- trim_text(patient_raw$patient_name)
  all_names <- unique(c(names_long[!is.na(names_long)], names_patient[!is.na(names_patient)]))
  data.frame(
    patient_id = sprintf("P%04d", seq_along(all_names)),
    patient_name = all_names,
    in_longitudinal_raw = all_names %in% names_long,
    in_patient_raw = all_names %in% names_patient,
    stringsAsFactors = FALSE
  )
}

build_longitudinal <- function(long_raw, key) {
  x <- long_raw
  x$patient_name <- trim_text(x$patient_name)
  x$sample_type <- toupper(trim_text(x$sample_type))
  x$lab_date <- parse_lab_date(x$lab_date_raw)
  x$treatment_start_date <- parse_treatment_start(x$treatment_start_raw)
  x$visit_label_months <- if ("visit_label_raw" %in% names(x)) parse_visit_label_months(x$visit_label_raw) else NA_real_
  x$bcr_abl_copy <- parse_number(x$bcr_abl_copy_raw)
  x$abl_copy <- parse_number(x$abl_copy_raw)
  x$ratio <- parse_number(x$ratio_raw)
  x$IS <- parse_number(x$is_raw)
  x$log_mrd <- parse_number(x$log_mrd_raw)

  start_by_name <- tapply(x$treatment_start_date, x$patient_name, first_non_missing)
  x$treatment_start_date_filled <- as.Date(start_by_name[x$patient_name], origin = "1970-01-01")
  x$t_months_from_dates <- as.numeric(x$lab_date - x$treatment_start_date_filled) / 30.4375
  x$t_months <- x$t_months_from_dates
  missing_time <- is.na(x$t_months) & !is.na(x$visit_label_months)
  x$t_months[missing_time] <- x$visit_label_months[missing_time]
  baseline_time <- !is.na(x$visit_label_months) & x$visit_label_months == 0
  x$t_months[baseline_time] <- 0
  x$t_months <- ifelse(!is.na(x$t_months) & abs(x$t_months) < 0.5, 0, x$t_months)

  x <- merge(x, key[, c("patient_id", "patient_name")], by = "patient_name", all.x = TRUE, sort = FALSE)

  required <- c("patient_id", "lab_date", "sample_type", "bcr_abl_copy", "abl_copy", "ratio", "IS", "log_mrd", "t_months")
  x$has_required_longitudinal <- stats::complete.cases(x[, required]) & x$t_months >= 0

  kept <- x[x$has_required_longitudinal, , drop = FALSE]
  kept <- kept[order(kept$patient_id, kept$t_months, kept$lab_date, kept$row_id), , drop = FALSE]
  kept$patient_num <- match(kept$patient_id, sort(unique(kept$patient_id)))
  kept$visit_index <- ave(kept$row_id, kept$patient_id, FUN = seq_along)
  kept$prev_t_months <- ave(kept$t_months, kept$patient_id, FUN = function(v) c(0, head(v, -1)))
  kept$gap_months <- pmax(0, kept$t_months - kept$prev_t_months)
  kept$log_gap <- log1p(kept$gap_months)
  kept$sample_bm <- as.integer(kept$sample_type == "BM")
  kept$log_bcrabl <- log1p(kept$bcr_abl_copy)
  kept$log_abl <- log1p(kept$abl_copy)
  kept$log_IS <- log1p(kept$IS)
  kept$cmr <- as.integer(kept$log_mrd <= -5.0)
  kept$dmr <- as.integer(kept$log_mrd <= -4.5)
  kept$dmr_from_log_mrd <- kept$dmr

  output <- kept[, c(
    "patient_id", "patient_num", "visit_index", "t_months", "gap_months",
    "sample_type", "sample_bm", "bcr_abl_copy", "abl_copy", "ratio", "IS",
    "log_bcrabl", "log_abl", "log_IS", "log_mrd", "cmr", "dmr",
    "dmr_from_log_mrd"
  )]
  row.names(output) <- NULL

  audit <- data.frame(
    item = c(
      "raw_longitudinal_rows",
      "raw_longitudinal_patients",
      "dropped_missing_or_invalid_core_fields",
      "model_ready_longitudinal_rows",
      "model_ready_longitudinal_patients"
    ),
    value = c(
      nrow(x),
      length(unique(x$patient_id[!is.na(x$patient_id)])),
      sum(!x$has_required_longitudinal),
      nrow(output),
      length(unique(output$patient_id))
    )
  )

  list(data = output, audit = audit, raw_augmented = x)
}

build_patient_covariates <- function(patient_raw, key) {
  x <- patient_raw
  x$patient_name <- trim_text(x$patient_name)
  x$sex <- trim_text(x$sex)
  x$age <- parse_number(x$age)
  x$duration_years <- parse_number(x$duration_years)
  pct_cols <- grep("^ph_", names(x), value = TRUE)
  for (nm in pct_cols) x[[nm]] <- parse_percent(x[[nm]])

  valid <- !is.na(x$patient_name) & x$sex %in% c("\u7537", "\u5973") & !is.na(x$age)
  x <- x[valid, , drop = FALSE]
  x <- merge(x, key[, c("patient_id", "patient_name")], by = "patient_name", all.x = TRUE, sort = FALSE)
  x$sex_male <- as.integer(x$sex == "\u7537")
  x <- x[order(x$patient_id), , drop = FALSE]

  out <- x[, c(
    "patient_id", "age", "duration_years", "ph_baseline_pct", "ph_3m_pct",
    "ph_6m_pct", "ph_9m_pct", "ph_12m_pct", "ph_18m_pct", "ph_2y_pct",
    "ph_3y_pct", "ph_4y_pct", "ph_5y_pct", "ph_6y_pct", "ph_7y_pct",
    "ph_8y_pct", "ph_9y_pct", "ph_10y_pct", "ph_11y_pct", "ph_12y_pct",
    "sex_male"
  )]
  row.names(out) <- NULL
  out
}

build_interval_and_patient <- function(longitudinal, covariates) {
  pieces <- vector("list", length(unique(longitudinal$patient_id)))
  names(pieces) <- unique(longitudinal$patient_id)
  for (pid in names(pieces)) {
    g <- longitudinal[longitudinal$patient_id == pid, , drop = FALSE]
    g <- g[order(g$t_months, g$visit_index), , drop = FALSE]
    previous_events <- c(0, head(cumsum(g$dmr_from_log_mrd), -1))
    g$at_risk <- previous_events == 0
    g$event_interval <- as.integer(g$dmr_from_log_mrd == 1 & previous_events == 0)
    pieces[[pid]] <- g[g$at_risk, , drop = FALSE]
  }
  intervals_raw <- do.call(rbind, pieces)
  intervals <- intervals_raw[, c(
    "patient_id", "patient_num", "visit_index", "t_months", "gap_months",
    "event_interval"
  )]
  intervals$log_gap <- log1p(intervals$gap_months)
  intervals$t_start <- pmax(0, intervals$t_months - intervals$gap_months)
  intervals$t_end <- intervals$t_months
  intervals <- intervals[, c("patient_id", "patient_num", "visit_index", "t_start", "t_end", "gap_months", "log_gap", "event_interval")]
  row.names(intervals) <- NULL

  patient_rows <- lapply(split(longitudinal, longitudinal$patient_id), function(g) {
    g <- g[order(g$t_months, g$visit_index), , drop = FALSE]
    event_rows <- g[g$dmr_from_log_mrd == 1, , drop = FALSE]
    event <- nrow(event_rows) > 0
    data.frame(
      patient_id = g$patient_id[1],
      patient_num = g$patient_num[1],
      dmr_event = as.integer(event),
      time_to_dmr_or_censor = if (event) event_rows$t_months[1] else max(g$t_months),
      followup_months = max(g$t_months),
      n_visits = nrow(g),
      n_intervals = sum(intervals$patient_id == g$patient_id[1]),
      baseline_log_mrd = g$log_mrd[1],
      last_log_mrd = g$log_mrd[nrow(g)],
      min_log_mrd = min(g$log_mrd, na.rm = TRUE),
      ever_cmr = as.integer(any(g$cmr == 1)),
      stringsAsFactors = FALSE
    )
  })
  patient <- do.call(rbind, patient_rows)
  row.names(patient) <- NULL
  patient <- merge(patient, covariates, by = "patient_id", all.x = TRUE, sort = FALSE)
  patient <- patient[order(patient$patient_num), , drop = FALSE]
  patient$has_complete_covariates <- as.integer(stats::complete.cases(patient[, c(
    "age", "duration_years", "ph_baseline_pct", "sex_male"
  )]))

  complete_covariates <- patient[patient$has_complete_covariates == 1, , drop = FALSE]
  row.names(complete_covariates) <- NULL

  list(intervals = intervals, patient = patient, complete_patient = complete_covariates)
}

write_data_dictionary <- function(cleaning_audit, longitudinal, intervals, patient) {
  lines <- c(
    "# GLW CML model-ready datasets",
    "",
    paste("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "The files in this folder are regenerated from `03_Data/Raw/glw.csv` and",
    "`03_Data/Raw/PH+染色体-Table 1.csv`. Public modeling files use coded",
    "patient identifiers only. The reversible name map is isolated under",
    "`03_Data/Processed/private/patient_key.csv` and should not be shared.",
    "",
    "## Core derived variables",
    "",
    "- `patient_id`: privacy-preserving coded identifier.",
    "- `patient_num`: integer patient index for modeling software.",
    "- `t_months`: months from imatinib start to molecular monitoring date.",
    "- `gap_months`: months since previous visit for that patient.",
    "- `log_mrd`: observed LOG-MRD.",
    "- `dmr` and `dmr_from_log_mrd`: 1 when `log_mrd <= -4.5`.",
    "- `cmr`: 1 when `log_mrd <= -5.0`.",
    "- `event_interval`: first at-risk interval ending in DMR.",
    "",
    "## Output files",
    "",
    "- `real_longitudinal_analysis.csv`: complete longitudinal molecular visits.",
    "- `real_interval_survival_analysis.csv`: at-risk intervals through first DMR.",
    "- `real_patient_level_analysis.csv`: patient outcomes plus baseline covariates.",
    "- `real_patient_level_complete_covariates.csv`: subset with complete core covariates.",
    "- `data_cleaning_audit.csv`: row counts and exclusion counts.",
    "",
    "## Current regenerated counts",
    "",
    paste("- Longitudinal observations:", nrow(longitudinal)),
    paste("- Patients:", length(unique(longitudinal$patient_id))),
    paste("- At-risk survival intervals:", nrow(intervals)),
    paste("- DMR event patients:", sum(patient$dmr_event, na.rm = TRUE)),
    paste("- Censored patients:", sum(patient$dmr_event == 0, na.rm = TRUE)),
    paste("- Complete core baseline covariates:", sum(patient$has_complete_covariates == 1, na.rm = TRUE)),
    "",
    "## Cleaning audit",
    "",
    paste(cleaning_audit$item, cleaning_audit$value, sep = ": ")
  )
  writeLines(lines, file.path(PROCESSED_DIR, "README_data_dictionary.md"), useBytes = TRUE)
}

write_outputs <- function() {
  long_raw <- read_longitudinal_raw()
  patient_raw <- read_patient_raw()
  key <- make_patient_key(long_raw, patient_raw)
  longitudinal_build <- build_longitudinal(long_raw, key)
  covariates <- build_patient_covariates(patient_raw, key)
  survival_build <- build_interval_and_patient(longitudinal_build$data, covariates)

  patient_key <- key[order(key$patient_id), , drop = FALSE]
  utils::write.csv(patient_key, file.path(PRIVATE_DIR, "patient_key.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  utils::write.csv(longitudinal_build$data, file.path(PROCESSED_DIR, "real_longitudinal_analysis.csv"), row.names = FALSE)
  utils::write.csv(survival_build$intervals, file.path(PROCESSED_DIR, "real_interval_survival_analysis.csv"), row.names = FALSE)
  utils::write.csv(survival_build$patient, file.path(PROCESSED_DIR, "real_patient_level_analysis.csv"), row.names = FALSE)
  utils::write.csv(survival_build$complete_patient, file.path(PROCESSED_DIR, "real_patient_level_complete_covariates.csv"), row.names = FALSE)
  utils::write.csv(covariates, file.path(PROCESSED_DIR, "clean_patient_covariates.csv"), row.names = FALSE)

  patient_audit <- data.frame(
    item = c(
      "raw_patient_table_rows",
      "usable_patient_covariate_rows",
      "patient_level_analysis_rows",
      "complete_core_covariate_rows",
      "interval_survival_rows",
      "dmr_event_patients",
      "censored_patients"
    ),
    value = c(
      nrow(patient_raw),
      nrow(covariates),
      nrow(survival_build$patient),
      nrow(survival_build$complete_patient),
      nrow(survival_build$intervals),
      sum(survival_build$patient$dmr_event, na.rm = TRUE),
      sum(survival_build$patient$dmr_event == 0, na.rm = TRUE)
    )
  )
  cleaning_audit <- rbind(longitudinal_build$audit, patient_audit)
  utils::write.csv(cleaning_audit, file.path(PROCESSED_DIR, "data_cleaning_audit.csv"), row.names = FALSE)

  missing_longitudinal <- data.frame(
    variable = c("patient_name", "lab_date", "sample_type", "bcr_abl_copy", "abl_copy", "ratio", "IS", "log_mrd", "t_months"),
    missing_or_invalid = c(
      sum(is.na(longitudinal_build$raw_augmented$patient_name)),
      sum(is.na(longitudinal_build$raw_augmented$lab_date)),
      sum(is.na(longitudinal_build$raw_augmented$sample_type)),
      sum(is.na(longitudinal_build$raw_augmented$bcr_abl_copy)),
      sum(is.na(longitudinal_build$raw_augmented$abl_copy)),
      sum(is.na(longitudinal_build$raw_augmented$ratio)),
      sum(is.na(longitudinal_build$raw_augmented$IS)),
      sum(is.na(longitudinal_build$raw_augmented$log_mrd)),
      sum(is.na(longitudinal_build$raw_augmented$t_months) | longitudinal_build$raw_augmented$t_months < 0)
    )
  )
  utils::write.csv(missing_longitudinal, file.path(PROCESSED_DIR, "longitudinal_missingness.csv"), row.names = FALSE)

  write_data_dictionary(cleaning_audit, longitudinal_build$data, survival_build$intervals, survival_build$patient)

  invisible(list(
    longitudinal = longitudinal_build$data,
    intervals = survival_build$intervals,
    patient = survival_build$patient,
    cleaning_audit = cleaning_audit
  ))
}

if (sys.nframe() == 0) {
  result <- write_outputs()
  cat("Wrote model-ready datasets to:", PROCESSED_DIR, "\n")
  print(result$cleaning_audit)
}
