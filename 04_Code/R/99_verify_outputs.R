options(stringsAsFactors = FALSE)

source(file.path("04_Code", "R", "01_clean_generate_model_data.R"))

TABLE_DIR <- file.path(PROJECT_ROOT, "06_Tables")
FIGURE_DIR <- file.path(PROJECT_ROOT, "05_Figures")
LATEX_DIR <- file.path(PROJECT_ROOT, "02_LaTeX")

long <- read.csv(file.path(PROCESSED_DIR, "real_longitudinal_analysis.csv"), stringsAsFactors = FALSE)
intervals <- read.csv(file.path(PROCESSED_DIR, "real_interval_survival_analysis.csv"), stringsAsFactors = FALSE)
patient <- read.csv(file.path(PROCESSED_DIR, "real_patient_level_analysis.csv"), stringsAsFactors = FALSE)

required_long <- c(
  "patient_id", "patient_num", "visit_index", "t_months", "gap_months",
  "sample_type", "sample_bm", "bcr_abl_copy", "abl_copy", "ratio", "IS",
  "log_bcrabl", "log_abl", "log_IS", "log_mrd", "cmr", "dmr", "dmr_from_log_mrd"
)
required_interval <- c("patient_id", "patient_num", "visit_index", "t_start", "t_end", "gap_months", "log_gap", "event_interval")
required_patient <- c("patient_id", "patient_num", "dmr_event", "time_to_dmr_or_censor", "followup_months", "n_visits", "n_intervals")

checks <- list()
checks[["longitudinal_required_columns_present"]] <- all(required_long %in% names(long))
checks[["longitudinal_complete_required_fields"]] <- all(stats::complete.cases(long[, required_long]))
checks[["longitudinal_nonnegative_time"]] <- all(long$t_months >= 0 & long$gap_months >= 0)
checks[["interval_required_columns_present"]] <- all(required_interval %in% names(intervals))
checks[["interval_complete_required_fields"]] <- all(stats::complete.cases(intervals[, required_interval]))
checks[["interval_time_order"]] <- all(intervals$t_start <= intervals$t_end)
checks[["one_event_interval_per_patient_max"]] <- all(tapply(intervals$event_interval, intervals$patient_id, sum) <= 1)
checks[["patient_required_columns_present"]] <- all(required_patient %in% names(patient))
checks[["patient_complete_outcome_fields"]] <- all(stats::complete.cases(patient[, required_patient]))
checks[["patient_ids_are_coded"]] <- all(grepl("^P[0-9]{4}$", unique(c(long$patient_id, intervals$patient_id, patient$patient_id))))
checks[["no_name_columns_in_public_model_files"]] <- !any(grepl("name", names(long), ignore.case = TRUE)) &&
  !any(grepl("name", names(intervals), ignore.case = TRUE)) &&
  !any(grepl("name", names(patient), ignore.case = TRUE))

public_csvs <- list.files(PROCESSED_DIR, pattern = "\\.csv$", full.names = TRUE, recursive = FALSE)
public_tables <- list.files(TABLE_DIR, pattern = "\\.(csv|tex|md)$", full.names = TRUE, recursive = FALSE)
public_latex <- list.files(LATEX_DIR, pattern = "\\.tex$", full.names = TRUE, recursive = FALSE)
public_text_files <- c(public_csvs, public_tables, public_latex)
cjk_hits <- character()
for (f in public_text_files) {
  text <- paste(readLines(f, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  if (grepl("[\u4e00-\u9fff]", text)) cjk_hits <- c(cjk_hits, f)
}
checks[["no_cjk_text_in_public_csv_table_latex_outputs"]] <- length(cjk_hits) == 0

figure_files <- list.files(FIGURE_DIR, pattern = "\\.(pdf|png)$", full.names = TRUE)
table_files <- list.files(TABLE_DIR, pattern = "\\.(csv|tex|md)$", full.names = TRUE)
latex_pdf <- file.path(LATEX_DIR, "glw_data_analysis_report.pdf")
stan_rds <- file.path(PROCESSED_DIR, "stan_data_real_molecular_only.rds")
checks[["figures_exist_and_nonempty"]] <- length(figure_files) >= 10 && all(file.info(figure_files)$size > 0)
checks[["tables_exist_and_nonempty"]] <- length(table_files) >= 9 && all(file.info(table_files)$size > 0)
checks[["latex_pdf_exists_and_nonempty"]] <- file.exists(latex_pdf) && file.info(latex_pdf)$size > 0
checks[["stan_rds_exists_and_nonempty"]] <- file.exists(stan_rds) && file.info(stan_rds)$size > 0

verification <- data.frame(
  check = names(checks),
  passed = unlist(checks, use.names = FALSE),
  stringsAsFactors = FALSE
)

utils::write.csv(verification, file.path(TABLE_DIR, "verification_checks.csv"), row.names = FALSE)

report <- c(
  "# Verification checks",
  "",
  paste("- Longitudinal rows:", nrow(long)),
  paste("- Patients:", length(unique(long$patient_id))),
  paste("- Interval rows:", nrow(intervals)),
  paste("- Patient-level rows:", nrow(patient)),
  paste("- DMR event patients:", sum(patient$dmr_event, na.rm = TRUE)),
  "",
  "## Check results",
  "",
  paste0("- ", verification$check, ": ", ifelse(verification$passed, "PASS", "FAIL"))
)
if (length(cjk_hits) > 0) {
  report <- c(report, "", "## CJK privacy scan hits", "", cjk_hits)
}
writeLines(report, file.path(TABLE_DIR, "verification_report.md"), useBytes = TRUE)

print(verification)
if (!all(verification$passed)) {
  stop("One or more verification checks failed.")
}
