source(file.path("04_Code", "R", "01_clean_generate_model_data.R"))

long_raw <- read_longitudinal_raw()
patient_raw <- read_patient_raw()
key <- make_patient_key(long_raw, patient_raw)
build <- build_longitudinal(long_raw, key)

excluded <- build$raw_augmented[!build$raw_augmented$has_required_longitudinal, c(
  "row_id", "patient_name", "treatment_start_raw", "lab_date_raw",
  "sample_type", "bcr_abl_copy_raw", "abl_copy_raw", "ratio_raw",
  "is_raw", "log_mrd_raw", "visit_label_raw", "lab_date",
  "treatment_start_date_filled", "visit_label_months", "t_months"
)]

cat("Excluded longitudinal rows:", nrow(excluded), "\n")
print(excluded)

cat("\nPatients in raw longitudinal:", length(unique(build$raw_augmented$patient_id)), "\n")
cat("Patients in model longitudinal:", length(unique(build$data$patient_id)), "\n")
cat("Rows with missing time only:", sum(is.na(build$raw_augmented$t_months)), "\n")
