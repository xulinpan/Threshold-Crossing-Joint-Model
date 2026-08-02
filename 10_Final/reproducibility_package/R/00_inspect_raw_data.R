options(stringsAsFactors = FALSE)

project_root <- normalizePath(file.path(getwd()), winslash = "/", mustWork = TRUE)

read_utf8_csv <- function(path) {
  read.csv(
    file.path(project_root, path),
    fileEncoding = "UTF-8-BOM",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

files <- c("03_Data/Raw/glw.csv", "03_Data/Raw/PH+染色体-Table 1.csv")

for (f in files) {
  x <- read_utf8_csv(f)
  cat("\nFILE:", f, "\n")
  cat("DIM:", paste(dim(x), collapse = " x "), "\n")
  cat("COLUMNS:\n")
  print(names(x))
  cat("MISSING BY COLUMN:\n")
  print(colSums(is.na(x) | x == ""))
  cat("FIRST 3 ROWS:\n")
  print(utils::head(x, 3))
}

for (f in c(
  "03_Data/Processed/real_longitudinal_analysis.csv",
  "03_Data/Processed/real_patient_level_analysis.csv",
  "03_Data/Processed/real_interval_survival_analysis.csv"
)) {
  if (file.exists(file.path(project_root, f))) {
    x <- read.csv(file.path(project_root, f), stringsAsFactors = FALSE)
    cat("\nEXISTING:", f, "\n")
    cat("DIM:", paste(dim(x), collapse = " x "), "\n")
    if ("patient_id" %in% names(x)) {
      cat("PATIENTS:", length(unique(x$patient_id)), "\n")
    }
  }
}

cat("\nXLSX CONTENTS:\n")
print(utils::unzip(file.path(project_root, "03_Data/Raw/glw_data.xlsx"), list = TRUE)[, c("Name", "Length")])
