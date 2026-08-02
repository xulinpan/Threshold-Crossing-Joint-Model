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

write_model_csv <- function(x, file_name) {
  write.csv(x, file.path(MODEL_DIR, file_name), row.names = FALSE)
}

make_patient_interval_data <- function(interval, patient) {
  by_patient <- split(interval, interval$patient_num)
  rows <- lapply(by_patient, function(dat) {
    dat <- dat[order(dat$t_end), , drop = FALSE]
    event_row <- dat[dat$event_interval == 1, , drop = FALSE]
    if (nrow(event_row) == 1) {
      data.frame(
        patient_num = dat$patient_num[1],
        left = event_row$t_start[1],
        right = event_row$t_end[1],
        interval_type = "interval",
        stringsAsFactors = FALSE
      )
    } else {
      data.frame(
        patient_num = dat$patient_num[1],
        left = max(dat$t_end, na.rm = TRUE),
        right = Inf,
        interval_type = "right_censored",
        stringsAsFactors = FALSE
      )
    }
  })
  out <- do.call(rbind, rows)
  merge(out, patient, by = "patient_num", all.x = TRUE, sort = FALSE)
}

tidy_survreg <- function(fit, model_name) {
  tab <- as.data.frame(summary(fit)$table)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  names(tab) <- sub("Std\\. Error", "std_error", names(tab))
  names(tab) <- sub("Value", "estimate", names(tab))
  names(tab) <- sub("z", "z_value", names(tab))
  names(tab) <- sub("p", "p_value", names(tab))
  data.frame(model = model_name, tab[, c("term", setdiff(names(tab), "term"))], check.names = FALSE)
}

tidy_lme <- function(fit, model_name) {
  tab <- as.data.frame(summary(fit)$tTable)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  data.frame(
    model = model_name,
    term = tab$term,
    estimate = tab$Value,
    std_error = tab$Std.Error,
    df = tab$DF,
    t_value = tab$`t-value`,
    p_value = tab$`p-value`,
    check.names = FALSE
  )
}

fit_benchmark_models <- function() {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("The survival package is required for interval-censored benchmark models.")
  }
  if (!requireNamespace("nlme", quietly = TRUE)) {
    stop("The nlme package is required for longitudinal benchmark models.")
  }
  if (!requireNamespace("splines", quietly = TRUE)) {
    stop("The splines package is required for longitudinal benchmark models.")
  }

  long <- read_processed("real_longitudinal_analysis.csv")
  interval <- read_processed("real_interval_survival_analysis.csv")
  patient <- read_processed("real_patient_level_analysis.csv")

  patient_interval <- make_patient_interval_data(interval, patient)
  patient_interval$right_for_model <- patient_interval$right
  patient_interval$t_event <- with(
    patient_interval,
    survival::Surv(left, right_for_model, type = "interval2")
  )

  molecular_fit <- survival::survreg(
    t_event ~ baseline_log_mrd,
    data = patient_interval,
    dist = "weibull"
  )

  clinical_dat <- patient_interval[
    stats::complete.cases(
      patient_interval[, c("baseline_log_mrd", "age", "sex_male", "duration_years")]
    ),
    ,
    drop = FALSE
  ]
  clinical_fit <- survival::survreg(
    t_event ~ baseline_log_mrd + age + sex_male + duration_years,
    data = clinical_dat,
    dist = "weibull"
  )

  long$t_years <- pmax(long$t_months / 12, 0)
  long$is_floor <- as.integer(long$log_mrd <= -5)

  longitudinal_fit <- try(
    nlme::lme(
      log_mrd ~ splines::ns(t_years, df = 3) + sample_bm,
      random = ~ 1 + t_years | patient_id,
      data = long,
      method = "REML",
      control = nlme::lmeControl(opt = "optim", msMaxIter = 200)
    ),
    silent = TRUE
  )

  longitudinal_model <- "random_intercept_slope"
  if (inherits(longitudinal_fit, "try-error")) {
    longitudinal_fit <- nlme::lme(
      log_mrd ~ splines::ns(t_years, df = 3) + sample_bm,
      random = ~ 1 | patient_id,
      data = long,
      method = "REML",
      control = nlme::lmeControl(opt = "optim", msMaxIter = 200)
    )
    longitudinal_model <- "random_intercept"
  }

  survreg_out <- rbind(
    tidy_survreg(molecular_fit, "interval_weibull_molecular"),
    tidy_survreg(clinical_fit, "interval_weibull_complete_case")
  )
  lme_out <- tidy_lme(longitudinal_fit, paste0("longitudinal_lme_", longitudinal_model))

  model_fit_summary <- data.frame(
    metric = c(
      "patients_in_interval_model",
      "patients_in_complete_case_interval_model",
      "longitudinal_observations",
      "longitudinal_patients",
      "longitudinal_model_random_effect",
      "molecular_interval_model_aic",
      "clinical_interval_model_aic",
      "longitudinal_model_aic"
    ),
    value = c(
      nrow(patient_interval),
      nrow(clinical_dat),
      nrow(long),
      length(unique(long$patient_id)),
      longitudinal_model,
      round(stats::AIC(molecular_fit), 2),
      round(stats::AIC(clinical_fit), 2),
      round(stats::AIC(longitudinal_fit), 2)
    )
  )

  write_model_csv(patient_interval, "patient_interval_outcome_dataset.csv")
  write_model_csv(survreg_out, "benchmark_interval_survival_coefficients.csv")
  write_model_csv(lme_out, "benchmark_longitudinal_coefficients.csv")
  write_model_csv(model_fit_summary, "benchmark_model_fit_summary.csv")

  sink(file.path(MODEL_DIR, "benchmark_model_summary.txt"))
  cat("GLW benchmark model summary\n")
  cat("Generated:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n\n")
  cat("Interval-censored Weibull model: molecular-only\n")
  print(summary(molecular_fit))
  cat("\nInterval-censored Weibull model: complete-case clinical adjustment\n")
  print(summary(clinical_fit))
  cat("\nLongitudinal mixed model:", longitudinal_model, "\n")
  print(summary(longitudinal_fit))
  sink()

  invisible(
    list(
      molecular_fit = molecular_fit,
      clinical_fit = clinical_fit,
      longitudinal_fit = longitudinal_fit,
      patient_interval = patient_interval
    )
  )
}

if (sys.nframe() == 0) {
  fit_benchmark_models()
  cat("Saved benchmark model outputs to:", MODEL_DIR, "\n")
}
