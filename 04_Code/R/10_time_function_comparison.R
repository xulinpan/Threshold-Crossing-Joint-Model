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
    "04_Code/R/10_time_function_comparison.R"
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
        dir.exists(file.path(candidate, "08_Model"))) {
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
LOCAL_R_LIB <- file.path(ROOT_DIR, "04_Code", "R", "library")
if (dir.exists(LOCAL_R_LIB)) {
  .libPaths(unique(c(LOCAL_R_LIB, .libPaths())))
}

required_packages <- c("nlme", "splines", "RColorBrewer")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

read_processed <- function(name) {
  read.csv(file.path(PROCESSED_DIR, name), stringsAsFactors = FALSE)
}

write_model_csv <- function(x, file_name) {
  write.csv(x, file.path(MODEL_DIR, file_name), row.names = FALSE)
}

safe_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

safe_logit <- function(p) {
  p <- safe_prob(p)
  log(p / (1 - p))
}

brier_score <- function(obs, pred) {
  mean((safe_prob(pred) - obs)^2)
}

calibration_intercept <- function(obs, pred) {
  obs <- as.integer(obs)
  pred <- safe_prob(pred)
  if (length(unique(obs)) < 2 || length(obs) < 5) return(NA_real_)
  fit <- try(stats::glm(obs ~ 1 + offset(safe_logit(pred)), family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  out <- unname(stats::coef(fit)[1])
  ifelse(is.finite(out), out, NA_real_)
}

calibration_slope <- function(obs, pred) {
  obs <- as.integer(obs)
  pred <- safe_prob(pred)
  if (length(unique(obs)) < 2 || length(obs) < 5 || stats::sd(pred) == 0) return(NA_real_)
  fit <- try(stats::glm(obs ~ safe_logit(pred), family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  out <- unname(stats::coef(fit)[2])
  ifelse(is.finite(out), out, NA_real_)
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(as.numeric(x), digits = digits, format = "f"))
}

safe_name <- function(x) {
  gsub("[^A-Za-z0-9_]", "_", x)
}

knot_term_names <- function(knots_years) {
  paste0("hinge_", safe_name(formatC(knots_years, digits = 2, format = "f")))
}

make_hinge_basis <- function(data, knots_years) {
  data$time_linear_years <- data$t_years
  for (idx in seq_along(knots_years)) {
    k <- knots_years[idx]
    nm <- knot_term_names(knots_years)[idx]
    data[[nm]] <- pmax(data$t_years - k, 0)
  }
  data
}

make_ns_basis <- function(data, ns_spec) {
  basis <- splines::ns(
    data$t_years,
    knots = ns_spec$knots,
    Boundary.knots = ns_spec$boundary
  )
  for (j in seq_len(ncol(basis))) {
    data[[paste0("ns", j)]] <- basis[, j]
  }
  data
}

fit_lme_safe <- function(data, fixed_formula, random_formula, fallback_random = ~ 1 | patient_id) {
  ctrl <- nlme::lmeControl(
    opt = "optim",
    maxIter = 200,
    msMaxIter = 200,
    returnObject = TRUE,
    tolerance = 1e-6
  )
  fit <- try(
    nlme::lme(
      fixed = fixed_formula,
      random = random_formula,
      data = data,
      method = "ML",
      control = ctrl,
      na.action = stats::na.omit
    ),
    silent = TRUE
  )
  if (!inherits(fit, "try-error")) {
    return(list(fit = fit, random_used = deparse(random_formula), fallback = FALSE))
  }

  fallback_fit <- nlme::lme(
    fixed = fixed_formula,
    random = fallback_random,
    data = data,
    method = "ML",
    control = ctrl,
    na.action = stats::na.omit
  )
  list(fit = fallback_fit, random_used = deparse(fallback_random), fallback = TRUE)
}

long <- read_processed("real_longitudinal_analysis.csv")
interval <- read_processed("real_interval_survival_analysis.csv")
patient <- read_processed("real_patient_level_analysis.csv")

long$patient_id <- factor(long$patient_id)
interval$patient_id <- factor(interval$patient_id, levels = levels(long$patient_id))
long$t_years <- pmax(as.numeric(long$t_months) / 12, 0)
long$lt <- log1p(long$t_years)
long$lt2 <- long$lt^2
long$floor_observation <- as.integer(long$log_mrd <= -5)

interval$t_start_years <- pmax(as.numeric(interval$t_start) / 12, 0)
interval$t_end_years <- pmax(as.numeric(interval$t_end) / 12, 0)
interval$t_years <- 0.5 * (interval$t_start_years + interval$t_end_years)
interval$gap_years <- pmax(as.numeric(interval$gap_months) / 12, 1e-6)
interval$lt <- log1p(interval$t_years)
interval$lt2 <- interval$lt^2
interval$sample_bm <- 0

patient <- patient[order(patient$patient_num), , drop = FALSE]

knots_months <- c(3, 6, 12, 24, 60)
knots_years <- knots_months / 12
long <- make_hinge_basis(long, knots_years)
interval <- make_hinge_basis(interval, knots_years)

ns_train <- splines::ns(long$t_years, df = 4)
ns_spec <- list(
  knots = attr(ns_train, "knots"),
  boundary = attr(ns_train, "Boundary.knots")
)
long <- make_ns_basis(long, ns_spec)
interval <- make_ns_basis(interval, ns_spec)

hinge_terms <- c(
  "time_linear_years",
  knot_term_names(knots_years)
)
ns_terms <- paste0("ns", seq_len(ncol(ns_train)))

time_models <- list(
  log_quadratic = list(
    label = "Log-quadratic in log(1+t)",
    fixed = stats::as.formula("log_mrd ~ lt + lt2 + sample_bm"),
    random = stats::as.formula("~ lt | patient_id"),
    terms = c("lt", "lt2")
  ),
  piecewise_linear = list(
    label = "Piecewise linear knots 3, 6, 12, 24, 60 months",
    fixed = stats::as.formula(paste("log_mrd ~", paste(c(hinge_terms, "sample_bm"), collapse = " + "))),
    random = stats::as.formula("~ time_linear_years | patient_id"),
    terms = hinge_terms
  ),
  penalized_spline_proxy = list(
    label = "Low-rank spline proxy (4 df)",
    fixed = stats::as.formula(paste("log_mrd ~", paste(c(ns_terms, "sample_bm"), collapse = " + "))),
    random = stats::as.formula("~ lt | patient_id"),
    terms = ns_terms
  )
)

fit_records <- list()
prediction_records <- list()

for (key in names(time_models)) {
  spec <- time_models[[key]]
  fitted <- fit_lme_safe(long, spec$fixed, spec$random)
  fit <- fitted$fit

  pred_long <- as.numeric(stats::predict(fit, newdata = long, level = 1))
  pred_long_fixed <- as.numeric(stats::predict(fit, newdata = long, level = 0))
  nonfloor <- long$floor_observation == 0

  interval$predicted_mrd <- as.numeric(stats::predict(fit, newdata = interval, level = 1))
  interval_fit <- stats::glm(
    event_interval ~ predicted_mrd + lt + log1p(gap_years),
    family = stats::binomial(),
    data = interval
  )
  interval$predicted_probability <- safe_prob(stats::predict(interval_fit, type = "response"))

  patient_probs <- vapply(patient$patient_num, function(id) {
    rows <- interval[interval$patient_num == id, , drop = FALSE]
    if (nrow(rows) == 0) return(NA_real_)
    1 - prod(1 - safe_prob(rows$predicted_probability))
  }, numeric(1))

  patient_eval <- data.frame(
    patient_id = patient$patient_id,
    patient_num = patient$patient_num,
    observed_event = as.integer(patient$dmr_event),
    predicted_probability = patient_probs,
    stringsAsFactors = FALSE
  )
  patient_eval <- patient_eval[is.finite(patient_eval$predicted_probability), , drop = FALSE]

  fit_records[[key]] <- data.frame(
    model_key = key,
    model = spec$label,
    n_fixed_time_parameters = length(spec$terms),
    random_effect_structure = fitted$random_used,
    random_slope_fallback = fitted$fallback,
    aic = stats::AIC(fit),
    bic = stats::BIC(fit),
    rmse_all_observed = sqrt(mean((long$log_mrd - pred_long)^2, na.rm = TRUE)),
    mae_all_observed = mean(abs(long$log_mrd - pred_long), na.rm = TRUE),
    rmse_nonfloor = sqrt(mean((long$log_mrd[nonfloor] - pred_long[nonfloor])^2, na.rm = TRUE)),
    mae_nonfloor = mean(abs(long$log_mrd[nonfloor] - pred_long[nonfloor]), na.rm = TRUE),
    floor_rate_observed = mean(long$floor_observation),
    mean_fixed_prediction_at_12m = mean(pred_long_fixed[abs(long$t_months - 12) < 3], na.rm = TRUE),
    interval_brier = brier_score(interval$event_interval, interval$predicted_probability),
    interval_calibration_intercept = calibration_intercept(interval$event_interval, interval$predicted_probability),
    interval_calibration_slope = calibration_slope(interval$event_interval, interval$predicted_probability),
    patient_brier = brier_score(patient_eval$observed_event, patient_eval$predicted_probability),
    patient_calibration_intercept = calibration_intercept(patient_eval$observed_event, patient_eval$predicted_probability),
    patient_calibration_slope = calibration_slope(patient_eval$observed_event, patient_eval$predicted_probability),
    note = "Screening mixed-model comparison; floor observations are treated as observed -5.0 here. Final claims require Bayesian refit with censoring.",
    stringsAsFactors = FALSE
  )

  prediction_records[[key]] <- data.frame(
    model_key = key,
    model = spec$label,
    patient_num = interval$patient_num,
    visit_index = interval$visit_index,
    t_mid_months = interval$t_years * 12,
    gap_months = interval$gap_years * 12,
    observed_event = interval$event_interval,
    predicted_mrd = interval$predicted_mrd,
    predicted_probability = interval$predicted_probability,
    stringsAsFactors = FALSE
  )

  time_models[[key]]$fit <- fit
}

summary_tab <- do.call(rbind, fit_records)
monotone_summary <- summary_tab[summary_tab$model_key == "penalized_spline_proxy", , drop = FALSE]
monotone_summary$model_key <- "monotone_soft_constraint"
monotone_summary$model <- "Soft monotone-decline constrained spline"
monotone_summary$random_effect_structure <- "~lt | patient_id; positive population increments penalized"
monotone_summary$random_slope_fallback <- NA
metric_cols <- c(
  "aic", "bic", "rmse_all_observed", "mae_all_observed",
  "rmse_nonfloor", "mae_nonfloor", "mean_fixed_prediction_at_12m",
  "interval_brier", "interval_calibration_intercept",
  "interval_calibration_slope", "patient_brier",
  "patient_calibration_intercept", "patient_calibration_slope"
)
monotone_summary[, metric_cols] <- NA_real_
monotone_summary$note <- paste(
  "Stan-only sensitivity candidate. Numeric performance requires a full Bayesian",
  "refit with the left-censoring likelihood; the plotted curve is a population",
  "projection of the spline proxy."
)
summary_tab <- rbind(summary_tab, monotone_summary)
write_model_csv(summary_tab, "time_function_comparison_summary.csv")
write_model_csv(do.call(rbind, prediction_records), "time_function_interval_predictions.csv")

grid <- data.frame(
  t_years = seq(0, max(long$t_years, na.rm = TRUE), length.out = 300),
  sample_bm = 0
)
grid$patient_id <- long$patient_id[1]
grid$lt <- log1p(grid$t_years)
grid$lt2 <- grid$lt^2
grid <- make_hinge_basis(grid, knots_years)
grid <- make_ns_basis(grid, ns_spec)

trend_records <- lapply(names(time_models), function(key) {
  pred <- as.numeric(stats::predict(time_models[[key]]$fit, newdata = grid, level = 0))
  data.frame(
    model_key = key,
    model = time_models[[key]]$label,
    t_months = grid$t_years * 12,
    predicted_log_mrd = pred,
    stringsAsFactors = FALSE
  )
})
trends <- do.call(rbind, trend_records)
spline_trend <- trends[trends$model_key == "penalized_spline_proxy", , drop = FALSE]
monotone_trend <- spline_trend[order(spline_trend$t_months), , drop = FALSE]
monotone_trend$model_key <- "monotone_soft_constraint"
monotone_trend$model <- "Soft monotone-decline spline projection"
monotone_trend$predicted_log_mrd <- cummin(monotone_trend$predicted_log_mrd)
trends <- rbind(trends, monotone_trend)
write_model_csv(trends, "time_function_population_trends.csv")

palette <- RColorBrewer::brewer.pal(4, "Dark2")
names(palette) <- c(
  "log_quadratic",
  "piecewise_linear",
  "penalized_spline_proxy",
  "monotone_soft_constraint"
)

plot_time_trends <- function(path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(path, width = 7.2, height = 4.8, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = 7.2, height = 4.8, units = "in", res = 600)
  }
  op <- graphics::par(
    mar = c(4.3, 4.5, 1.1, 0.6),
    mgp = c(2.5, 0.8, 0),
    las = 1
  )
  on.exit({
    graphics::par(op)
    grDevices::dev.off()
  }, add = TRUE)

  yrange <- range(c(trends$predicted_log_mrd, -4.5, -5), finite = TRUE)
  graphics::plot(
    NA,
    xlim = range(trends$t_months, finite = TRUE),
    ylim = yrange,
    xlab = "Treatment time (months)",
    ylab = "Predicted latent log-MRD",
    bty = "l"
  )
  graphics::abline(h = -4.5, col = "gray45", lty = 2, lwd = 1.2)
  graphics::abline(h = -5.0, col = "gray65", lty = 3, lwd = 1.2)
  for (key in names(palette)) {
    dat <- trends[trends$model_key == key, , drop = FALSE]
    graphics::lines(
      dat$t_months,
      dat$predicted_log_mrd,
      col = palette[[key]],
      lwd = ifelse(key == "monotone_soft_constraint", 2.0, 2.5),
      lty = ifelse(key == "monotone_soft_constraint", 5, 1)
    )
  }
  graphics::legend(
    "topright",
    legend = c("Log quadratic", "Piecewise linear", "Spline proxy", "Monotone projection", "DMR", "CMR/floor"),
    col = c(palette, "gray45", "gray65"),
    lwd = c(rep(2.4, 3), 2.0, 1.2, 1.2),
    lty = c(1, 1, 1, 5, 2, 3),
    bty = "n",
    cex = 0.83
  )
}

plot_time_trends(file.path(MODEL_DIR, "figure_12_time_function_population_trends.pdf"), "pdf")
plot_time_trends(file.path(MODEL_DIR, "figure_12_time_function_population_trends.png"), "png")

tex_rows <- vapply(seq_len(nrow(summary_tab)), function(i) {
  paste(
    summary_tab$model[i],
    summary_tab$n_fixed_time_parameters[i],
    fmt(summary_tab$aic[i], 1),
    fmt(summary_tab$rmse_nonfloor[i], 3),
    fmt(summary_tab$interval_brier[i], 3),
    fmt(summary_tab$interval_calibration_slope[i], 2),
    fmt(summary_tab$patient_brier[i], 3),
    fmt(summary_tab$patient_calibration_slope[i], 2),
    sep = " & "
  )
}, character(1))

time_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Screening comparison of longitudinal time functions. This table is an approximate mixed-model benchmark for time-function selection; floor observations are treated as observed at \\(-5.0\\) and final inference should be based on the Bayesian joint model with censoring.}",
  "\\label{tab:time-function-comparison}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{lrrrrrrr}",
  "\\toprule",
  "Time function & Time parameters & AIC & Non-floor RMSE & Interval Brier & Interval calibration slope & Patient Brier & Patient calibration slope \\\\",
  "\\midrule",
  paste0(tex_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "}%",
  "\\end{table}",
  ""
)
writeLines(time_table, file.path(MODEL_DIR, "table_08_time_function_comparison.tex"), useBytes = TRUE)

recommendation <- c(
  "# Time-Function Comparison and Recommendation",
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "## Scope",
  "",
  "This file compares candidate longitudinal time functions for the CML molecular-monitoring model. The empirical comparison is a screening benchmark using Gaussian mixed models; it is not a replacement for the Bayesian joint model, because assay-floor observations are treated as observed at -5.0 in this screening step. Any final manuscript claim about posterior predictive performance or DMR calibration requires refitting the selected Bayesian joint model with the left-censoring likelihood.",
  "",
  "## Empirical Screening Results",
  "",
  paste0(
    "- Best AIC in this screening benchmark: ",
    summary_tab$model[which.min(ifelse(is.finite(summary_tab$aic), summary_tab$aic, Inf))],
    " (AIC ", fmt(min(summary_tab$aic, na.rm = TRUE), 1), ")."
  ),
  paste0(
    "- Lowest non-floor RMSE: ",
    summary_tab$model[which.min(ifelse(is.finite(summary_tab$rmse_nonfloor), summary_tab$rmse_nonfloor, Inf))],
    " (RMSE ", fmt(min(summary_tab$rmse_nonfloor, na.rm = TRUE), 3), ")."
  ),
  paste0(
    "- Lowest interval-level Brier score in the two-stage calibration proxy: ",
    summary_tab$model[which.min(ifelse(is.finite(summary_tab$interval_brier), summary_tab$interval_brier, Inf))],
    " (Brier ", fmt(min(summary_tab$interval_brier, na.rm = TRUE), 3), ")."
  ),
  paste0(
    "- Lowest patient-level Brier score in the two-stage calibration proxy: ",
    summary_tab$model[which.min(ifelse(is.finite(summary_tab$patient_brier), summary_tab$patient_brier, Inf))],
    " (Brier ", fmt(min(summary_tab$patient_brier, na.rm = TRUE), 3), ")."
  ),
  "",
  "## Recommendation",
  "",
  "Recommended fitted primary time function for the current manuscript version: retain the current log-quadratic function unless the Bayesian joint model is fully refitted with the spline time trend. It is parsimonious, interpretable, has already passed HMC diagnostics in the renewed Stan fit, and performed best for non-floor RMSE and patient-level Brier score in this screening benchmark.",
  "",
  "Recommended renewal candidate if the model is refitted: a low-rank penalized spline for the population time trend, with independent patient-level random intercept and random log-time slope. The spline is the best flexibility candidate because it had the lowest AIC in the screening benchmark and can represent rapid early decline and later plateau without forcing one global quadratic curvature. It should replace the current primary model only if the full Bayesian refit preserves convergence, posterior predictive performance, and interval- and patient-level DMR calibration.",
  "",
  "The clinically knotted piecewise-linear model should be reported as an interpretability sensitivity analysis. A monotone-decline-constrained spline should be framed as a biological sensitivity analysis, preferably with a soft rather than hard population-level constraint so that genuine late increases in MRD are not masked.",
  "",
  "## Output Files",
  "",
  "- `time_function_comparison_summary.csv`",
  "- `time_function_interval_predictions.csv`",
  "- `time_function_population_trends.csv`",
  "- `table_08_time_function_comparison.tex`",
  "- `figure_12_time_function_population_trends.pdf`",
  "- `figure_12_time_function_population_trends.png`",
  "",
  "## Required Reanalysis Before Submission",
  "",
  "1. Refit the full Bayesian joint model using the selected penalized spline time function.",
  "2. Keep the assay-floor likelihood as left-censored at log-MRD <= -5.0.",
  "3. Compare posterior predictive checks against the current log-quadratic model, especially early follow-up, floor probability, and late follow-up behavior.",
  "4. Recompute interval-level and patient-level DMR calibration using posterior predicted event probabilities.",
  "5. Prefer the spline as primary only if convergence diagnostics are acceptable and the calibration/PPC checks are at least as good as the current primary model.",
  ""
)
writeLines(recommendation, file.path(MODEL_DIR, "time_function_modeling_recommendations.md"), useBytes = TRUE)

cat("Time-function comparison complete.\n")
cat("Summary:", file.path(MODEL_DIR, "time_function_comparison_summary.csv"), "\n")
cat("Figure:", file.path(MODEL_DIR, "figure_12_time_function_population_trends.png"), "\n")
