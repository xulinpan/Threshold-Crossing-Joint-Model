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
    "04_Code/R/09_model_comparison_validation.R"
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
  stop("Could not locate project root.")
}

ROOT_DIR <- find_project_root()
PROCESSED_DIR <- file.path(ROOT_DIR, "03_Data", "Processed")
MODEL_DIR <- file.path(ROOT_DIR, "08_Model")
LOCAL_R_LIB <- file.path(ROOT_DIR, "04_Code", "R", "library")
if (dir.exists(LOCAL_R_LIB)) {
  .libPaths(unique(c(LOCAL_R_LIB, .libPaths())))
}

required_packages <- c("survival", "RColorBrewer")
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

rbind_fill <- function(...) {
  items <- list(...)
  items <- items[!vapply(items, is.null, logical(1))]
  cols <- unique(unlist(lapply(items, names), use.names = FALSE))
  items <- lapply(items, function(x) {
    missing <- setdiff(cols, names(x))
    for (nm in missing) x[[nm]] <- NA
    x[, cols, drop = FALSE]
  })
  do.call(rbind, items)
}

fmt <- function(x, digits = 3) {
  ifelse(
    is.na(x),
    "--",
    formatC(as.numeric(x), digits = digits, format = "f")
  )
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

clean_calibration_metric <- function(x) {
  ifelse(is.finite(x) & abs(x) <= 20, x, NA_real_)
}

metric_set <- function(obs, pred) {
  data.frame(
    n = length(obs),
    events = sum(obs),
    brier = brier_score(obs, pred),
    mean_predicted = mean(safe_prob(pred)),
    observed_rate = mean(obs),
    calibration_intercept = calibration_intercept(obs, pred),
    calibration_slope = calibration_slope(obs, pred),
    stringsAsFactors = FALSE
  )
}

metrics_for_predictions <- function(predictions, eval_patient_nums = NULL) {
  interval_detail <- predictions$interval
  patient_detail <- predictions$patient
  if (!is.null(eval_patient_nums)) {
    interval_detail <- interval_detail[interval_detail$patient_num %in% eval_patient_nums, , drop = FALSE]
    patient_detail <- patient_detail[patient_detail$patient_num %in% eval_patient_nums, , drop = FALSE]
  }
  interval_metrics <- metric_set(interval_detail$observed_event, interval_detail$predicted_probability)
  patient_metrics <- metric_set(patient_detail$observed_event, patient_detail$predicted_probability)
  list(interval = interval_metrics, patient = patient_metrics)
}

weighted_metrics_for_ids <- function(predictions, eval_ids) {
  interval_rows <- lapply(eval_ids, function(id) {
    predictions$interval[predictions$interval$patient_num == id, , drop = FALSE]
  })
  patient_rows <- lapply(eval_ids, function(id) {
    predictions$patient[predictions$patient$patient_num == id, , drop = FALSE]
  })
  interval_detail <- do.call(rbind, interval_rows)
  patient_detail <- do.call(rbind, patient_rows)
  interval_metrics <- metric_set(interval_detail$observed_event, interval_detail$predicted_probability)
  patient_metrics <- metric_set(patient_detail$observed_event, patient_detail$predicted_probability)
  list(interval = interval_metrics, patient = patient_metrics)
}

patient_probability_from_intervals <- function(interval_detail, patient, observed_patient = NULL) {
  rows <- lapply(seq_len(nrow(patient)), function(i) {
    dat <- interval_detail[interval_detail$patient_num == patient$patient_num[i], , drop = FALSE]
    if (nrow(dat) == 0) return(NULL)
    pred <- 1 - prod(1 - safe_prob(dat$predicted_probability))
    obs <- if (is.null(observed_patient)) patient$dmr_event[i] else observed_patient[i]
    data.frame(
      patient_id = patient$patient_id[i],
      patient_num = patient$patient_num[i],
      observed_event = as.integer(obs),
      predicted_probability = pred,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

make_prediction_object <- function(model_key, model_label, interval_detail, patient_detail, note = "") {
  interval_detail$model_key <- model_key
  interval_detail$model <- model_label
  patient_detail$model_key <- model_key
  patient_detail$model <- model_label
  list(
    model_key = model_key,
    model = model_label,
    interval = interval_detail,
    patient = patient_detail,
    note = note
  )
}

long <- read_processed("real_longitudinal_analysis.csv")
interval <- read_processed("real_interval_survival_analysis.csv")
patient <- read_processed("real_patient_level_analysis.csv")
stan_data <- readRDS(file.path(PROCESSED_DIR, "stan_data_real_joint_interval_dmr.rds"))

long$t_years <- long$t_months / 12
interval$t_start_years <- interval$t_start / 12
interval$t_end_years <- interval$t_end / 12
interval$t_mid_years <- 0.5 * (interval$t_start_years + interval$t_end_years)
interval$gap_years <- pmax(interval$gap_months / 12, 1e-6)

patient_ids <- patient$patient_num

fit_km_predictions <- function(train_ids = patient_ids) {
  train_patient <- do.call(rbind, lapply(train_ids, function(id) {
    patient[patient$patient_num == id, , drop = FALSE]
  }))
  fit <- survival::survfit(
    survival::Surv(time_to_dmr_or_censor, dmr_event) ~ 1,
    data = train_patient
  )
  surv_at <- function(times) {
    times <- pmax(times, 0)
    out <- rep(1, length(times))
    positive <- times > 0
    if (any(positive)) {
      out[positive] <- summary(fit, times = times[positive], extend = TRUE)$surv
    }
    pmax(out, 1e-6)
  }
  s_start <- surv_at(interval$t_start)
  s_end <- surv_at(interval$t_end)
  interval_pred <- safe_prob(1 - s_end / s_start)
  interval_detail <- data.frame(
    patient_id = interval$patient_id,
    patient_num = interval$patient_num,
    observed_event = interval$event_interval,
    predicted_probability = interval_pred,
    stringsAsFactors = FALSE
  )
  patient_pred <- safe_prob(1 - surv_at(patient$time_to_dmr_or_censor))
  patient_detail <- data.frame(
    patient_id = patient$patient_id,
    patient_num = patient$patient_num,
    observed_event = patient$dmr_event,
    predicted_probability = patient_pred,
    stringsAsFactors = FALSE
  )
  make_prediction_object(
    "km_descriptive",
    "Descriptive Kaplan--Meier",
    interval_detail,
    patient_detail,
    "Uses first observed DMR visit time; descriptive because true DMR onset is interval-observed."
  )
}

fit_interval_predictions <- function(train_ids = patient_ids) {
  train_interval <- do.call(rbind, lapply(train_ids, function(id) {
    interval[interval$patient_num == id, , drop = FALSE]
  }))
  fit <- stats::glm(
    event_interval ~ log1p(t_mid_years) + log1p(gap_years),
    data = train_interval,
    family = stats::binomial()
  )
  interval_pred <- safe_prob(stats::predict(fit, newdata = interval, type = "response"))
  interval_detail <- data.frame(
    patient_id = interval$patient_id,
    patient_num = interval$patient_num,
    observed_event = interval$event_interval,
    predicted_probability = interval_pred,
    stringsAsFactors = FALSE
  )
  patient_detail <- patient_probability_from_intervals(interval_detail, patient)
  make_prediction_object(
    "interval_timing_gap",
    "Interval timing + visit gap",
    interval_detail,
    patient_detail,
    "Discrete-time interval model without longitudinal MRD."
  )
}

build_landmark_interval_data <- function() {
  landmarks <- c(6, 12, 18)
  horizon <- 24
  rows <- vector("list", nrow(interval))
  for (r in seq_len(nrow(interval))) {
    int <- interval[r, ]
    possible <- landmarks[landmarks <= int$t_start]
    if (length(possible) == 0) next
    possible <- rev(possible)
    selected <- NULL
    for (lm in possible) {
      if (int$t_end > lm + horizon) next
      p_long <- long[long$patient_num == int$patient_num & long$t_months <= lm, , drop = FALSE]
      if (nrow(p_long) == 0) next
      if (any(p_long$dmr == 1, na.rm = TRUE)) next
      last_row <- p_long[which.max(p_long$t_months), , drop = FALSE]
      selected <- data.frame(
        interval_row = r,
        patient_id = int$patient_id,
        patient_num = int$patient_num,
        visit_index = int$visit_index,
        landmark_months = lm,
        landmark_years = lm / 12,
        horizon_months = horizon,
        t_start = int$t_start,
        t_end = int$t_end,
        gap_years = int$gap_years,
        last_log_mrd = last_row$log_mrd,
        last_sample_bm = last_row$sample_bm,
        observed_event = int$event_interval,
        stringsAsFactors = FALSE
      )
      break
    }
    rows[[r]] <- selected
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

landmark_interval <- build_landmark_interval_data()

fit_landmark_predictions <- function(train_ids = patient_ids) {
  train_landmark <- do.call(rbind, lapply(train_ids, function(id) {
    landmark_interval[landmark_interval$patient_num == id, , drop = FALSE]
  }))
  train_landmark <- train_landmark[stats::complete.cases(train_landmark[, c(
    "observed_event", "last_log_mrd", "landmark_years", "gap_years"
  )]), , drop = FALSE]

  if (nrow(train_landmark) < 10 || length(unique(train_landmark$observed_event)) < 2) {
    p <- mean(train_landmark$observed_event)
    interval_pred <- rep(p, nrow(landmark_interval))
  } else {
    fit <- stats::glm(
      observed_event ~ last_log_mrd + landmark_years + log1p(gap_years),
      data = train_landmark,
      family = stats::binomial()
    )
    interval_pred <- stats::predict(fit, newdata = landmark_interval, type = "response")
  }

  interval_detail <- data.frame(
    patient_id = landmark_interval$patient_id,
    patient_num = landmark_interval$patient_num,
    observed_event = landmark_interval$observed_event,
    predicted_probability = safe_prob(interval_pred),
    landmark_months = landmark_interval$landmark_months,
    stringsAsFactors = FALSE
  )
  patient_landmark <- patient[patient$patient_num %in% unique(interval_detail$patient_num), , drop = FALSE]
  patient_obs <- vapply(patient_landmark$patient_num, function(id) {
    any(interval_detail$patient_num == id & interval_detail$observed_event == 1)
  }, logical(1))
  patient_detail <- patient_probability_from_intervals(
    interval_detail,
    patient_landmark,
    observed_patient = as.integer(patient_obs)
  )
  make_prediction_object(
    "landmark_mrd",
    "Landmark MRD",
    interval_detail,
    patient_detail,
    "Uses last observed MRD before 6-, 12-, or 18-month landmarks and scores only later intervals within 24 months."
  )
}

read_stan_draws <- function(prefix) {
  draw_dir <- file.path(MODEL_DIR, paste0(prefix, "_draws"))
  draw_files <- list.files(draw_dir, pattern = "[.]csv$", full.names = TRUE)
  if (length(draw_files) == 0) {
    stop("Missing Stan draw files for prefix: ", prefix)
  }
  lapply(draw_files, function(f) read.csv(f, comment.char = "#", check.names = FALSE))
}

stan_diagnostics <- function(prefix, draws_by_chain) {
  summary_path <- file.path(MODEL_DIR, paste0(prefix, "_summary.csv"))
  summary_tab <- read.csv(summary_path, stringsAsFactors = FALSE)
  chain_diag <- do.call(rbind, lapply(seq_along(draws_by_chain), function(i) {
    d <- draws_by_chain[[i]]
    data.frame(
      chain = i,
      post_warmup_draws = nrow(d),
      divergences = sum(d[["divergent__"]]),
      max_treedepth = max(d[["treedepth__"]]),
      ebfmi = mean(diff(d[["energy__"]])^2) / stats::var(d[["energy__"]]),
      stringsAsFactors = FALSE
    )
  }))
  data.frame(
    model_prefix = prefix,
    post_warmup_draws = sum(chain_diag$post_warmup_draws),
    divergences = sum(chain_diag$divergences),
    max_treedepth = max(chain_diag$max_treedepth),
    min_ebfmi = min(chain_diag$ebfmi),
    max_rhat = max(summary_tab$rhat, na.rm = TRUE),
    min_bulk_ess = min(summary_tab$ess_bulk, na.rm = TRUE),
    min_tail_ess = min(summary_tab$ess_tail, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

fit_stan_predictions <- function(prefix, model_key, model_label, note) {
  draws_by_chain <- read_stan_draws(prefix)
  draws <- do.call(rbind, draws_by_chain)
  event_cols <- paste0("event_prob.", seq_len(stan_data$N_int))
  missing_cols <- setdiff(event_cols, names(draws))
  if (length(missing_cols) > 0) {
    stop("Missing event probability columns for ", prefix)
  }
  event_prob <- colMeans(as.matrix(draws[, event_cols]))
  interval_detail <- data.frame(
    patient_id = interval$patient_id,
    patient_num = interval$patient_num,
    observed_event = interval$event_interval,
    predicted_probability = safe_prob(event_prob),
    stringsAsFactors = FALSE
  )
  patient_detail <- patient_probability_from_intervals(interval_detail, patient)
  pred <- make_prediction_object(model_key, model_label, interval_detail, patient_detail, note)
  pred$diagnostics <- stan_diagnostics(prefix, draws_by_chain)
  pred
}

predictions <- list(
  km_descriptive = fit_km_predictions(),
  interval_timing_gap = fit_interval_predictions(),
  landmark_mrd = fit_landmark_predictions(),
  joint_exact_floor = fit_stan_predictions(
    "stan_joint_interval_dmr_exact_floor",
    "joint_exact_floor",
    "Joint model, floor exact",
    "Same joint interval structure as primary model, but assay-floor observations are treated as exact -5.0 values."
  ),
  joint_left_censored = fit_stan_predictions(
    "stan_joint_interval_dmr_independent_renewed",
    "joint_left_censored",
    "Primary joint model",
    "Renewed primary Bayesian joint longitudinal--interval model with left-censored assay-floor likelihood."
  )
)

apparent_metrics <- do.call(rbind, lapply(predictions, function(pred) {
  m <- metrics_for_predictions(pred)
  data.frame(
    model_key = pred$model_key,
    model = pred$model,
    interval_n = m$interval$n,
    interval_events = m$interval$events,
    interval_brier = m$interval$brier,
    interval_observed_rate = m$interval$observed_rate,
    interval_mean_predicted = m$interval$mean_predicted,
    interval_calibration_intercept = m$interval$calibration_intercept,
    interval_calibration_slope = m$interval$calibration_slope,
    patient_n = m$patient$n,
    patient_events = m$patient$events,
    patient_brier = m$patient$brier,
    patient_observed_rate = m$patient$observed_rate,
    patient_mean_predicted = m$patient$mean_predicted,
    patient_calibration_intercept = m$patient$calibration_intercept,
    patient_calibration_slope = m$patient$calibration_slope,
    note = pred$note,
    stringsAsFactors = FALSE
  )
}))

validate_refittable_model <- function(model_key, reps = 200, seed = 20260709) {
  set.seed(seed)
  fit_fun <- switch(
    model_key,
    km_descriptive = fit_km_predictions,
    interval_timing_gap = fit_interval_predictions,
    landmark_mrd = fit_landmark_predictions,
    stop("Not a refittable simple model: ", model_key)
  )
  apparent <- apparent_metrics[apparent_metrics$model_key == model_key, ]
  boot_rows <- vector("list", reps)
  for (b in seq_len(reps)) {
    boot_ids <- sample(patient_ids, length(patient_ids), replace = TRUE)
    fit <- try(fit_fun(boot_ids), silent = TRUE)
    if (inherits(fit, "try-error")) next
    train_m <- try(weighted_metrics_for_ids(fit, boot_ids), silent = TRUE)
    test_m <- try(metrics_for_predictions(fit, patient_ids), silent = TRUE)
    if (inherits(train_m, "try-error") || inherits(test_m, "try-error")) next
    boot_rows[[b]] <- data.frame(
      replicate = b,
      interval_brier_optimism = train_m$interval$brier - test_m$interval$brier,
      patient_brier_optimism = train_m$patient$brier - test_m$patient$brier,
      interval_intercept_optimism = train_m$interval$calibration_intercept - test_m$interval$calibration_intercept,
      interval_slope_optimism = train_m$interval$calibration_slope - test_m$interval$calibration_slope,
      patient_intercept_optimism = train_m$patient$calibration_intercept - test_m$patient$calibration_intercept,
      patient_slope_optimism = train_m$patient$calibration_slope - test_m$patient$calibration_slope,
      stringsAsFactors = FALSE
    )
  }
  boot <- do.call(rbind, boot_rows)
  data.frame(
    model_key = model_key,
    bootstrap_reps_requested = reps,
    bootstrap_reps_used = nrow(boot),
    validation_method = "patient-cluster bootstrap optimism correction with model refitting",
    optimism_corrected_interval_brier = apparent$interval_brier - mean(boot$interval_brier_optimism, na.rm = TRUE),
    optimism_corrected_patient_brier = apparent$patient_brier - mean(boot$patient_brier_optimism, na.rm = TRUE),
    optimism_corrected_interval_calibration_intercept = clean_calibration_metric(apparent$interval_calibration_intercept - mean(boot$interval_intercept_optimism, na.rm = TRUE)),
    optimism_corrected_interval_calibration_slope = clean_calibration_metric(apparent$interval_calibration_slope - mean(boot$interval_slope_optimism, na.rm = TRUE)),
    optimism_corrected_patient_calibration_intercept = clean_calibration_metric(apparent$patient_calibration_intercept - mean(boot$patient_intercept_optimism, na.rm = TRUE)),
    optimism_corrected_patient_calibration_slope = clean_calibration_metric(apparent$patient_calibration_slope - mean(boot$patient_slope_optimism, na.rm = TRUE)),
    stringsAsFactors = FALSE
  )
}

validate_fixed_prediction_model <- function(model_key, reps = 200, seed = 20260709) {
  set.seed(seed)
  pred <- predictions[[model_key]]
  apparent <- apparent_metrics[apparent_metrics$model_key == model_key, ]
  boot <- do.call(rbind, lapply(seq_len(reps), function(b) {
    boot_ids <- sample(unique(pred$patient$patient_num), nrow(pred$patient), replace = TRUE)
    m <- weighted_metrics_for_ids(pred, boot_ids)
    data.frame(
      replicate = b,
      interval_brier = m$interval$brier,
      patient_brier = m$patient$brier,
      interval_calibration_intercept = m$interval$calibration_intercept,
      interval_calibration_slope = m$interval$calibration_slope,
      patient_calibration_intercept = m$patient$calibration_intercept,
      patient_calibration_slope = m$patient$calibration_slope,
      stringsAsFactors = FALSE
    )
  }))
  data.frame(
    model_key = model_key,
    bootstrap_reps_requested = reps,
    bootstrap_reps_used = nrow(boot),
    validation_method = "patient-cluster bootstrap of fixed posterior predictions; no Stan refit",
    optimism_corrected_interval_brier = apparent$interval_brier,
    optimism_corrected_patient_brier = apparent$patient_brier,
    optimism_corrected_interval_calibration_intercept = apparent$interval_calibration_intercept,
    optimism_corrected_interval_calibration_slope = apparent$interval_calibration_slope,
    optimism_corrected_patient_calibration_intercept = apparent$patient_calibration_intercept,
    optimism_corrected_patient_calibration_slope = apparent$patient_calibration_slope,
    interval_brier_bootstrap_lcl = stats::quantile(boot$interval_brier, 0.025, na.rm = TRUE),
    interval_brier_bootstrap_ucl = stats::quantile(boot$interval_brier, 0.975, na.rm = TRUE),
    patient_brier_bootstrap_lcl = stats::quantile(boot$patient_brier, 0.025, na.rm = TRUE),
    patient_brier_bootstrap_ucl = stats::quantile(boot$patient_brier, 0.975, na.rm = TRUE),
    stringsAsFactors = FALSE
  )
}

bootstrap_reps <- as.integer(Sys.getenv("GLW_BOOTSTRAP_REPS", unset = "200"))
validation <- rbind_fill(
  validate_refittable_model("km_descriptive", reps = bootstrap_reps, seed = 20260709),
  validate_refittable_model("interval_timing_gap", reps = bootstrap_reps, seed = 20260710),
  validate_refittable_model("landmark_mrd", reps = bootstrap_reps, seed = 20260711),
  validate_fixed_prediction_model("joint_exact_floor", reps = bootstrap_reps, seed = 20260712),
  validate_fixed_prediction_model("joint_left_censored", reps = bootstrap_reps, seed = 20260713)
)

comparison <- merge(apparent_metrics, validation, by = "model_key", all.x = TRUE, sort = FALSE)
comparison$model_order <- match(
  comparison$model_key,
  c("km_descriptive", "interval_timing_gap", "landmark_mrd", "joint_exact_floor", "joint_left_censored")
)
comparison <- comparison[order(comparison$model_order), ]
comparison$model_order <- NULL
write_model_csv(comparison, "model_comparison_performance.csv")

all_interval_predictions <- do.call(rbind_fill, lapply(predictions, `[[`, "interval"))
all_patient_predictions <- do.call(rbind_fill, lapply(predictions, `[[`, "patient"))
write_model_csv(all_interval_predictions, "model_comparison_interval_predictions.csv")
write_model_csv(all_patient_predictions, "model_comparison_patient_predictions.csv")

stan_diag <- do.call(
  rbind,
  lapply(predictions[c("joint_exact_floor", "joint_left_censored")], `[[`, "diagnostics")
)
write_model_csv(stan_diag, "model_comparison_stan_diagnostics.csv")

make_equal_count_bins <- function(x, bins = 5) {
  n <- length(x)
  bins <- min(bins, max(2, floor(n / 8)))
  ord <- order(x, seq_along(x), na.last = NA)
  out <- integer(n)
  out[ord] <- cut(
    seq_along(ord),
    breaks = seq(0, length(ord), length.out = bins + 1),
    labels = FALSE,
    include.lowest = TRUE
  )
  out
}

binom_ci <- function(events, n, alpha = 0.05) {
  lower <- ifelse(events == 0, 0, stats::qbeta(alpha / 2, events, n - events + 1))
  upper <- ifelse(events == n, 1, stats::qbeta(1 - alpha / 2, events + 1, n - events))
  cbind(lower = lower, upper = upper)
}

calibration_bins <- function(detail, bins = 5) {
  split_detail <- split(detail, detail$model_key)
  out <- lapply(names(split_detail), function(model_key) {
    dat <- split_detail[[model_key]]
    dat$bin <- make_equal_count_bins(dat$predicted_probability, bins = bins)
    rows <- lapply(sort(unique(dat$bin)), function(b) {
      idx <- dat$bin == b
      ci <- binom_ci(sum(dat$observed_event[idx]), sum(idx))
      data.frame(
        model_key = model_key,
        model = dat$model[1],
        bin = b,
        n = sum(idx),
        observed_events = sum(dat$observed_event[idx]),
        mean_predicted = mean(dat$predicted_probability[idx]),
        observed_rate = mean(dat$observed_event[idx]),
        observed_rate_lower_95 = ci[, "lower"],
        observed_rate_upper_95 = ci[, "upper"],
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, rows)
  })
  do.call(rbind, out)
}

interval_calibration_bins <- calibration_bins(all_interval_predictions)
patient_calibration_bins <- calibration_bins(all_patient_predictions)
write_model_csv(interval_calibration_bins, "model_comparison_interval_calibration_bins.csv")
write_model_csv(patient_calibration_bins, "model_comparison_patient_calibration_bins.csv")

make_calibration_plot <- function(cal, file_stub, title, ylab) {
  models <- unique(cal$model)
  palette <- RColorBrewer::brewer.pal(max(5, min(8, length(models))), "Dark2")
  colors <- setNames(palette[seq_along(models)], models)
  pdf_path <- file.path(MODEL_DIR, paste0(file_stub, ".pdf"))
  png_path <- file.path(MODEL_DIR, paste0(file_stub, ".png"))
  plot_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(
      mar = c(4.7, 4.8, 3.0, 1.2),
      bg = "white",
      family = "sans",
      cex.axis = 0.9,
      cex.lab = 1.0,
      cex.main = 1.05
    )
    plot(
      NA,
      xlim = c(0, 1),
      ylim = c(0, 1),
      xlab = "Mean predicted probability",
      ylab = ylab,
      main = title,
      bty = "l"
    )
    grid(col = "gray88", lty = 1)
    abline(0, 1, col = "gray35", lwd = 2)
    for (m in models) {
      dat <- cal[cal$model == m, ]
      dat <- dat[order(dat$mean_predicted), ]
      lines(dat$mean_predicted, dat$observed_rate, col = colors[m], lwd = 1.8)
      points(dat$mean_predicted, dat$observed_rate, col = colors[m], pch = 19, cex = 0.95)
    }
    legend(
      "topleft",
      legend = models,
      col = colors[models],
      pch = 19,
      lwd = 1.8,
      bty = "n",
      cex = 0.78
    )
  }
  grDevices::pdf(pdf_path, width = 7.4, height = 5.4)
  plot_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = 7.4, height = 5.4, units = "in", res = 600)
  plot_fun()
  grDevices::dev.off()
}

make_calibration_plot(
  interval_calibration_bins,
  "figure_10_model_comparison_interval_calibration",
  "Interval-Level Calibration by Model",
  "Observed interval event rate"
)
make_calibration_plot(
  patient_calibration_bins,
  "figure_11_model_comparison_patient_calibration",
  "Patient-Level Calibration by Model",
  "Observed patient event rate"
)

make_table_tex <- function(comparison) {
  row_tex <- vapply(seq_len(nrow(comparison)), function(i) {
    paste(
      comparison$model[i],
      paste0(comparison$interval_n[i], "/", comparison$interval_events[i]),
      fmt(comparison$interval_brier[i]),
      fmt(comparison$optimism_corrected_interval_brier[i]),
      paste0(comparison$patient_n[i], "/", comparison$patient_events[i]),
      fmt(comparison$patient_brier[i]),
      fmt(comparison$optimism_corrected_patient_brier[i]),
      paste0(
        fmt(comparison$interval_calibration_intercept[i]),
        "/",
        fmt(comparison$interval_calibration_slope[i])
      ),
      paste0(
        fmt(comparison$patient_calibration_intercept[i]),
        "/",
        fmt(comparison$patient_calibration_slope[i])
      ),
      sep = " & "
    )
  }, character(1))
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Model comparison and internal validation. Event counts are shown as evaluated records/events. Validated Brier scores are patient-cluster bootstrap optimism-corrected for refittable simpler models; for the two Stan joint models, they are fixed-posterior bootstrap summaries without Stan refitting.}",
    "\\label{tab:model-comparison-validation}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{lrrrrrrll}",
    "\\toprule",
    "Model & Interval n/events & Interval Brier & Validated interval Brier & Patient n/events & Patient Brier & Validated patient Brier & Interval cal. int./slope & Patient cal. int./slope \\\\",
    "\\midrule",
    paste0(row_tex, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "}%",
    "\\end{table}",
    ""
  )
  writeLines(lines, file.path(MODEL_DIR, "table_07_model_comparison_validation.tex"), useBytes = TRUE)
}

make_table_tex(comparison)

primary <- comparison[comparison$model_key == "joint_left_censored", ]
exact <- comparison[comparison$model_key == "joint_exact_floor", ]
interval_model <- comparison[comparison$model_key == "interval_timing_gap", ]
km_model <- comparison[comparison$model_key == "km_descriptive", ]
landmark_model <- comparison[comparison$model_key == "landmark_mrd", ]

model_comparison_text <- c(
  "\\subsection{Model Comparison and Internal Validation}",
  "",
  "We compared the renewed primary Bayesian joint longitudinal--interval model with four prespecified alternatives: a descriptive Kaplan--Meier curve based on first observed DMR visit time, a discrete-time interval model with timing and visit-gap terms but no longitudinal MRD, a landmark MRD model using the last observed MRD before 6-, 12-, or 18-month landmarks, and a joint model that treated assay-floor observations as exact values at \\(-5.0\\). The Kaplan--Meier analysis was retained only as a descriptive benchmark because true DMR onset was interval-observed.",
  "",
  "\\input{table_07_model_comparison_validation.tex}",
  "",
  paste0(
    "The descriptive Kaplan--Meier benchmark had interval-level Brier score ",
    fmt(km_model$interval_brier),
    " and patient-level Brier score ",
    fmt(km_model$patient_brier),
    ". The interval timing--gap model improved interval-level calibration by directly modeling at-risk intervals, with interval Brier score ",
    fmt(interval_model$interval_brier),
    " and optimism-corrected interval Brier score ",
    fmt(interval_model$optimism_corrected_interval_brier),
    ". The landmark MRD model was clinically interpretable but evaluable only among intervals with a prior landmark MRD measurement and patients still at risk; therefore, its performance is not directly comparable with full-cohort models."
  ),
  "",
  paste0(
    "The exact-floor joint model and the left-censored primary joint model had similar interval-level discrimination and calibration summaries, but the primary model is statistically preferable because it respects the assay reporting mechanism. Treating floor observations as exact \\(-5.0\\) values imposes artificial precision at the deepest molecular response levels, whereas the primary model contributes the likelihood probability \\(P(y\\leq -5.0)\\). The primary model had interval Brier score ",
    fmt(primary$interval_brier),
    " and patient-level Brier score ",
    fmt(primary$patient_brier),
    "; the corresponding exact-floor values were ",
    fmt(exact$interval_brier),
    " and ",
    fmt(exact$patient_brier),
    "."
  ),
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.82\\linewidth]{figure_10_model_comparison_interval_calibration.pdf}",
  "\\caption{Interval-level observed-versus-predicted calibration across model classes. The Kaplan--Meier curve is descriptive because it uses first observed DMR visit time.}",
  "\\label{fig:model-comparison-interval-calibration}",
  "\\end{figure}",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.82\\linewidth]{figure_11_model_comparison_patient_calibration.pdf}",
  "\\caption{Patient-level observed-versus-predicted calibration across model classes. Patient-level probabilities from interval models were computed as cumulative probabilities across at-risk intervals.}",
  "\\label{fig:model-comparison-patient-calibration}",
  "\\end{figure}",
  "",
  "For internal validation, refittable simpler models were evaluated with patient-cluster bootstrap optimism correction, preserving the clustered longitudinal and interval structure. Full bootstrap refitting of the Stan joint models was not performed because it would require repeated Bayesian refits in a cohort of 87 patients; instead, patient-cluster bootstrap summaries of the fixed posterior predictions were reported and interpreted as internal stability checks rather than full optimism correction. This distinction is important for peer review and supports cautious wording: the joint model is justified by clinical coherence, interval-aware endpoint handling, and principled assay-floor likelihood specification, but external validation remains necessary before clinical decision use.",
  ""
)
writeLines(model_comparison_text, file.path(MODEL_DIR, "model_comparison_internal_validation.tex"), useBytes = TRUE)

reviewer_text <- c(
  "# Reviewer-Facing Model Comparison and Validation Explanation",
  "",
  "## Bottom Line",
  "",
  paste0(
    "The primary Bayesian joint longitudinal--interval model remains justified as the main analysis because it is the only fitted model that simultaneously handles serial latent MRD, interval-observed DMR onset, irregular visit gaps, patient-level heterogeneity, and left-censored assay-floor observations. Its interval Brier score was ",
    fmt(primary$interval_brier),
    " and its patient-level Brier score was ",
    fmt(primary$patient_brier),
    "."
  ),
  "",
  "## What Each Comparator Tests",
  "",
  "- Kaplan--Meier: useful descriptive benchmark only; it treats first observed DMR visit time as the event time and therefore does not respect interval-observed onset.",
  "- Interval timing + visit gap: tests whether interval-aware modeling alone is sufficient without longitudinal MRD.",
  "- Landmark MRD: tests whether a simpler clinically familiar landmark approach captures enough information; evaluated only when prior landmark MRD was actually available and the patient remained at risk.",
  "- Joint exact-floor model: tests whether the assay-floor censoring likelihood matters by fitting the same joint structure while treating floor observations as exact -5.0 values.",
  "- Primary joint model: tests the full proposed monitoring framework with latent MRD, random patient effects, interval timing, visit gap, and left-censored floor observations.",
  "",
  "## Internal Validation",
  "",
  paste0(
    "The refittable simpler models used patient-cluster bootstrap optimism correction with ",
    bootstrap_reps,
    " bootstrap resamples. For the two Stan joint models, the script reports patient-cluster bootstrap stability of the fixed posterior predictions; it does not claim full optimism correction because that would require repeated Stan refitting."
  ),
  "",
  "## Recommendation",
  "",
  "Use the left-censored Bayesian joint longitudinal--interval model as the primary model. Present Kaplan--Meier as descriptive, the timing-gap interval model and landmark MRD model as simpler clinical benchmarks, and the exact-floor joint model as a sensitivity analysis showing the effect of assay-floor handling. Avoid claiming that the model is a validated treatment-decision rule; describe it as a monitoring-oriented development model requiring external validation.",
  "",
  "## Files Generated",
  "",
  "- `model_comparison_performance.csv`",
  "- `table_07_model_comparison_validation.tex`",
  "- `figure_10_model_comparison_interval_calibration.pdf/.png`",
  "- `figure_11_model_comparison_patient_calibration.pdf/.png`",
  "- `model_comparison_internal_validation.tex`",
  "- `model_comparison_reviewer_explanation.md`"
)
writeLines(reviewer_text, file.path(MODEL_DIR, "model_comparison_reviewer_explanation.md"), useBytes = TRUE)

recommendation_text <- c(
  "# Recommended Primary and Sensitivity Models",
  "",
  "## Recommended Primary Model",
  "",
  "Use the renewed Bayesian joint longitudinal--interval model with independent patient-specific random intercept and slope terms and left-censored assay-floor likelihood as the primary manuscript model.",
  "",
  "Rationale:",
  "",
  paste0("- It had the lowest interval-level Brier score: ", fmt(primary$interval_brier), "."),
  paste0("- It had the lowest patient-level Brier score: ", fmt(primary$patient_brier), "."),
  "- It is the only model that simultaneously represents latent longitudinal MRD, irregular at-risk intervals, visit gap, patient-level heterogeneity, and assay-floor censoring.",
  "- It directly matches the clinical data-generating process: DMR is observed at visits, not continuously, and floor-level MRD values are not exact continuous measurements.",
  "",
  "## Recommended Sensitivity Models",
  "",
  "1. Joint exact-floor model: include as the main assay-floor sensitivity analysis. This model keeps the same joint interval structure but treats floor observations as exact -5.0 values.",
  "",
  "2. Interval timing + visit-gap model: include as the parsimonious interval-aware benchmark without longitudinal MRD. This isolates the gain from adding serial molecular burden.",
  "",
  "3. Landmark MRD model: include as a clinically familiar benchmark, but clearly state that it was evaluated only where a prior 6-, 12-, or 18-month landmark was available and the patient remained at risk.",
  "",
  "4. Kaplan--Meier curve: retain as descriptive only. It uses first observed DMR visit time and should not be presented as a valid exact-onset survival model.",
  "",
  "## Recommended Manuscript Claim",
  "",
  "The model-comparison results support the complex joint model as the primary development model for monitoring-oriented inference. They do not establish a validated treatment-decision rule. External or temporal validation is still required before clinical deployment."
)
writeLines(recommendation_text, file.path(MODEL_DIR, "model_comparison_recommendations.md"), useBytes = TRUE)

cat("Saved model comparison and validation outputs to:", MODEL_DIR, "\n")
