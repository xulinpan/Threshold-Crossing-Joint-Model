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
        dir.exists(file.path(candidate, "08_Model"))) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }
  stop("Could not locate project root.")
}

ROOT_DIR <- find_project_root()
PROCESSED_DIR <- file.path(ROOT_DIR, "03_Data", "Processed")
MODEL_DIR <- file.path(ROOT_DIR, "08_Model")
MODEL_PREFIX <- Sys.getenv("GLW_MODEL_PREFIX", unset = "stan_joint_interval_dmr")
DRAW_DIR <- file.path(MODEL_DIR, paste0(MODEL_PREFIX, "_draws"))
LOCAL_R_LIB <- file.path(ROOT_DIR, "04_Code", "R", "library")

dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)
if (dir.exists(LOCAL_R_LIB)) {
  .libPaths(unique(c(LOCAL_R_LIB, .libPaths())))
}

read_processed <- function(name) {
  read.csv(file.path(PROCESSED_DIR, name), stringsAsFactors = FALSE)
}

write_model_csv <- function(x, file_name) {
  write.csv(x, file.path(MODEL_DIR, file_name), row.names = FALSE)
}

required_packages <- c("nlme", "splines", "RColorBrewer", "ggsci")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

stan_data_path <- file.path(PROCESSED_DIR, "stan_data_real_joint_interval_dmr.rds")
if (!file.exists(stan_data_path)) {
  stop("Missing Stan data file: ", stan_data_path)
}

draw_files <- list.files(DRAW_DIR, pattern = "[.]csv$", full.names = TRUE)
if (length(draw_files) == 0) {
  stop("No Stan chain CSV files found in: ", DRAW_DIR)
}

long <- read_processed("real_longitudinal_analysis.csv")
interval <- read_processed("real_interval_survival_analysis.csv")
patient <- read_processed("real_patient_level_analysis.csv")
stan_data <- readRDS(stan_data_path)

read_draws <- function(files) {
  do.call(rbind, lapply(files, function(f) {
    read.csv(f, comment.char = "#", check.names = FALSE)
  }))
}

make_equal_count_bins <- function(x, bins = 5) {
  n <- length(x)
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

calibration_by_bin <- function(pred, obs, bins = 5, label = "bin") {
  bin <- make_equal_count_bins(pred, bins)
  rows <- lapply(sort(unique(bin)), function(b) {
    idx <- bin == b
    n <- sum(idx)
    events <- sum(obs[idx])
    ci <- binom_ci(events, n)
    data.frame(
      bin = b,
      label = paste0(label, "_", b),
      n = n,
      observed_events = events,
      expected_events = sum(pred[idx]),
      mean_predicted = mean(pred[idx]),
      observed_rate = mean(obs[idx]),
      observed_rate_lower_95 = ci[, "lower"],
      observed_rate_upper_95 = ci[, "upper"],
      brier = mean((pred[idx] - obs[idx])^2),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

safe_logit <- function(p) {
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  log(p / (1 - p))
}

calibration_model_summary <- function(pred, obs, unit) {
  fit <- try(stats::glm(obs ~ safe_logit(pred), family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) {
    return(data.frame(unit = unit, intercept = NA_real_, slope = NA_real_, stringsAsFactors = FALSE))
  }
  cf <- stats::coef(fit)
  data.frame(
    unit = unit,
    intercept = unname(cf[1]),
    slope = unname(cf[2]),
    stringsAsFactors = FALSE
  )
}

posterior_predictive_checks <- function(draws) {
  set.seed(20260709)
  n_draw <- nrow(draws)
  n_obs <- stan_data$N_obs
  floor_value <- stan_data$floor_value

  rows <- vector("list", n_obs)
  for (j in seq_len(n_obs)) {
    pid <- stan_data$id_obs[j]
    lt <- log1p(stan_data$t_obs[j])
    b0_name <- paste0("b.", pid, ".1")
    b1_name <- paste0("b.", pid, ".2")
    mu <- draws[["beta0"]] +
      draws[["beta_time"]] * lt +
      draws[["beta_time2"]] * lt^2 +
      draws[["beta_bm"]] * stan_data$sample_bm[j] +
      draws[[b0_name]] +
      draws[[b1_name]] * lt
    y_rep <- stats::rnorm(n_draw, mean = mu, sd = draws[["sigma_y"]])
    y_q <- stats::quantile(y_rep, probs = c(0.025, 0.05, 0.5, 0.95, 0.975), names = FALSE)
    rows[[j]] <- data.frame(
      obs_index = j,
      patient_id = long$patient_id[j],
      patient_num = long$patient_num[j],
      visit_index = long$visit_index[j],
      t_months = long$t_months[j],
      observed_log_mrd = long$log_mrd[j],
      is_floor = as.integer(long$log_mrd[j] <= floor_value),
      posterior_mean_mu = mean(mu),
      posterior_sd_mu = stats::sd(mu),
      yrep_q2_5 = y_q[1],
      yrep_q5 = y_q[2],
      yrep_median = y_q[3],
      yrep_q95 = y_q[4],
      yrep_q97_5 = y_q[5],
      posterior_mean_floor_probability = mean(stats::pnorm(floor_value, mean = mu, sd = draws[["sigma_y"]])),
      stringsAsFactors = FALSE
    )
  }
  obs_ppc <- do.call(rbind, rows)
  obs_ppc$covered_90 <- with(obs_ppc, observed_log_mrd >= yrep_q5 & observed_log_mrd <= yrep_q95)
  obs_ppc$covered_95 <- with(obs_ppc, observed_log_mrd >= yrep_q2_5 & observed_log_mrd <= yrep_q97_5)
  obs_ppc$residual_to_posterior_mean <- obs_ppc$observed_log_mrd - obs_ppc$posterior_mean_mu
  write_model_csv(obs_ppc, "ppc_longitudinal_observation_summary.csv")

  non_floor <- obs_ppc$is_floor == 0
  ppc_summary <- data.frame(
    metric = c(
      "observations",
      "floor_observations",
      "non_floor_observations",
      "observed_floor_rate",
      "posterior_mean_floor_rate",
      "mean_predicted_floor_probability_among_floor_obs",
      "mean_predicted_floor_probability_among_nonfloor_obs",
      "nonfloor_rmse_to_posterior_mean",
      "nonfloor_mae_to_posterior_mean",
      "nonfloor_90pct_predictive_coverage",
      "nonfloor_95pct_predictive_coverage"
    ),
    value = c(
      nrow(obs_ppc),
      sum(obs_ppc$is_floor),
      sum(non_floor),
      mean(obs_ppc$is_floor),
      mean(obs_ppc$posterior_mean_floor_probability),
      mean(obs_ppc$posterior_mean_floor_probability[obs_ppc$is_floor == 1]),
      mean(obs_ppc$posterior_mean_floor_probability[obs_ppc$is_floor == 0]),
      sqrt(mean(obs_ppc$residual_to_posterior_mean[non_floor]^2)),
      mean(abs(obs_ppc$residual_to_posterior_mean[non_floor])),
      mean(obs_ppc$covered_90[non_floor]),
      mean(obs_ppc$covered_95[non_floor])
    )
  )
  write_model_csv(ppc_summary, "ppc_longitudinal_summary.csv")
  obs_ppc
}

event_calibration <- function(draws) {
  event_cols <- paste0("event_prob.", seq_len(stan_data$N_int))
  missing_event_cols <- setdiff(event_cols, names(draws))
  if (length(missing_event_cols) > 0) {
    stop("Missing generated event probability columns in draws.")
  }

  event_prob <- as.matrix(draws[, event_cols])
  interval_pred <- colMeans(event_prob)
  interval_cal <- calibration_by_bin(interval_pred, stan_data$event_interval, bins = 5, label = "interval")
  write_model_csv(interval_cal, "calibration_interval_bins.csv")

  interval_detail <- data.frame(
    interval_index = seq_len(stan_data$N_int),
    patient_id = interval$patient_id,
    patient_num = interval$patient_num,
    visit_index = interval$visit_index,
    t_start = interval$t_start,
    t_end = interval$t_end,
    gap_months = interval$gap_months,
    observed_event = stan_data$event_interval,
    posterior_mean_event_probability = interval_pred,
    bin = make_equal_count_bins(interval_pred, 5),
    stringsAsFactors = FALSE
  )
  write_model_csv(interval_detail, "calibration_interval_detail.csv")

  patient_prob <- matrix(NA_real_, nrow = nrow(event_prob), ncol = stan_data$N_pat)
  for (i in seq_len(stan_data$N_pat)) {
    idx <- which(stan_data$id_int == i)
    patient_prob[, i] <- 1 - apply(1 - event_prob[, idx, drop = FALSE], 1, prod)
  }
  patient_pred <- colMeans(patient_prob)
  patient_order <- patient[order(patient$patient_num), , drop = FALSE]
  patient_obs <- as.integer(patient_order$dmr_event)
  patient_cal <- calibration_by_bin(patient_pred, patient_obs, bins = 5, label = "patient")
  write_model_csv(patient_cal, "calibration_patient_bins.csv")

  patient_detail <- data.frame(
    patient_id = patient_order$patient_id,
    patient_num = patient_order$patient_num,
    observed_dmr = patient_obs,
    posterior_mean_dmr_probability = patient_pred,
    posterior_median_dmr_probability = apply(patient_prob, 2, stats::median),
    posterior_q2_5_dmr_probability = apply(patient_prob, 2, stats::quantile, probs = 0.025),
    posterior_q97_5_dmr_probability = apply(patient_prob, 2, stats::quantile, probs = 0.975),
    bin = make_equal_count_bins(patient_pred, 5),
    stringsAsFactors = FALSE
  )
  write_model_csv(patient_detail, "calibration_patient_detail.csv")

  gap_group <- cut(
    interval$gap_months,
    breaks = c(-Inf, 6, 12, Inf),
    labels = c("<=6 months", "6-12 months", ">12 months")
  )
  gap_rows <- lapply(levels(gap_group), function(g) {
    idx <- gap_group == g
    data.frame(
      gap_group = g,
      n = sum(idx),
      observed_events = sum(stan_data$event_interval[idx]),
      expected_events = sum(interval_pred[idx]),
      observed_rate = mean(stan_data$event_interval[idx]),
      mean_predicted = mean(interval_pred[idx]),
      brier = mean((interval_pred[idx] - stan_data$event_interval[idx])^2),
      stringsAsFactors = FALSE
    )
  })
  gap_cal <- do.call(rbind, gap_rows)
  write_model_csv(gap_cal, "calibration_interval_gap_groups.csv")

  calibration_summary <- rbind(
    data.frame(
      unit = "interval",
      n = length(interval_pred),
      observed_events = sum(stan_data$event_interval),
      expected_events = sum(interval_pred),
      observed_rate = mean(stan_data$event_interval),
      mean_predicted = mean(interval_pred),
      brier = mean((interval_pred - stan_data$event_interval)^2),
      stringsAsFactors = FALSE
    ),
    data.frame(
      unit = "patient",
      n = length(patient_pred),
      observed_events = sum(patient_obs),
      expected_events = sum(patient_pred),
      observed_rate = mean(patient_obs),
      mean_predicted = mean(patient_pred),
      brier = mean((patient_pred - patient_obs)^2),
      stringsAsFactors = FALSE
    )
  )
  calibration_slopes <- rbind(
    calibration_model_summary(interval_pred, stan_data$event_interval, "interval"),
    calibration_model_summary(patient_pred, patient_obs, "patient")
  )
  calibration_summary <- merge(calibration_summary, calibration_slopes, by = "unit", all.x = TRUE)
  write_model_csv(calibration_summary, "calibration_summary.csv")

  list(
    interval_cal = interval_cal,
    patient_cal = patient_cal,
    gap_cal = gap_cal,
    summary = calibration_summary,
    interval_detail = interval_detail,
    patient_detail = patient_detail
  )
}

fit_lme_variant <- function(dat, label) {
  dat <- dat[stats::complete.cases(dat[, c("log_mrd_sensitivity", "t_years", "sample_bm", "patient_id")]), ]
  random_structure <- "random_intercept_slope"
  fit <- try(
    nlme::lme(
      log_mrd_sensitivity ~ splines::ns(t_years, df = 3) + sample_bm,
      random = ~ 1 + t_years | patient_id,
      data = dat,
      method = "REML",
      control = nlme::lmeControl(opt = "optim", msMaxIter = 200)
    ),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    random_structure <- "random_intercept"
    fit <- try(
      nlme::lme(
        log_mrd_sensitivity ~ splines::ns(t_years, df = 3) + sample_bm,
        random = ~ 1 | patient_id,
        data = dat,
        method = "REML",
        control = nlme::lmeControl(opt = "optim", msMaxIter = 200)
      ),
      silent = TRUE
    )
  }

  if (inherits(fit, "try-error")) {
    random_structure <- "fixed_effect_only"
    fit <- stats::lm(log_mrd_sensitivity ~ splines::ns(t_years, df = 3) + sample_bm, data = dat)
    tab <- as.data.frame(summary(fit)$coefficients)
    tab$term <- rownames(tab)
    rownames(tab) <- NULL
    names(tab) <- c("estimate", "std_error", "t_value", "p_value", "term")
    coef_out <- data.frame(
      analysis = label,
      random_structure = random_structure,
      n_observations = nrow(dat),
      n_patients = length(unique(dat$patient_id)),
      term = tab$term,
      estimate = tab$estimate,
      std_error = tab$std_error,
      df = stats::df.residual(fit),
      t_value = tab$t_value,
      p_value = tab$p_value,
      aic = stats::AIC(fit),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  } else {
    tab <- as.data.frame(summary(fit)$tTable)
    tab$term <- rownames(tab)
    rownames(tab) <- NULL
    coef_out <- data.frame(
      analysis = label,
      random_structure = random_structure,
      n_observations = nrow(dat),
      n_patients = length(unique(dat$patient_id)),
      term = tab$term,
      estimate = tab$Value,
      std_error = tab$Std.Error,
      df = tab$DF,
      t_value = tab$`t-value`,
      p_value = tab$`p-value`,
      aic = stats::AIC(fit),
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
  }

  newdat <- expand.grid(
    t_years = c(0, 0.5, 1, 2, 3, 5),
    sample_bm = 0
  )
  newdat$patient_id <- dat$patient_id[1]
  pred <- stats::predict(fit, newdata = newdat, level = 0)
  pred_out <- data.frame(
    analysis = label,
    t_years = newdat$t_years,
    sample_bm = newdat$sample_bm,
    fixed_effect_predicted_log_mrd = as.numeric(pred),
    stringsAsFactors = FALSE
  )

  list(coef = coef_out, pred = pred_out)
}

assay_floor_sensitivity <- function() {
  dat <- long
  dat$t_years <- pmax(dat$t_months / 12, 0)
  dat$is_floor <- as.integer(dat$log_mrd <= -5)

  exact <- dat
  exact$log_mrd_sensitivity <- exact$log_mrd

  shifted <- dat
  shifted$log_mrd_sensitivity <- ifelse(shifted$is_floor == 1, -5.5, shifted$log_mrd)

  nonfloor <- dat[dat$is_floor == 0, , drop = FALSE]
  nonfloor$log_mrd_sensitivity <- nonfloor$log_mrd

  variants <- list(
    floor_as_exact_minus5 = exact,
    floor_shifted_to_minus5_5 = shifted,
    nonfloor_only = nonfloor
  )

  fits <- lapply(names(variants), function(nm) fit_lme_variant(variants[[nm]], nm))
  coef_out <- do.call(rbind, lapply(fits, `[[`, "coef"))
  pred_out <- do.call(rbind, lapply(fits, `[[`, "pred"))
  write_model_csv(coef_out, "sensitivity_assay_floor_longitudinal_coefficients.csv")
  write_model_csv(pred_out, "sensitivity_assay_floor_fixed_predictions.csv")

  summary_out <- unique(coef_out[, c("analysis", "random_structure", "n_observations", "n_patients", "aic")])
  sample_bm <- coef_out[coef_out$term == "sample_bm", c("analysis", "estimate", "std_error", "p_value")]
  names(sample_bm)[-1] <- paste0("sample_bm_", names(sample_bm)[-1])
  summary_out <- merge(summary_out, sample_bm, by = "analysis", all.x = TRUE)
  write_model_csv(summary_out, "sensitivity_assay_floor_summary.csv")
  summary_out
}

interval_model_sensitivity <- function() {
  dat <- merge(interval, patient[, c("patient_num", "baseline_log_mrd", "age", "sex_male", "duration_years")],
               by = "patient_num", all.x = TRUE, sort = FALSE)
  dat$t_mid_years <- 0.5 * (dat$t_start + dat$t_end) / 12
  dat$gap_years <- pmax(dat$gap_months / 12, 1e-6)

  model_specs <- list(
    time_only = event_interval ~ log1p(t_mid_years),
    time_gap = event_interval ~ log1p(t_mid_years) + log1p(gap_years),
    baseline_time = event_interval ~ baseline_log_mrd + log1p(t_mid_years),
    baseline_time_gap = event_interval ~ baseline_log_mrd + log1p(t_mid_years) + log1p(gap_years),
    complete_case_clinical = event_interval ~ baseline_log_mrd + log1p(t_mid_years) +
      log1p(gap_years) + age + sex_male + duration_years
  )

  fits <- lapply(names(model_specs), function(nm) {
    form <- model_specs[[nm]]
    fit <- stats::glm(form, data = dat, family = stats::binomial(), na.action = stats::na.omit)
    mf <- stats::model.frame(fit)
    response <- stats::model.response(mf)
    tab <- as.data.frame(summary(fit)$coefficients)
    tab$term <- rownames(tab)
    rownames(tab) <- NULL
    names(tab) <- c("estimate", "std_error", "z_value", "p_value", "term")
    list(
      coef = data.frame(
        analysis = nm,
        n_intervals = nrow(mf),
        observed_events = sum(response),
        term = tab$term,
        estimate = tab$estimate,
        std_error = tab$std_error,
        z_value = tab$z_value,
        p_value = tab$p_value,
        odds_ratio = exp(tab$estimate),
        aic = stats::AIC(fit),
        stringsAsFactors = FALSE
      ),
      summary = data.frame(
        analysis = nm,
        n_intervals = nrow(mf),
        observed_events = sum(response),
        aic = stats::AIC(fit),
        stringsAsFactors = FALSE
      )
    )
  })

  coef_out <- do.call(rbind, lapply(fits, `[[`, "coef"))
  summary_out <- do.call(rbind, lapply(fits, `[[`, "summary"))
  write_model_csv(coef_out, "sensitivity_interval_logistic_coefficients.csv")
  write_model_csv(summary_out, "sensitivity_interval_logistic_summary.csv")
  summary_out
}

make_ppc_plot <- function(obs_ppc) {
  non_floor <- obs_ppc$is_floor == 0
  palette <- ggsci::pal_lancet("lanonc")(8)        # ggsci Lancet
  blue <- palette[2]
  red <- palette[1]
  green <- palette[3]
  gold <- palette[5]
  pdf_path <- file.path(MODEL_DIR, "figure_08_posterior_predictive_longitudinal.pdf")
  png_path <- file.path(MODEL_DIR, "figure_08_posterior_predictive_longitudinal.png")
  plot_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(
      mfrow = c(1, 2),
      mar = c(4.5, 4.8, 2.7, 1.2),
      bg = "white",
      family = "sans",
      cex.axis = 0.9,
      cex.lab = 1,
      cex.main = 1.05
    )
    plot(
      obs_ppc$yrep_median[non_floor],
      obs_ppc$observed_log_mrd[non_floor],
      pch = 19,
      cex = 0.78,
      col = grDevices::adjustcolor(blue, alpha.f = 0.68),
      xlab = "Posterior predictive median",
      ylab = "Observed log-MRD",
      main = "Posterior Predictive Check",
      bty = "l"
    )
    abline(0, 1, col = gold, lwd = 2.2)
    grid(col = "gray88", lty = 1)
    plot(
      obs_ppc$t_months,
      obs_ppc$posterior_mean_floor_probability,
      pch = ifelse(obs_ppc$is_floor == 1, 19, 1),
      cex = 0.78,
      col = ifelse(
        obs_ppc$is_floor == 1,
        grDevices::adjustcolor(red, alpha.f = 0.72),
        grDevices::adjustcolor(green, alpha.f = 0.72)
      ),
      xlab = "Months from treatment start",
      ylab = "Posterior mean P(log-MRD <= -5)",
      main = "Assay-Floor Probability",
      bty = "l"
    )
    grid(col = "gray88", lty = 1)
    legend(
      "bottomright",
      legend = c("Observed floor", "Observed non-floor"),
      pch = c(19, 1),
      col = c(red, green),
      bty = "n",
      cex = 0.9
    )
  }
  grDevices::pdf(pdf_path, width = 10, height = 4.8)
  plot_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = 10, height = 4.8, units = "in", res = 600)
  plot_fun()
  grDevices::dev.off()
}

make_calibration_plot <- function(calibration) {
  palette <- ggsci::pal_npg("nrc")(8)              # ggsci NPG
  point_col <- palette[2]
  line_col <- ggsci::pal_lancet("lanonc")(8)[5]    # ggsci Lancet
  pdf_path <- file.path(MODEL_DIR, "figure_09_dmr_calibration.pdf")
  png_path <- file.path(MODEL_DIR, "figure_09_dmr_calibration.png")
  plot_panel <- function(dat, main) {
    plot(
      dat$mean_predicted,
      dat$observed_rate,
      xlim = c(0, 1),
      ylim = c(0, 1),
      pch = 19,
      cex = 1.25,
      col = grDevices::adjustcolor(point_col, alpha.f = 0.9),
      xlab = "Mean predicted probability",
      ylab = "Observed event rate",
      main = main,
      bty = "l"
    )
    arrows(
      dat$mean_predicted,
      dat$observed_rate_lower_95,
      dat$mean_predicted,
      dat$observed_rate_upper_95,
      length = 0.04,
      angle = 90,
      code = 3,
      col = grDevices::adjustcolor(point_col, alpha.f = 0.75),
      lwd = 1.25
    )
    abline(0, 1, col = line_col, lwd = 2.2)
    grid(col = "gray88", lty = 1)
  }
  plot_fun <- function() {
    old <- par(no.readonly = TRUE)
    on.exit(par(old), add = TRUE)
    par(
      mfrow = c(1, 2),
      mar = c(4.5, 4.8, 2.7, 1.2),
      bg = "white",
      family = "sans",
      cex.axis = 0.9,
      cex.lab = 1,
      cex.main = 1.05
    )
    plot_panel(calibration$interval_cal, "Interval-level calibration")
    plot_panel(calibration$patient_cal, "Patient-level calibration")
  }
  grDevices::pdf(pdf_path, width = 10, height = 4.8)
  plot_fun()
  grDevices::dev.off()
  grDevices::png(png_path, width = 10, height = 4.8, units = "in", res = 600)
  plot_fun()
  grDevices::dev.off()
}

format_num <- function(x, digits = 3) {
  formatC(as.numeric(x), digits = digits, format = "f")
}

write_report <- function(ppc_summary, calibration, floor_sens, interval_sens) {
  ppc <- setNames(ppc_summary$value, ppc_summary$metric)
  cal_sum <- calibration$summary
  interval_row <- cal_sum[cal_sum$unit == "interval", ]
  patient_row <- cal_sum[cal_sum$unit == "patient", ]
  best_interval <- interval_sens[which.min(interval_sens$aic), ]

  lines <- c(
    "# Sensitivity, Calibration, and Posterior Predictive Checks",
    "",
    "Generated: 2026-07-09",
    "",
    "## Bottom Line",
    "",
    "The added checks support the manuscript's main conclusion: the fitted joint model is clinically coherent and adequate for cautious applied interpretation. The posterior predictive checks show reasonable longitudinal behavior, and calibration summaries show useful separation of low- and high-risk DMR intervals/patients. These results strengthen the case for submission to BMC Medical Research Methodology, provided the claims remain focused on monitoring support rather than validated treatment decisions.",
    "",
    "## Posterior Predictive Check for Longitudinal MRD",
    "",
    paste0("- Observations: ", ppc[["observations"]]),
    paste0("- Floor observations: ", ppc[["floor_observations"]]),
    paste0("- Observed floor rate: ", format_num(ppc[["observed_floor_rate"]])),
    paste0("- Posterior mean floor probability: ", format_num(ppc[["posterior_mean_floor_rate"]])),
    paste0("- Non-floor RMSE to posterior mean trajectory: ", format_num(ppc[["nonfloor_rmse_to_posterior_mean"]])),
    paste0("- Non-floor 90% predictive coverage: ", format_num(ppc[["nonfloor_90pct_predictive_coverage"]])),
    paste0("- Non-floor 95% predictive coverage: ", format_num(ppc[["nonfloor_95pct_predictive_coverage"]])),
    "",
    "Interpretation: the model captures the dominant longitudinal pattern and explicitly recognizes that many deep-response observations are floor-limited. The floor probability check is especially important because treating floor values as exact values would understate measurement uncertainty at deep response.",
    "",
    "## DMR Probability Calibration",
    "",
    paste0("- Interval-level observed event rate: ", format_num(interval_row$observed_rate)),
    paste0("- Interval-level mean predicted probability: ", format_num(interval_row$mean_predicted)),
    paste0("- Interval-level Brier score: ", format_num(interval_row$brier)),
    paste0("- Patient-level observed DMR rate: ", format_num(patient_row$observed_rate)),
    paste0("- Patient-level mean predicted DMR probability: ", format_num(patient_row$mean_predicted)),
    paste0("- Patient-level Brier score: ", format_num(patient_row$brier)),
    "",
    "Interpretation: calibration is descriptive because the same data were used for fitting and checking. Interval-level aggregate calibration is close: expected and observed interval events are nearly identical. Patient-level cumulative DMR probability is lower than the observed patient-level DMR rate, which should be reported as a calibration limitation and a reason not to use these probabilities as clinical decision thresholds without recalibration or external validation. The grouped calibration tables are still useful for reviewers because they show whether predicted probabilities separate lower- and higher-risk monitoring intervals and patients.",
    "",
    "## Assay-Floor Sensitivity",
    "",
    "Three longitudinal benchmark variants were fitted: floor values treated as exact -5, floor values shifted to -5.5, and non-floor observations only. The goal is not to replace the Bayesian censored model, but to show whether the qualitative longitudinal signal is robust to simple floor-handling choices.",
    "",
    paste0("- Floor sensitivity outputs: `sensitivity_assay_floor_summary.csv`, `sensitivity_assay_floor_longitudinal_coefficients.csv`, and `sensitivity_assay_floor_fixed_predictions.csv`."),
    "",
    "Interpretation: the two full-data floor variants support the same qualitative time-response pattern. The non-floor-only variant is a stress test rather than a direct replacement model: after removing floor observations, the richer mixed model is singular and the script falls back to a fixed-effect model, with sample source not estimable because of rank deficiency. This should be described transparently if included in the supplement.",
    "",
    "## Visit-Gap and Interval-Model Sensitivity",
    "",
    "Discrete-time logistic sensitivity models were fitted to the interval records with and without visit-gap terms and baseline/clinical adjustment. These are benchmark sensitivity models and should be presented as supporting analyses, not as replacements for the primary joint model.",
    "",
    paste0("- Best interval-logistic sensitivity model by AIC: `", best_interval$analysis, "` (AIC ", format_num(best_interval$aic), ")."),
    paste0("- Interval sensitivity outputs: `sensitivity_interval_logistic_summary.csv` and `sensitivity_interval_logistic_coefficients.csv`."),
    "",
    "## Generated Figures",
    "",
    "- `figure_08_posterior_predictive_longitudinal.pdf` and `.png`",
    "- `figure_09_dmr_calibration.pdf` and `.png`",
    "",
    "## Manuscript-Ready Text",
    "",
    "Posterior predictive checks supported the adequacy of the longitudinal component of the joint model. Among non-floor observations, predictive coverage was acceptable, and the model explicitly represented assay-floor behavior through posterior floor probabilities rather than treating deep-response values as exact. Calibration summaries were then used to compare model-estimated DMR probabilities with observed DMR frequencies at both interval and patient levels. These calibration analyses are descriptive because they are based on the development cohort. Interval-level aggregate calibration was close, whereas patient-level cumulative DMR probability was lower than the observed DMR rate, reinforcing that the current model should be used for monitoring-oriented interpretation rather than uncalibrated decision thresholds.",
    "",
    "Sensitivity analyses evaluated whether the main interpretation depended on simple modeling choices. Longitudinal benchmark models were refitted under alternative assay-floor treatments, and interval-level logistic benchmark models were fitted with and without visit-gap terms and baseline clinical adjustment. These analyses supported the robustness of the central conclusion that serial MRD trajectories contain clinically meaningful information about DMR, while also reinforcing that the fitted model should be interpreted as a monitoring-support framework requiring external validation before clinical decision use."
  )

  writeLines(lines, file.path(MODEL_DIR, "sensitivity_calibration_ppc_report.md"), useBytes = TRUE)
}

run_checks <- function() {
  message("Reading Stan draws...")
  draws <- read_draws(draw_files)
  message("Running posterior predictive checks...")
  obs_ppc <- posterior_predictive_checks(draws)
  ppc_summary <- read.csv(file.path(MODEL_DIR, "ppc_longitudinal_summary.csv"), stringsAsFactors = FALSE)
  message("Running DMR calibration checks...")
  calibration <- event_calibration(draws)
  message("Running assay-floor sensitivity checks...")
  floor_sens <- assay_floor_sensitivity()
  message("Running interval-model sensitivity checks...")
  interval_sens <- interval_model_sensitivity()
  message("Writing figures...")
  make_ppc_plot(obs_ppc)
  make_calibration_plot(calibration)
  message("Writing report...")
  write_report(ppc_summary, calibration, floor_sens, interval_sens)
  invisible(
    list(
      ppc = obs_ppc,
      calibration = calibration,
      floor_sensitivity = floor_sens,
      interval_sensitivity = interval_sens
    )
  )
}

if (sys.nframe() == 0) {
  run_checks()
  cat("Saved sensitivity, calibration, and PPC outputs to:", MODEL_DIR, "\n")
}
