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
    "04_Code/R/12_dynamic_prediction_recalibration.R"
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
DRAW_DIR <- file.path(MODEL_DIR, "stan_joint_interval_dmr_independent_renewed_draws")
LOCAL_R_LIB <- file.path(ROOT_DIR, "04_Code", "R", "library")
if (dir.exists(LOCAL_R_LIB)) {
  .libPaths(unique(c(LOCAL_R_LIB, .libPaths())))
}

required_packages <- c("RColorBrewer", "ggsci")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}

dir.create(MODEL_DIR, showWarnings = FALSE, recursive = TRUE)

safe_prob <- function(p, eps = 1e-6) {
  pmin(pmax(as.numeric(p), eps), 1 - eps)
}

safe_logit <- function(p) {
  p <- safe_prob(p)
  log(p / (1 - p))
}

inv_logit <- function(x) {
  1 / (1 + exp(-x))
}

fmt <- function(x, digits = 3) {
  ifelse(is.na(x), "--", formatC(as.numeric(x), digits = digits, format = "f"))
}

brier_score <- function(obs, pred) {
  mean((safe_prob(pred) - obs)^2)
}

calibration_intercept <- function(obs, pred) {
  obs <- as.integer(obs)
  pred <- safe_prob(pred)
  if (length(unique(obs)) < 2 || length(obs) < 10) return(NA_real_)
  fit <- try(suppressWarnings(stats::glm(obs ~ 1 + offset(safe_logit(pred)), family = stats::binomial())), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  out <- unname(stats::coef(fit)[1])
  ifelse(is.finite(out), out, NA_real_)
}

calibration_slope <- function(obs, pred) {
  obs <- as.integer(obs)
  pred <- safe_prob(pred)
  if (length(unique(obs)) < 2 || length(obs) < 10 || stats::sd(pred) == 0) return(NA_real_)
  fit <- try(suppressWarnings(stats::glm(obs ~ safe_logit(pred), family = stats::binomial())), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  out <- unname(stats::coef(fit)[2])
  ifelse(is.finite(out), out, NA_real_)
}

metric_set <- function(obs, pred) {
  data.frame(
    n = length(obs),
    events = sum(obs),
    observed_rate = mean(obs),
    mean_predicted = mean(safe_prob(pred)),
    brier = brier_score(obs, pred),
    calibration_intercept = calibration_intercept(obs, pred),
    calibration_slope = calibration_slope(obs, pred),
    stringsAsFactors = FALSE
  )
}

read_processed <- function(name) {
  read.csv(file.path(PROCESSED_DIR, name), stringsAsFactors = FALSE)
}

read_draw_subset <- function(draw_dir, n_pat, max_draws = Inf) {
  draw_files <- list.files(draw_dir, pattern = "[.]csv$", full.names = TRUE)
  if (length(draw_files) == 0) stop("No draw CSV files found in: ", draw_dir)

  needed <- c(
    "beta0", "beta_time", "beta_time2",
    "gamma0", "gamma_time", "gamma_gap", "alpha_mrd",
    paste0("b.", seq_len(n_pat), ".1"),
    paste0("b.", seq_len(n_pat), ".2")
  )

  parts <- lapply(draw_files, function(file) {
    header <- names(read.csv(
      file,
      comment.char = "#",
      nrows = 0,
      check.names = FALSE
    ))
    missing <- setdiff(needed, header)
    if (length(missing) > 0) {
      stop("Draw file is missing needed columns: ", paste(missing, collapse = ", "))
    }
    col_classes <- rep("NULL", length(header))
    names(col_classes) <- header
    col_classes[needed] <- "numeric"
    read.csv(
      file,
      comment.char = "#",
      check.names = FALSE,
      colClasses = col_classes
    )
  })
  draws <- do.call(rbind, parts)
  if (is.finite(max_draws) && nrow(draws) > max_draws) {
    set.seed(20260709)
    draws <- draws[sort(sample(seq_len(nrow(draws)), max_draws)), , drop = FALSE]
  }
  draws
}

make_landmark_dataset <- function(patient, long, landmarks_months, horizons_months) {
  rows <- list()
  idx <- 0L
  for (s in landmarks_months) {
    history_counts <- aggregate(
      t_months ~ patient_num,
      long[long$t_months <= s, , drop = FALSE],
      length
    )
    names(history_counts)[2] <- "n_history_mrd"
    last_history <- do.call(rbind, lapply(split(long[long$t_months <= s, , drop = FALSE], long$patient_num[long$t_months <= s]), function(dat) {
      dat <- dat[order(dat$t_months), , drop = FALSE]
      tail(dat, 1)
    }))
    last_history <- last_history[, c("patient_num", "t_months", "log_mrd", "sample_bm"), drop = FALSE]
    names(last_history) <- c("patient_num", "last_mrd_month", "last_log_mrd", "last_sample_bm")

    for (H in horizons_months) {
      tmp <- merge(patient, history_counts, by = "patient_num", all.x = TRUE, sort = FALSE)
      tmp <- merge(tmp, last_history, by = "patient_num", all.x = TRUE, sort = FALSE)
      tmp$n_history_mrd[is.na(tmp$n_history_mrd)] <- 0

      event_by_horizon <- tmp$dmr_event == 1 & tmp$time_to_dmr_or_censor <= s + H
      outcome_observable <- event_by_horizon | tmp$followup_months >= s + H
      eligible <- tmp$n_history_mrd > 0 &
        tmp$time_to_dmr_or_censor > s &
        outcome_observable

      tmp <- tmp[eligible, , drop = FALSE]
      if (nrow(tmp) == 0) next

      idx <- idx + 1L
      tmp$landmark_month <- s
      tmp$horizon_month <- H
      tmp$prediction_time_month <- s + H
      tmp$outcome_dmr_by_horizon <- as.integer(
        tmp$dmr_event == 1 & tmp$time_to_dmr_or_censor <= s + H
      )
      tmp$support_label <- ifelse(nrow(tmp) >= 20 && sum(tmp$outcome_dmr_by_horizon) >= 5,
        "supported",
        "sparse"
      )
      rows[[idx]] <- tmp
    }
  }
  out <- do.call(rbind, rows)
  out[order(out$horizon_month, out$landmark_month, out$patient_num), , drop = FALSE]
}

predict_dynamic_probabilities <- function(landmark_data, draws) {
  pred <- numeric(nrow(landmark_data))
  pred_q05 <- numeric(nrow(landmark_data))
  pred_q95 <- numeric(nrow(landmark_data))

  for (r in seq_len(nrow(landmark_data))) {
    patient_idx <- landmark_data$patient_num[r]
    horizon_years <- landmark_data$horizon_month[r] / 12
    t_mid_years <- (landmark_data$landmark_month[r] + 0.5 * landmark_data$horizon_month[r]) / 12
    lt <- log1p(t_mid_years)
    b0 <- draws[[paste0("b.", patient_idx, ".1")]]
    b1 <- draws[[paste0("b.", patient_idx, ".2")]]
    mu_mid <- draws$beta0 + draws$beta_time * lt + draws$beta_time2 * lt^2 +
      b0 + b1 * lt
    eta <- draws$gamma0 + draws$gamma_time * lt +
      draws$gamma_gap * log1p(horizon_years) +
      draws$alpha_mrd * mu_mid
    p <- 1 - exp(-exp(eta) * horizon_years)
    p <- safe_prob(p)
    pred[r] <- mean(p)
    pred_q05[r] <- stats::quantile(p, 0.05, names = FALSE)
    pred_q95[r] <- stats::quantile(p, 0.95, names = FALSE)
  }

  landmark_data$pred_original <- safe_prob(pred)
  landmark_data$pred_original_q05 <- safe_prob(pred_q05)
  landmark_data$pred_original_q95 <- safe_prob(pred_q95)
  landmark_data
}

fit_recalibration <- function(dat) {
  fit <- try(
    suppressWarnings(stats::glm(
      outcome_dmr_by_horizon ~ safe_logit(pred_original),
      data = dat,
      family = stats::binomial()
    )),
    silent = TRUE
  )
  if (inherits(fit, "try-error")) {
    return(c(a = NA_real_, b = NA_real_))
  }
  co <- stats::coef(fit)
  c(a = unname(co[1]), b = unname(co[2]))
}

apply_recalibration <- function(p, a, b) {
  if (!is.finite(a) || !is.finite(b)) return(rep(NA_real_, length(p)))
  safe_prob(inv_logit(a + b * safe_logit(p)))
}

bootstrap_recalibration_optimism <- function(dat, B = 200) {
  patient_ids <- unique(dat$patient_num)
  if (length(patient_ids) < 10 || length(unique(dat$outcome_dmr_by_horizon)) < 2) {
    return(data.frame(
      bootstrap_reps = 0,
      brier_optimism = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  out <- numeric(0)
  set.seed(20260709)
  for (b in seq_len(B)) {
    sampled_ids <- sample(patient_ids, length(patient_ids), replace = TRUE)
    boot_rows <- do.call(rbind, lapply(seq_along(sampled_ids), function(j) {
      x <- dat[dat$patient_num == sampled_ids[j], , drop = FALSE]
      x$bootstrap_cluster <- j
      x
    }))
    if (length(unique(boot_rows$outcome_dmr_by_horizon)) < 2) next
    co <- fit_recalibration(boot_rows)
    if (!all(is.finite(co))) next
    p_boot <- apply_recalibration(boot_rows$pred_original, co[["a"]], co[["b"]])
    p_test <- apply_recalibration(dat$pred_original, co[["a"]], co[["b"]])
    apparent <- brier_score(boot_rows$outcome_dmr_by_horizon, p_boot)
    test <- brier_score(dat$outcome_dmr_by_horizon, p_test)
    out <- c(out, test - apparent)
  }

  data.frame(
    bootstrap_reps = length(out),
    brier_optimism = ifelse(length(out) > 0, mean(out), NA_real_),
    stringsAsFactors = FALSE
  )
}

calibration_bins <- function(dat, pred_col, model_label, n_bins = 4) {
  p <- safe_prob(dat[[pred_col]])
  q <- unique(stats::quantile(p, probs = seq(0, 1, length.out = n_bins + 1), na.rm = TRUE))
  if (length(q) < 3) {
    dat$calibration_bin <- "All"
  } else {
    q[1] <- min(q[1], min(p, na.rm = TRUE)) - 1e-8
    q[length(q)] <- max(q[length(q)], max(p, na.rm = TRUE)) + 1e-8
    dat$calibration_bin <- cut(p, breaks = q, include.lowest = TRUE, labels = FALSE)
  }
  parts <- split(dat, dat$calibration_bin)
  out <- do.call(rbind, lapply(names(parts), function(bin) {
    x <- parts[[bin]]
    data.frame(
      model = model_label,
      horizon_month = unique(x$horizon_month)[1],
      bin = bin,
      n = nrow(x),
      observed_rate = mean(x$outcome_dmr_by_horizon),
      mean_predicted = mean(safe_prob(x[[pred_col]])),
      min_predicted = min(safe_prob(x[[pred_col]])),
      max_predicted = max(safe_prob(x[[pred_col]])),
      stringsAsFactors = FALSE
    )
  }))
  row.names(out) <- NULL
  out
}

decision_curve <- function(dat, pred_col, model_label, thresholds = seq(0.05, 0.80, by = 0.05)) {
  y <- dat$outcome_dmr_by_horizon
  p <- safe_prob(dat[[pred_col]])
  prevalence <- mean(y)
  do.call(rbind, lapply(thresholds, function(pt) {
    treat <- p >= pt
    tp <- sum(treat & y == 1)
    fp <- sum(treat & y == 0)
    n <- length(y)
    nb_model <- tp / n - fp / n * pt / (1 - pt)
    nb_all <- prevalence - (1 - prevalence) * pt / (1 - pt)
    data.frame(
      model = model_label,
      horizon_month = unique(dat$horizon_month)[1],
      threshold = pt,
      net_benefit = nb_model,
      treat_all_net_benefit = nb_all,
      treat_none_net_benefit = 0,
      stringsAsFactors = FALSE
    )
  }))
}

write_validation_table_tex <- function(tab, path) {
  shown <- tab[tab$validation_level == "pooled_by_horizon", , drop = FALSE]
  shown <- shown[order(shown$horizon_month, shown$model), , drop = FALSE]
  rows <- vapply(seq_len(nrow(shown)), function(i) {
    paste(
      paste0(shown$horizon_month[i], " months"),
      shown$model[i],
      shown$n[i],
      shown$events[i],
      fmt(shown$observed_rate[i], 3),
      fmt(shown$mean_predicted[i], 3),
      fmt(shown$brier[i], 3),
      fmt(shown$optimism_corrected_brier[i], 3),
      fmt(shown$calibration_intercept[i], 2),
      fmt(shown$calibration_slope[i], 2),
      sep = " & "
    )
  }, character(1))

  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    "\\caption{Development-cohort dynamic prediction and recalibration performance. Predictions are landmark-level probabilities of documented DMR by the stated horizon among patients who were DMR-free at the landmark and had at least one prior MRD measurement. Optimism correction applies to the recalibration layer only; the Bayesian joint model itself was not refitted in bootstrap samples.}",
    "\\label{tab:dynamic-prediction-validation}",
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{llrrrrrrrr}",
    "\\toprule",
    "Horizon & Model & n & Events & Observed & Mean predicted & Brier & Opt.-corrected Brier & Cal. intercept & Cal. slope \\\\",
    "\\midrule",
    paste0(rows, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "}%",
    "\\end{table}",
    ""
  )
  writeLines(lines, path, useBytes = TRUE)
}

plot_calibration <- function(bins, path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(path, width = 7.4, height = 4.8, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = 7.4, height = 4.8, units = "in", res = 600)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  horizons <- sort(unique(bins$horizon_month))
  cols <- ggsci::pal_lancet("lanonc")(max(3, length(horizons)))[seq_along(horizons)]  # ggsci Lancet
  names(cols) <- horizons
  pchs <- c(16, 17)
  names(pchs) <- c("Original posterior", "Recalibrated")

  graphics::par(mfrow = c(1, 2), mar = c(4.2, 4.2, 2.2, 0.7), mgp = c(2.5, 0.8, 0), las = 1)
  for (model_name in c("Original posterior", "Recalibrated")) {
    graphics::plot(
      NA,
      xlim = c(0, 1),
      ylim = c(0, 1),
      xlab = "Mean predicted probability",
      ylab = "Observed DMR proportion",
      main = model_name,
      bty = "l"
    )
    graphics::abline(0, 1, col = "gray55", lty = 2, lwd = 1.2)
    for (h in horizons) {
      dat <- bins[bins$model == model_name & bins$horizon_month == h, , drop = FALSE]
      graphics::lines(dat$mean_predicted, dat$observed_rate, col = cols[as.character(h)], lwd = 1.8)
      graphics::points(
        dat$mean_predicted,
        dat$observed_rate,
        col = cols[as.character(h)],
        pch = pchs[[model_name]],
        cex = pmax(0.7, sqrt(dat$n) / 3)
      )
    }
    graphics::legend(
      "topleft",
      legend = paste0(horizons, " mo"),
      col = cols,
      pch = 16,
      lwd = 1.8,
      bty = "n",
      cex = 0.85
    )
  }
}

plot_decision_curve <- function(dca, path, device = c("pdf", "png")) {
  device <- match.arg(device)
  if (device == "pdf") {
    grDevices::pdf(path, width = 7.2, height = 4.8, useDingbats = FALSE)
  } else {
    grDevices::png(path, width = 7.2, height = 4.8, units = "in", res = 600)
  }
  on.exit(grDevices::dev.off(), add = TRUE)

  horizons <- sort(unique(dca$horizon_month))
  cols <- ggsci::pal_nejm("default")(max(3, length(horizons)))[seq_along(horizons)]   # ggsci NEJM
  names(cols) <- horizons
  yr <- range(c(dca$net_benefit, dca$treat_all_net_benefit, 0), finite = TRUE)

  graphics::par(mar = c(4.2, 4.5, 1.0, 0.7), mgp = c(2.6, 0.8, 0), las = 1)
  graphics::plot(
    NA,
    xlim = range(dca$threshold),
    ylim = yr,
    xlab = "Risk threshold",
    ylab = "Net benefit",
    bty = "l"
  )
  graphics::abline(h = 0, col = "gray55", lty = 2)
  for (h in horizons) {
    dat <- dca[dca$horizon_month == h, , drop = FALSE]
    dat <- dat[order(dat$threshold), , drop = FALSE]
    graphics::lines(dat$threshold, dat$net_benefit, col = cols[as.character(h)], lwd = 2.2)
    graphics::lines(dat$threshold, dat$treat_all_net_benefit, col = cols[as.character(h)], lwd = 1.1, lty = 3)
  }
  base <- dca[dca$horizon_month == horizons[1], , drop = FALSE]
  graphics::lines(base$threshold, base$treat_none_net_benefit, col = "gray30", lty = 2, lwd = 1.5)
  graphics::legend(
    "topright",
    legend = c(paste0(horizons, " mo model"), "Treat all", "Treat none"),
    col = c(cols, "gray45", "gray30"),
    lwd = c(rep(2.2, length(horizons)), 1.1, 1.5),
    lty = c(rep(1, length(horizons)), 3, 2),
    bty = "n",
    cex = 0.85
  )
}

main <- function() {
  patient <- read_processed("real_patient_level_analysis.csv")
  long <- read_processed("real_longitudinal_analysis.csv")
  patient <- patient[order(patient$patient_num), , drop = FALSE]
  long <- long[order(long$patient_num, long$visit_index), , drop = FALSE]

  landmarks_months <- c(6, 12, 18, 24)
  horizons_months <- c(6, 12, 24)
  max_draws <- as.numeric(Sys.getenv("GLW_DYNAMIC_MAX_DRAWS", unset = "Inf"))
  bootstrap_reps <- as.integer(Sys.getenv("GLW_DYNAMIC_BOOTSTRAP_REPS", unset = "200"))

  draws <- read_draw_subset(DRAW_DIR, n_pat = nrow(patient), max_draws = max_draws)
  landmark_data <- make_landmark_dataset(patient, long, landmarks_months, horizons_months)
  landmark_data <- predict_dynamic_probabilities(landmark_data, draws)

  coef_rows <- list()
  pred_parts <- list()
  for (H in horizons_months) {
    dat <- landmark_data[landmark_data$horizon_month == H, , drop = FALSE]
    co <- fit_recalibration(dat)
    dat$pred_recalibrated <- apply_recalibration(dat$pred_original, co[["a"]], co[["b"]])
    coef_rows[[as.character(H)]] <- data.frame(
      horizon_month = H,
      n = nrow(dat),
      events = sum(dat$outcome_dmr_by_horizon),
      intercept_a = co[["a"]],
      slope_b = co[["b"]],
      stringsAsFactors = FALSE
    )
    pred_parts[[as.character(H)]] <- dat
  }
  predictions <- do.call(rbind, pred_parts)

  write.csv(
    predictions,
    file.path(MODEL_DIR, "dynamic_prediction_landmark_predictions.csv"),
    row.names = FALSE
  )
  coef_tab <- do.call(rbind, coef_rows)
  write.csv(
    coef_tab,
    file.path(MODEL_DIR, "dynamic_prediction_recalibration_coefficients.csv"),
    row.names = FALSE
  )

  validation_rows <- list()
  idx <- 0L
  for (level in c("pooled_by_horizon", "by_landmark_and_horizon")) {
    if (level == "pooled_by_horizon") {
      split_key <- split(predictions, predictions$horizon_month)
    } else {
      split_key <- split(predictions, paste(predictions$landmark_month, predictions$horizon_month, sep = "_"))
    }
    for (key in names(split_key)) {
      dat <- split_key[[key]]
      for (model in c("Original posterior", "Recalibrated")) {
        pred_col <- ifelse(model == "Original posterior", "pred_original", "pred_recalibrated")
        met <- metric_set(dat$outcome_dmr_by_horizon, dat[[pred_col]])
        idx <- idx + 1L
        validation_rows[[idx]] <- cbind(
          data.frame(
            validation_level = level,
            landmark_month = ifelse(level == "pooled_by_horizon", NA_real_, unique(dat$landmark_month)),
            horizon_month = unique(dat$horizon_month),
            model = model,
            stringsAsFactors = FALSE
          ),
          met
        )
      }
    }
  }
  validation <- do.call(rbind, validation_rows)

  optimism_rows <- lapply(horizons_months, function(H) {
    dat <- predictions[predictions$horizon_month == H, , drop = FALSE]
    opt <- bootstrap_recalibration_optimism(dat, B = bootstrap_reps)
    data.frame(
      horizon_month = H,
      opt,
      stringsAsFactors = FALSE
    )
  })
  optimism <- do.call(rbind, optimism_rows)
  write.csv(optimism, file.path(MODEL_DIR, "dynamic_prediction_bootstrap_optimism.csv"), row.names = FALSE)

  validation$bootstrap_reps <- NA_integer_
  validation$brier_optimism <- NA_real_
  validation$optimism_corrected_brier <- validation$brier
  for (H in horizons_months) {
    opt <- optimism[optimism$horizon_month == H, , drop = FALSE]
    row_idx <- validation$validation_level == "pooled_by_horizon" &
      validation$horizon_month == H &
      validation$model == "Recalibrated"
    validation$bootstrap_reps[row_idx] <- opt$bootstrap_reps
    validation$brier_optimism[row_idx] <- opt$brier_optimism
    validation$optimism_corrected_brier[row_idx] <- validation$brier[row_idx] + opt$brier_optimism
  }
  validation$note <- ifelse(
    validation$model == "Recalibrated" & validation$validation_level == "pooled_by_horizon",
    "Bootstrap optimism correction applies only to the recalibration layer with fixed posterior predictions.",
    "Apparent development-cohort estimate; full Bayesian model was not refitted."
  )
  validation <- validation[order(validation$validation_level, validation$horizon_month, validation$landmark_month, validation$model), ]
  write.csv(validation, file.path(MODEL_DIR, "dynamic_prediction_validation.csv"), row.names = FALSE)

  bins <- do.call(rbind, lapply(horizons_months, function(H) {
    dat <- predictions[predictions$horizon_month == H, , drop = FALSE]
    rbind(
      calibration_bins(dat, "pred_original", "Original posterior"),
      calibration_bins(dat, "pred_recalibrated", "Recalibrated")
    )
  }))
  write.csv(bins, file.path(MODEL_DIR, "dynamic_prediction_calibration_bins.csv"), row.names = FALSE)

  dca <- do.call(rbind, lapply(horizons_months, function(H) {
    dat <- predictions[predictions$horizon_month == H, , drop = FALSE]
    decision_curve(dat, "pred_recalibrated", "Recalibrated")
  }))
  write.csv(dca, file.path(MODEL_DIR, "dynamic_prediction_decision_curve.csv"), row.names = FALSE)

  write_validation_table_tex(validation, file.path(MODEL_DIR, "table_09_dynamic_prediction_validation.tex"))
  plot_calibration(bins, file.path(MODEL_DIR, "figure_13_dynamic_prediction_calibration.pdf"), "pdf")
  plot_calibration(bins, file.path(MODEL_DIR, "figure_13_dynamic_prediction_calibration.png"), "png")
  plot_decision_curve(dca, file.path(MODEL_DIR, "figure_14_dynamic_prediction_decision_curve.pdf"), "pdf")
  plot_decision_curve(dca, file.path(MODEL_DIR, "figure_14_dynamic_prediction_decision_curve.png"), "png")

  results_md <- c(
    "# Dynamic DMR Prediction and Recalibration",
    "",
    paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "",
    "## Scope",
    "",
    "This analysis uses posterior draws from the renewed Bayesian joint longitudinal--interval model to generate landmark-level probabilities of documented DMR by 6, 12, and 24 months after each landmark. Predictions are evaluated among patients who were DMR-free at the landmark, had at least one MRD measurement at or before the landmark, and had observable outcome status at the prediction horizon.",
    "",
    "Important limitation: these are development-cohort dynamic predictions using the already fitted posterior. The patient-specific random effects in the saved posterior were estimated using the full longitudinal dataset, so the numeric results should be interpreted as an apparent/prototype analysis. A strict prospective dynamic prediction analysis requires landmark-specific posterior updating using only MRD history available at time s, or leave-one-patient-out refitting.",
    "",
    "## Landmark and Horizon Definitions",
    "",
    paste0("- Landmarks: ", paste(landmarks_months, collapse = ", "), " months."),
    paste0("- Horizons: ", paste(horizons_months, collapse = ", "), " months."),
    "- Eligible records: DMR-free at landmark, at least one prior MRD measurement, and documented DMR by horizon or follow-up beyond the horizon.",
    "",
    "## Pooled Validation Summary",
    ""
  )
  pooled <- validation[validation$validation_level == "pooled_by_horizon", , drop = FALSE]
  for (H in horizons_months) {
    sub <- pooled[pooled$horizon_month == H, , drop = FALSE]
    orig <- sub[sub$model == "Original posterior", , drop = FALSE]
    cal <- sub[sub$model == "Recalibrated", , drop = FALSE]
    results_md <- c(
      results_md,
      paste0(
        "- ", H, "-month horizon: n=", orig$n,
        ", events=", orig$events,
        ", original mean predicted=", fmt(orig$mean_predicted, 3),
        " vs observed=", fmt(orig$observed_rate, 3),
        ", original Brier=", fmt(orig$brier, 3),
        "; recalibrated mean predicted=", fmt(cal$mean_predicted, 3),
        ", recalibrated Brier=", fmt(cal$brier, 3),
        ", optimism-corrected recalibrated Brier=", fmt(cal$optimism_corrected_brier, 3),
        "."
      )
    )
  }
  results_md <- c(
    results_md,
    "",
    "## Output Files",
    "",
    "- `dynamic_prediction_landmark_predictions.csv`",
    "- `dynamic_prediction_recalibration_coefficients.csv`",
    "- `dynamic_prediction_validation.csv`",
    "- `dynamic_prediction_bootstrap_optimism.csv`",
    "- `dynamic_prediction_calibration_bins.csv`",
    "- `dynamic_prediction_decision_curve.csv`",
    "- `table_09_dynamic_prediction_validation.tex`",
    "- `figure_13_dynamic_prediction_calibration.pdf/.png`",
    "- `figure_14_dynamic_prediction_decision_curve.pdf/.png`",
    ""
  )
  writeLines(results_md, file.path(MODEL_DIR, "dynamic_prediction_recalibration_report.md"), useBytes = TRUE)

  cat("Dynamic prediction analysis complete.\n")
  cat("Landmark prediction rows:", nrow(predictions), "\n")
  cat("Validation table:", file.path(MODEL_DIR, "dynamic_prediction_validation.csv"), "\n")
}

main()
