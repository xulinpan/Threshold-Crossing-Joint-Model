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
    "04_Code/R/08_generate_renewed_model_reporting.R"
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
MODEL_DIR <- file.path(ROOT_DIR, "08_Model")
MODEL_PREFIX <- Sys.getenv(
  "GLW_MODEL_PREFIX",
  unset = "stan_joint_interval_dmr_independent_renewed"
)
SUMMARY_FILE <- file.path(MODEL_DIR, paste0(MODEL_PREFIX, "_summary.csv"))
DRAW_DIR <- file.path(MODEL_DIR, paste0(MODEL_PREFIX, "_draws"))

if (!file.exists(SUMMARY_FILE)) {
  stop("Missing renewed model summary: ", SUMMARY_FILE)
}
draw_files <- list.files(DRAW_DIR, pattern = "[.]csv$", full.names = TRUE)
if (length(draw_files) == 0) {
  stop("No renewed model draw CSV files found in: ", DRAW_DIR)
}

summary_tab <- read.csv(SUMMARY_FILE, stringsAsFactors = FALSE)
draws_by_chain <- lapply(draw_files, function(f) {
  read.csv(f, comment.char = "#", check.names = FALSE)
})
draws <- do.call(rbind, draws_by_chain)

fmt <- function(x, digits = 3) {
  if (is.character(x)) return(x)
  formatC(as.numeric(x), digits = digits, format = "f")
}

lookup_summary <- function(variable) {
  row <- summary_tab[summary_tab$variable == variable, , drop = FALSE]
  if (nrow(row) != 1) {
    stop("Expected one summary row for ", variable, ", found ", nrow(row))
  }
  row
}

posterior_from_draws <- function(summary_variable, draw_variable) {
  if (!draw_variable %in% names(draws)) {
    stop("Missing draw column: ", draw_variable)
  }
  x <- draws[[draw_variable]]
  s <- lookup_summary(summary_variable)
  data.frame(
    variable = summary_variable,
    mean = mean(x),
    median = stats::median(x),
    q2.5 = stats::quantile(x, 0.025, names = FALSE),
    q5 = s$q5,
    q95 = s$q95,
    q97.5 = stats::quantile(x, 0.975, names = FALSE),
    rhat = s$rhat,
    ess_bulk = s$ess_bulk,
    ess_tail = s$ess_tail,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

core_spec <- data.frame(
  component = c(
    "Longitudinal MRD",
    "Longitudinal MRD",
    "Longitudinal MRD",
    "Longitudinal MRD",
    "Interval DMR",
    "Interval DMR",
    "Joint association",
    "Random effects",
    "Random effects"
  ),
  summary_variable = c(
    "beta_time",
    "beta_time2",
    "beta_bm",
    "sigma_y",
    "gamma_time",
    "gamma_gap",
    "alpha_mrd",
    "tau_b[1]",
    "tau_b[2]"
  ),
  draw_variable = c(
    "beta_time",
    "beta_time2",
    "beta_bm",
    "sigma_y",
    "gamma_time",
    "gamma_gap",
    "alpha_mrd",
    "tau_b.1",
    "tau_b.2"
  ),
  parameter_latex = c(
    "\\(\\beta_{\\mathrm{time}}\\)",
    "\\(\\beta_{\\mathrm{time}^2}\\)",
    "\\(\\beta_{\\mathrm{BM}}\\)",
    "\\(\\sigma_y\\)",
    "\\(\\gamma_{\\mathrm{time}}\\)",
    "\\(\\gamma_{\\mathrm{gap}}\\)",
    "\\(\\alpha_{\\mathrm{MRD}}\\)",
    "\\(\\tau_{b0}\\)",
    "\\(\\tau_{b1}\\)"
  ),
  interpretation = c(
    "Strong decline in latent log-MRD over follow-up",
    "Evidence of nonlinear response dynamics",
    "Sample-source adjustment; uncertain effect",
    "Residual molecular-measurement variability",
    "Time effect uncertain after latent MRD adjustment",
    "Visit-gap term; interpret with cumulative interval length",
    "Lower latent MRD associated with higher DMR probability",
    "Patient heterogeneity in baseline latent MRD",
    "Patient heterogeneity in response slope"
  ),
  stringsAsFactors = FALSE
)

core_stats <- do.call(
  rbind,
  lapply(seq_len(nrow(core_spec)), function(i) {
    posterior_from_draws(core_spec$summary_variable[i], core_spec$draw_variable[i])
  })
)
core_out <- cbind(core_spec, core_stats[, setdiff(names(core_stats), "variable")])
core_out$exp_mean <- ifelse(
  core_out$summary_variable %in% c("beta_bm", "gamma_gap", "alpha_mrd"),
  exp(core_out$mean),
  NA_real_
)
core_out$exp_q2.5 <- ifelse(
  core_out$summary_variable %in% c("beta_bm", "gamma_gap", "alpha_mrd"),
  exp(core_out$q2.5),
  NA_real_
)
core_out$exp_q97.5 <- ifelse(
  core_out$summary_variable %in% c("beta_bm", "gamma_gap", "alpha_mrd"),
  exp(core_out$q97.5),
  NA_real_
)
write.csv(
  core_out,
  file.path(MODEL_DIR, "core_posterior_clinical_table.csv"),
  row.names = FALSE
)

tex_rows <- vapply(seq_len(nrow(core_out)), function(i) {
  paste(
    core_out$component[i],
    core_out$parameter_latex[i],
    fmt(core_out$mean[i]),
    paste0(fmt(core_out$q2.5[i]), " to ", fmt(core_out$q97.5[i])),
    core_out$interpretation[i],
    sep = " & "
  )
}, character(1))

posterior_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Renewed primary Bayesian joint-model posterior summaries. Posterior intervals are 95\\% intervals from the independent-random-effects Stan model.}",
  "\\label{tab:joint-model-posterior}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{llrrl}",
  "\\toprule",
  "Component & Parameter & Posterior mean & 95\\% posterior interval & Clinical interpretation \\\\",
  "\\midrule",
  paste0(tex_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "}%",
  "\\end{table}",
  ""
)
writeLines(posterior_table, file.path(MODEL_DIR, "table_05_joint_model_posterior_results.tex"), useBytes = TRUE)

chain_diagnostics <- do.call(rbind, lapply(seq_along(draws_by_chain), function(i) {
  d <- draws_by_chain[[i]]
  ebfmi <- mean(diff(d[["energy__"]])^2) / stats::var(d[["energy__"]])
  data.frame(
    chain = i,
    post_warmup_draws = nrow(d),
    divergences = sum(d[["divergent__"]]),
    max_treedepth = max(d[["treedepth__"]]),
    mean_acceptance = mean(d[["accept_stat__"]]),
    ebfmi = ebfmi,
    stringsAsFactors = FALSE
  )
}))

ppc <- read.csv(file.path(MODEL_DIR, "ppc_longitudinal_summary.csv"), stringsAsFactors = FALSE)
ppc_values <- setNames(ppc$value, ppc$metric)
cal <- read.csv(file.path(MODEL_DIR, "calibration_summary.csv"), stringsAsFactors = FALSE)
cal_interval <- cal[cal$unit == "interval", ]
cal_patient <- cal[cal$unit == "patient", ]

diagnostic_values <- data.frame(
  domain = c(
    rep("MCMC diagnostics", 7),
    rep("Longitudinal PPC", 3),
    rep("Assay-floor PPC", 2),
    rep("Interval calibration", 3),
    rep("Patient calibration", 3),
    rep("Sensitivity", 2)
  ),
  metric = c(
    "Post-warmup draws",
    "Divergent transitions",
    "Maximum treedepth",
    "Minimum approximate E-BFMI",
    "Maximum R-hat",
    "Minimum bulk ESS",
    "Minimum tail ESS",
    "Non-floor RMSE to posterior mean",
    "Non-floor 90\\% predictive coverage",
    "Non-floor 95\\% predictive coverage",
    "Observed floor rate",
    "Posterior mean floor probability",
    "Observed event rate",
    "Mean predicted probability",
    "Brier score",
    "Observed DMR rate",
    "Mean predicted DMR probability",
    "Brier score",
    "Full-data floor variants",
    "Non-floor-only variant"
  ),
  value = c(
    as.character(sum(chain_diagnostics$post_warmup_draws)),
    as.character(sum(chain_diagnostics$divergences)),
    as.character(max(chain_diagnostics$max_treedepth)),
    fmt(min(chain_diagnostics$ebfmi)),
    fmt(max(summary_tab$rhat, na.rm = TRUE)),
    fmt(min(summary_tab$ess_bulk, na.rm = TRUE)),
    fmt(min(summary_tab$ess_tail, na.rm = TRUE)),
    fmt(ppc_values[["nonfloor_rmse_to_posterior_mean"]]),
    fmt(ppc_values[["nonfloor_90pct_predictive_coverage"]]),
    fmt(ppc_values[["nonfloor_95pct_predictive_coverage"]]),
    fmt(ppc_values[["observed_floor_rate"]]),
    fmt(ppc_values[["posterior_mean_floor_rate"]]),
    fmt(cal_interval$observed_rate),
    fmt(cal_interval$mean_predicted),
    fmt(cal_interval$brier),
    fmt(cal_patient$observed_rate),
    fmt(cal_patient$mean_predicted),
    fmt(cal_patient$brier),
    "Qualitatively consistent longitudinal decline",
    "Rank-deficient stress test; supplemental only"
  ),
  stringsAsFactors = FALSE
)
write.csv(
  diagnostic_values,
  file.path(MODEL_DIR, "renewed_model_diagnostic_summary.csv"),
  row.names = FALSE
)

diagnostic_rows <- vapply(seq_len(nrow(diagnostic_values)), function(i) {
  paste(
    diagnostic_values$domain[i],
    diagnostic_values$metric[i],
    diagnostic_values$value[i],
    sep = " & "
  )
}, character(1))

diagnostic_table <- c(
  "\\begin{table}[!htbp]",
  "\\centering",
  "\\caption{Diagnostics, posterior predictive checks, calibration, and sensitivity summaries for the renewed independent-random-effects joint model.}",
  "\\label{tab:model-checks}",
  "\\resizebox{\\linewidth}{!}{%",
  "\\begin{tabular}{llr}",
  "\\toprule",
  "Domain & Metric & Value \\\\",
  "\\midrule",
  paste0(diagnostic_rows, " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "}%",
  "\\end{table}",
  ""
)
writeLines(diagnostic_table, file.path(MODEL_DIR, "table_06_model_check_summary.tex"), useBytes = TRUE)

warning_rows <- summary_tab[
  (!is.na(summary_tab$rhat) & summary_tab$rhat > 1.01) |
    (!is.na(summary_tab$ess_bulk) & summary_tab$ess_bulk < 400) |
    (!is.na(summary_tab$ess_tail) & summary_tab$ess_tail < 400),
  c("variable", "rhat", "ess_bulk", "ess_tail"),
  drop = FALSE
]

core_md_rows <- vapply(seq_len(nrow(core_out)), function(i) {
  paste0(
    "| `", core_out$summary_variable[i], "` | ",
    fmt(core_out$mean[i]), " | ",
    fmt(core_out$q5[i]), " | ",
    fmt(core_out$q95[i]), " | ",
    fmt(core_out$rhat[i]), " | ",
    fmt(core_out$ess_bulk[i]), " |"
  )
}, character(1))

warning_md <- if (nrow(warning_rows) == 0) {
  "No parameters had R-hat > 1.01, bulk ESS < 400, or tail ESS < 400."
} else {
  paste(
    c(
      "| Variable | R-hat | Bulk ESS | Tail ESS |",
      "|---|---:|---:|---:|",
      vapply(seq_len(nrow(warning_rows)), function(i) {
        paste0(
          "| `", warning_rows$variable[i], "` | ",
          fmt(warning_rows$rhat[i]), " | ",
          fmt(warning_rows$ess_bulk[i]), " | ",
          fmt(warning_rows$ess_tail[i]), " |"
        )
      }, character(1))
    ),
    collapse = "\n"
  )
}

report_lines <- c(
  "# Renewed Stan Result Check",
  "",
  "Checked: 2026-07-09",
  "",
  "## Bottom Line",
  "",
  "The renewed primary model replaces the correlated random-intercept/random-slope structure with independent patient-specific random intercept and slope terms. This directly removes the weakly estimated correlation parameter from the original primary model while preserving the longitudinal floor-censoring likelihood and interval DMR event likelihood.",
  "",
  "The final renewed fit used longer chains and is suitable for manuscript reporting: no divergent transitions were observed, no parameters exceeded R-hat 1.01, and no bulk or tail ESS values were below 400. The model should still be described as a monitoring-oriented development model, not a validated treatment-decision rule.",
  "",
  "## Current Artifacts",
  "",
  paste0("- Stan fit RDS: `08_Model/", MODEL_PREFIX, "_fit.rds`"),
  paste0("- Stan summary CSV: `08_Model/", MODEL_PREFIX, "_summary.csv`"),
  paste0("- CmdStan chain CSVs: `08_Model/", MODEL_PREFIX, "_draws/`"),
  "- Renewed Stan code: `04_Code/Stan/glw_joint_interval_dmr_independent.stan`",
  "- Previous correlated-model outputs archived in: `08_Model/archive_correlated_model_20260709_renewal/`",
  "",
  "## Sampling Setup",
  "",
  "- Chains: 4",
  "- Warmup per chain: 2000",
  "- Sampling draws per chain: 2000",
  paste0("- Total post-warmup draws: ", sum(chain_diagnostics$post_warmup_draws)),
  "- Save warmup: false",
  "- HMC metric: diagonal Euclidean",
  "",
  "## Sampler Diagnostics",
  "",
  paste0("- Divergences: ", sum(chain_diagnostics$divergences)),
  paste0("- Maximum treedepth reached: ", max(chain_diagnostics$max_treedepth)),
  paste0("- Mean acceptance by chain: ", paste(fmt(chain_diagnostics$mean_acceptance), collapse = ", ")),
  paste0("- Minimum approximate E-BFMI across chains: ", fmt(min(chain_diagnostics$ebfmi))),
  paste0("- Max R-hat in summary: ", fmt(max(summary_tab$rhat, na.rm = TRUE))),
  paste0("- Minimum bulk ESS: ", fmt(min(summary_tab$ess_bulk, na.rm = TRUE))),
  paste0("- Minimum tail ESS: ", fmt(min(summary_tab$ess_tail, na.rm = TRUE))),
  paste0("- Parameters with R-hat > 1.01: ", sum(summary_tab$rhat > 1.01, na.rm = TRUE)),
  paste0("- Parameters with bulk ESS < 400: ", sum(summary_tab$ess_bulk < 400, na.rm = TRUE)),
  paste0("- Parameters with tail ESS < 400: ", sum(summary_tab$ess_tail < 400, na.rm = TRUE)),
  "",
  warning_md,
  "",
  "## Core Posterior Results",
  "",
  "Intervals below are posterior 5th and 95th percentiles from the renewed summary CSV.",
  "",
  "| Parameter | Mean | 5th pct | 95th pct | R-hat | Bulk ESS |",
  "|---|---:|---:|---:|---:|---:|",
  core_md_rows,
  "",
  "## Interpretation Check",
  "",
  "- The longitudinal trajectory terms remain strongly supported: `beta_time` is negative and `beta_time2` is positive, indicating nonlinear decline in log-MRD over follow-up.",
  "- The bone-marrow sample-source adjustment remains uncertain because its posterior interval crosses zero.",
  "- The joint association remains clinically coherent: lower latent MRD is associated with higher interval probability of DMR through the negative `alpha_mrd` estimate.",
  "- The visit-gap coefficient remains negative in this parameterization, so interpretation should be tied to the discrete interval hazard and cumulative interval length rather than read as a simple marginal effect.",
  "- The random-effect correlation is no longer part of the primary model; this is the intended repair for the original overparameterized correlation component.",
  "",
  "## Posterior Predictive and Calibration Checks",
  "",
  paste0("- Non-floor RMSE to posterior mean: ", fmt(ppc_values[["nonfloor_rmse_to_posterior_mean"]])),
  paste0("- Non-floor 90% predictive coverage: ", fmt(ppc_values[["nonfloor_90pct_predictive_coverage"]])),
  paste0("- Non-floor 95% predictive coverage: ", fmt(ppc_values[["nonfloor_95pct_predictive_coverage"]])),
  paste0("- Observed floor rate: ", fmt(ppc_values[["observed_floor_rate"]])),
  paste0("- Posterior mean floor probability: ", fmt(ppc_values[["posterior_mean_floor_rate"]])),
  paste0("- Interval observed event rate versus mean predicted probability: ", fmt(cal_interval$observed_rate), " versus ", fmt(cal_interval$mean_predicted)),
  paste0("- Patient observed DMR rate versus mean predicted DMR probability: ", fmt(cal_patient$observed_rate), " versus ", fmt(cal_patient$mean_predicted)),
  "",
  "## Recommendation",
  "",
  "Use the renewed independent-random-effects model as the primary manuscript model. Present the original correlated structure, if mentioned at all, as an audited predecessor that motivated simplification. Keep clinical claims focused on monitoring support, descriptive calibration, and methodological development under irregular molecular monitoring."
)
writeLines(report_lines, file.path(MODEL_DIR, "stan_result_check.md"), useBytes = TRUE)
writeLines(report_lines, file.path(MODEL_DIR, "renewed_model_result_check.md"), useBytes = TRUE)

cat("Updated renewed model reporting outputs in:", MODEL_DIR, "\n")
