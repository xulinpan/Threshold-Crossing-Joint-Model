options(stringsAsFactors = FALSE)

source(file.path("04_Code", "R", "01_clean_generate_model_data.R"))

LOCAL_R_LIB <- file.path(PROJECT_ROOT, "04_Code", "R", "library")
if (dir.exists(LOCAL_R_LIB)) {
  .libPaths(c(LOCAL_R_LIB, .libPaths()))
}
for (pkg in c("RColorBrewer","ggsci")) if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  stop("RColorBrewer is required for figure palettes. Install it with install.packages('RColorBrewer').")
}

TABLE_DIR <- file.path(PROJECT_ROOT, "06_Tables")
FIGURE_DIR <- file.path(PROJECT_ROOT, "05_Figures")
LATEX_DIR <- file.path(PROJECT_ROOT, "02_LaTeX")
PNG_RESOLUTION_DPI <- 600
dir.create(TABLE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIGURE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LATEX_DIR, recursive = TRUE, showWarnings = FALSE)

long <- read.csv(file.path(PROCESSED_DIR, "real_longitudinal_analysis.csv"), stringsAsFactors = FALSE)
intervals <- read.csv(file.path(PROCESSED_DIR, "real_interval_survival_analysis.csv"), stringsAsFactors = FALSE)
patient <- read.csv(file.path(PROCESSED_DIR, "real_patient_level_analysis.csv"), stringsAsFactors = FALSE)
audit <- read.csv(file.path(PROCESSED_DIR, "data_cleaning_audit.csv"), stringsAsFactors = FALSE)
missingness <- read.csv(file.path(PROCESSED_DIR, "longitudinal_missingness.csv"), stringsAsFactors = FALSE)

brewer_dark <- ggsci::pal_lancet("lanonc")(8)      # ggsci Lancet (categorical)
brewer_set  <- ggsci::pal_npg("nrc")(8)            # ggsci NPG (categorical)
brewer_blue <- RColorBrewer::brewer.pal(9, "Blues")
brewer_orange <- RColorBrewer::brewer.pal(9, "YlOrRd")
brewer_grey <- RColorBrewer::brewer.pal(9, "Greys")
plot_col <- list(
  grid = brewer_grey[3],
  axis = brewer_grey[8],
  trajectory = grDevices::adjustcolor(brewer_grey[6], alpha.f = 0.24),
  trajectory_legend = grDevices::adjustcolor(brewer_grey[7], alpha.f = 0.62),
  trend = brewer_dark[1],
  dmr = brewer_dark[2],
  cmr = brewer_dark[3],
  km = brewer_dark[4],
  box = grDevices::adjustcolor(brewer_set[2], alpha.f = 0.56),
  box_border = brewer_grey[7],
  point = grDevices::adjustcolor(brewer_set[7], alpha.f = 0.52),
  median = brewer_dark[3],
  visit_hist = brewer_blue[6],
  gap_hist = brewer_orange[5],
  missing = brewer_orange[6]
)

fmt_num <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", formatC(x, format = "f", digits = digits))
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%s (%s); n=%d", fmt_num(mean(x), digits), fmt_num(stats::sd(x), digits), length(x))
}

fmt_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
  sprintf("%s [%s, %s]; n=%d", fmt_num(q[2], digits), fmt_num(q[1], digits), fmt_num(q[3], digits), length(x))
}

fmt_n_pct <- function(x, value = 1) {
  denom <- sum(!is.na(x))
  if (denom == 0) return("NA")
  num <- sum(x == value, na.rm = TRUE)
  sprintf("%d (%.1f%%); n=%d", num, 100 * num / denom, denom)
}

latex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("&", "\\\\&", x)
  x <- gsub("%", "\\\\%", x)
  x <- gsub("\\$", "\\\\$", x)
  x <- gsub("#", "\\\\#", x)
  x <- gsub("_", "\\\\_", x)
  x <- gsub("\\{", "\\\\{", x)
  x <- gsub("\\}", "\\\\}", x)
  x <- gsub("<", "$<$", x, fixed = TRUE)
  x <- gsub(">", "$>$", x, fixed = TRUE)
  x
}

write_latex_table <- function(df, file, caption, label, align = NULL) {
  if (is.null(align)) {
    align <- paste0("l", paste(rep("r", ncol(df) - 1), collapse = ""))
  }
  header <- paste(latex_escape(names(df)), collapse = " & ")
  body <- apply(df, 1, function(row) paste(latex_escape(row), collapse = " & "))
  lines <- c(
    "\\begin{table}[!htbp]",
    "\\centering",
    sprintf("\\caption{%s}", latex_escape(caption)),
    sprintf("\\label{%s}", label),
    "\\resizebox{\\linewidth}{!}{%",
    sprintf("\\begin{tabular}{%s}", align),
    "\\toprule",
    paste0(header, " \\\\"),
    "\\midrule",
    paste0(body, " \\\\"),
    "\\bottomrule",
    "\\end{tabular}",
    "}%",
    "\\end{table}"
  )
  writeLines(lines, file, useBytes = TRUE)
}

save_plot <- function(stem, plot_fun, width = 7.0, height = 5.0) {
  pdf_file <- file.path(FIGURE_DIR, paste0(stem, ".pdf"))
  png_file <- file.path(FIGURE_DIR, paste0(stem, ".png"))
  grDevices::pdf(pdf_file, width = width, height = height, useDingbats = FALSE)
  plot_fun()
  grDevices::dev.off()
  grDevices::png(png_file, width = width, height = height, units = "in", res = PNG_RESOLUTION_DPI)
  plot_fun()
  grDevices::dev.off()
  invisible(c(pdf = pdf_file, png = png_file))
}

group_data <- list(
  All = patient,
  DMR = patient[patient$dmr_event == 1, , drop = FALSE],
  Censored = patient[patient$dmr_event == 0, , drop = FALSE]
)

baseline_rows <- list(
  "Patients, n" = function(d) as.character(nrow(d)),
  "Age, mean (SD)" = function(d) fmt_mean_sd(d$age),
  "Male sex, n (%)" = function(d) fmt_n_pct(d$sex_male),
  "Imatinib duration, years, median [IQR]" = function(d) fmt_median_iqr(d$duration_years),
  "Baseline PH+ %, median [IQR]" = function(d) fmt_median_iqr(d$ph_baseline_pct),
  "3-month PH+ %, median [IQR]" = function(d) fmt_median_iqr(d$ph_3m_pct),
  "6-month PH+ %, median [IQR]" = function(d) fmt_median_iqr(d$ph_6m_pct),
  "12-month PH+ %, median [IQR]" = function(d) fmt_median_iqr(d$ph_12m_pct),
  "Follow-up, months, median [IQR]" = function(d) fmt_median_iqr(d$followup_months),
  "Visits per patient, median [IQR]" = function(d) fmt_median_iqr(d$n_visits, digits = 0),
  "DMR observed, n (%)" = function(d) fmt_n_pct(d$dmr_event)
)

baseline_table <- data.frame(
  Characteristic = names(baseline_rows),
  All = vapply(baseline_rows, function(fn) fn(group_data$All), character(1)),
  DMR = vapply(baseline_rows, function(fn) fn(group_data$DMR), character(1)),
  Censored = vapply(baseline_rows, function(fn) fn(group_data$Censored), character(1)),
  check.names = FALSE
)
utils::write.csv(baseline_table, file.path(TABLE_DIR, "table_01_baseline_characteristics.csv"), row.names = FALSE)
write_latex_table(
  baseline_table,
  file.path(TABLE_DIR, "table_01_baseline_characteristics.tex"),
  "Baseline characteristics and follow-up summaries.",
  "tab:baseline",
  align = "llll"
)

windows <- cut(
  long$t_months,
  breaks = c(-0.001, 3, 6, 12, 24, 60, Inf),
  labels = c("0-3", ">3-6", ">6-12", ">12-24", ">24-60", ">60")
)
long$time_window <- as.character(windows)
time_levels <- levels(windows)
longitudinal_table <- do.call(rbind, lapply(time_levels, function(win) {
  d <- long[long$time_window == win, , drop = FALSE]
  data.frame(
    Window_months = win,
    Observations = nrow(d),
    Patients = length(unique(d$patient_id)),
    BM_samples = fmt_n_pct(d$sample_bm),
    LOG_MRD = fmt_median_iqr(d$log_mrd),
    DMR = fmt_n_pct(d$dmr),
    CMR = fmt_n_pct(d$cmr),
    stringsAsFactors = FALSE
  )
}))
utils::write.csv(longitudinal_table, file.path(TABLE_DIR, "table_02_longitudinal_by_time.csv"), row.names = FALSE)
write_latex_table(
  longitudinal_table,
  file.path(TABLE_DIR, "table_02_longitudinal_by_time.tex"),
  "Longitudinal molecular monitoring summaries by follow-up window.",
  "tab:longitudinal-time",
  align = "lrrrrrr"
)

audit_table <- audit
names(audit_table) <- c("Cleaning step", "Count")
utils::write.csv(audit_table, file.path(TABLE_DIR, "table_03_data_cleaning_audit.csv"), row.names = FALSE)
write_latex_table(
  audit_table,
  file.path(TABLE_DIR, "table_03_data_cleaning_audit.tex"),
  "Data cleaning audit and generated analysis datasets.",
  "tab:cleaning-audit",
  align = "lr"
)

missing_table <- missingness[order(-missingness$missing_or_invalid), , drop = FALSE]
names(missing_table) <- c("Variable", "Missing or invalid rows")
utils::write.csv(missing_table, file.path(TABLE_DIR, "table_04_raw_missingness.csv"), row.names = FALSE)
write_latex_table(
  missing_table,
  file.path(TABLE_DIR, "table_04_raw_missingness.tex"),
  "Missing or invalid fields in the raw longitudinal file.",
  "tab:missingness",
  align = "lr"
)

eda_lines <- c(
  "# GLW Exploratory Data Analysis Summary",
  "",
  paste("- Model-ready longitudinal observations:", nrow(long)),
  paste("- Patients with longitudinal data:", length(unique(long$patient_id))),
  paste("- At-risk survival intervals:", nrow(intervals)),
  paste("- DMR event patients:", sum(patient$dmr_event, na.rm = TRUE)),
  paste("- Censored patients:", sum(patient$dmr_event == 0, na.rm = TRUE)),
  paste("- Median follow-up months:", fmt_num(stats::median(patient$followup_months, na.rm = TRUE))),
  paste("- Median visits per patient:", fmt_num(stats::median(patient$n_visits, na.rm = TRUE), digits = 0)),
  paste("- Patients with complete core covariates:", sum(patient$has_complete_covariates == 1, na.rm = TRUE)),
  "",
  "Core modeling files contain no direct patient names. The reversible patient-name map is stored only in `03_Data/Processed/private/patient_key.csv`."
)
writeLines(eda_lines, file.path(TABLE_DIR, "eda_summary.md"), useBytes = TRUE)

plot_trajectories <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(4.3, 4.6, 1.5, 1.0), family = "sans", las = 1, bty = "l", fg = plot_col$axis)
  plot(
    range(long$t_months, na.rm = TRUE), range(long$log_mrd, na.rm = TRUE),
    type = "n", xlab = "Months since imatinib start", ylab = "LOG-MRD",
    xaxs = "i"
  )
  grid(col = plot_col$grid, lty = 1)
  for (pid in unique(long$patient_id)) {
    d <- long[long$patient_id == pid, , drop = FALSE]
    d <- d[order(d$t_months), , drop = FALSE]
    lines(d$t_months, d$log_mrd, col = plot_col$trajectory, lwd = 0.85)
  }
  sm <- stats::lowess(long$t_months, long$log_mrd, f = 0.30)
  lines(sm$x, sm$y, col = plot_col$trend, lwd = 3.2)
  abline(h = -4.5, col = plot_col$dmr, lty = 2, lwd = 1.6)
  abline(h = -5.0, col = plot_col$cmr, lty = 3, lwd = 1.6)
  legend(
    "topright",
    legend = c("Patient trajectory", "Smoothed trend", "DMR threshold", "CMR threshold"),
    col = c(plot_col$trajectory_legend, plot_col$trend, plot_col$dmr, plot_col$cmr),
    lty = c(1, 1, 2, 3), lwd = c(1, 3, 1.5, 1.5), bty = "n", cex = 0.9
  )
}

plot_km <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(4.3, 4.6, 1.5, 1.0), family = "sans", las = 1, bty = "l", fg = plot_col$axis)
  fit <- survival::survfit(survival::Surv(time_to_dmr_or_censor, dmr_event) ~ 1, data = patient)
  plot(
    fit, conf.int = FALSE, mark.time = TRUE, col = plot_col$km, lwd = 2.6,
    xlab = "Months since imatinib start", ylab = "Probability without observed DMR",
    xaxs = "i", yaxs = "i"
  )
  grid(col = plot_col$grid, lty = 1)
  legend(
    "topright",
    legend = sprintf("Events: %d / %d", sum(patient$dmr_event), nrow(patient)),
    bty = "n"
  )
}

plot_cytogenetic <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mar = c(4.3, 4.6, 1.5, 1.0), family = "sans", las = 1, bty = "l", fg = plot_col$axis)
  cols <- c("ph_baseline_pct", "ph_3m_pct", "ph_6m_pct", "ph_12m_pct", "ph_18m_pct", "ph_2y_pct")
  labels <- c("0", "3", "6", "12", "18", "24")
  values <- lapply(cols, function(nm) patient[[nm]][!is.na(patient[[nm]])])
  boxplot(
    values, names = labels, outline = FALSE, ylim = c(0, 105),
    xlab = "Months", ylab = "Philadelphia chromosome positive cells (%)",
    col = plot_col$box,
    border = plot_col$box_border
  )
  grid(nx = NA, ny = NULL, col = plot_col$grid, lty = 1)
  stripchart(values, vertical = TRUE, method = "jitter", pch = 16, cex = 0.55,
             col = plot_col$point, add = TRUE)
  med <- vapply(values, stats::median, numeric(1), na.rm = TRUE)
  lines(seq_along(med), med, type = "b", pch = 18, lwd = 2.3, col = plot_col$median)
}

plot_visit_distributions <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  par(mfrow = c(1, 2), mar = c(4.3, 4.5, 2.0, 1.0), family = "sans", las = 1, bty = "l", fg = plot_col$axis)
  hist(long$t_months, breaks = 24, col = plot_col$visit_hist, border = "white",
       xlab = "Months since start", main = "Visit timing")
  grid(col = plot_col$grid, lty = 1)
  hist(long$gap_months[long$gap_months > 0], breaks = 24, col = plot_col$gap_hist, border = "white",
       xlab = "Months since previous visit", main = "Visit gaps")
  grid(col = plot_col$grid, lty = 1)
}

plot_missingness <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  d <- missingness[order(missingness$missing_or_invalid), , drop = FALSE]
  par(mar = c(4.0, 8.5, 1.5, 1.0), family = "sans", las = 1, bty = "l", fg = plot_col$axis)
  barplot(
    d$missing_or_invalid, names.arg = d$variable, horiz = TRUE,
    col = plot_col$missing, border = NA, xlab = "Raw rows with missing or invalid values"
  )
  grid(nx = NULL, ny = NA, col = plot_col$grid, lty = 1)
}

save_plot("figure_01_log_mrd_trajectories", plot_trajectories, width = 7.2, height = 5.2)
save_plot("figure_02_time_to_dmr_km", plot_km, width = 6.5, height = 5.0)
save_plot("figure_03_cytogenetic_response", plot_cytogenetic, width = 6.8, height = 5.0)
save_plot("figure_04_visit_timing_and_gaps", plot_visit_distributions, width = 8.4, height = 4.3)
save_plot("figure_05_raw_missingness", plot_missingness, width = 7.0, height = 4.8)

report_pct <- function(num, den) {
  if (is.na(den) || den == 0) return("NA")
  sprintf("%.1f\\%%", 100 * num / den)
}

report_median_iqr <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  q <- stats::quantile(x, probs = c(0.25, 0.5, 0.75), names = FALSE)
  sprintf("%s [%s, %s]", fmt_num(q[2], digits), fmt_num(q[1], digits), fmt_num(q[3], digits))
}

report_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("%s (SD %s; n=%d)", fmt_num(mean(x), digits), fmt_num(stats::sd(x), digits), length(x))
}

report_n_pct <- function(x, value = 1) {
  denom <- sum(!is.na(x))
  if (denom == 0) return("NA")
  num <- sum(x == value, na.rm = TRUE)
  sprintf("%d (%s; n=%d)", num, report_pct(num, denom), denom)
}

get_audit_value <- function(item) {
  value <- audit$value[audit$item == item]
  if (length(value) == 0) return(NA_integer_)
  value[[1]]
}

n_patients <- length(unique(long$patient_id))
n_events <- sum(patient$dmr_event, na.rm = TRUE)
n_censored <- sum(patient$dmr_event == 0, na.rm = TRUE)
n_complete_covariates <- sum(patient$has_complete_covariates == 1, na.rm = TRUE)
n_dmr_visits <- sum(long$dmr == 1, na.rm = TRUE)
n_cmr_visits <- sum(long$cmr == 1, na.rm = TRUE)
first_visits <- long[long$visit_index == 1, , drop = FALSE]
last_visit_rows <- do.call(rbind, lapply(split(long, long$patient_id), function(d) d[which.max(d$visit_index), , drop = FALSE]))
event_times <- patient$time_to_dmr_or_censor[patient$dmr_event == 1]
ph_cols <- c("ph_baseline_pct", "ph_3m_pct", "ph_6m_pct", "ph_12m_pct", "ph_18m_pct", "ph_2y_pct")
ph_available <- vapply(patient[ph_cols], function(x) sum(!is.na(x)), integer(1))
top_missing <- missingness[order(-missingness$missing_or_invalid), , drop = FALSE]
top_missing <- top_missing[seq_len(min(3, nrow(top_missing))), , drop = FALSE]
top_missing_text <- paste(sprintf("%s (%d rows)", latex_escape(top_missing$variable), top_missing$missing_or_invalid), collapse = ", ")
sex_summary_text <- report_n_pct(patient$sex_male)

report_file <- file.path(LATEX_DIR, "glw_data_analysis_report.tex")
report_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{booktabs}",
  "\\usepackage{graphicx}",
  "\\usepackage{float}",
  "\\usepackage{caption}",
  "\\usepackage{hyperref}",
  "\\captionsetup{font=small,labelfont=bf}",
  "\\title{Exploratory Data Analysis Report for the GLW CML Dataset}",
  "\\author{}",
  "\\date{\\today}",
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\begin{abstract}",
  sprintf("This report summarizes data cleaning and exploratory analysis for a chronic myeloid leukemia molecular monitoring dataset. The regenerated public analysis files contain %d longitudinal molecular observations from %d coded patients, %d at-risk interval-survival records, and %d patient-level records. A total of %d patients (%s) had an observed deep molecular response (DMR), while %d patients were censored without observed DMR. The report emphasizes privacy-preserving preprocessing, completeness of modeling variables, longitudinal LOG-MRD trajectories, time to first DMR, cytogenetic response patterns, and data limitations relevant to subsequent joint longitudinal--interval survival modeling.",
          nrow(long), n_patients, nrow(intervals), nrow(patient), n_events, report_pct(n_events, nrow(patient)), n_censored),
  "\\end{abstract}",
  "",
  "\\section{Analysis Objectives}",
  "The EDA had four practical objectives. First, the raw longitudinal and patient-level files were converted into privacy-preserving public analysis tables. Second, complete-case modeling fields were audited so that downstream statistical modeling receives coherent time, outcome, and molecular-measurement variables. Third, cohort-level summaries and high-resolution figures were generated to understand response dynamics and missingness. Fourth, the final LaTeX report, tables, figures, and verification checks were made reproducible from the R scripts in \\texttt{04\\_Code/R}.",
  "",
  "\\section{Data Sources and Cleaning}",
  sprintf("The raw visit-level file contained %d rows. After removing records with missing or invalid essential modeling fields, the longitudinal model file contains %d complete observations. The analysis preserves all direct patient identifiers only in the private key file under \\texttt{03\\_Data/Processed/private}; public analysis files use coded identifiers of the form \\texttt{P0001}.",
          get_audit_value("raw_longitudinal_rows"), nrow(long)),
  "",
  sprintf("Treatment time was represented as months since imatinib start. Baseline rows marked as pretreatment measurements were set to time zero; follow-up rows used the treatment-start date and laboratory date when available, with visit-label fallback rules for records that had explicit month labels. DMR was defined as \\texttt{log\\_mrd <= -4.5}; complete molecular response (CMR) was defined as \\texttt{log\\_mrd <= -5.0}."),
  "",
  "\\input{../06_Tables/table_03_data_cleaning_audit.tex}",
  "\\input{../06_Tables/table_04_raw_missingness.tex}",
  "",
  "\\section{Cohort Description}",
  sprintf("The patient-level analysis table contains %d coded patients. Baseline covariates are available for a subset: %d patients (%s) have complete core covariates for age, treatment duration, baseline Philadelphia chromosome percentage, and sex. Among patients with non-missing age, the mean age was %s years. Male patients accounted for %s of those with recorded sex.",
          nrow(patient), n_complete_covariates, report_pct(n_complete_covariates, nrow(patient)),
          report_mean_sd(patient$age), sex_summary_text),
  "",
  sprintf("Follow-up was heterogeneous, with median follow-up of %s months and a median of %s visits per patient. Patients with observed DMR had longer median observed follow-up (%s months) than censored patients (%s months), which is expected because longer observation creates more opportunity to document response.",
          report_median_iqr(patient$followup_months), report_median_iqr(patient$n_visits, digits = 0),
          report_median_iqr(patient$followup_months[patient$dmr_event == 1]),
          report_median_iqr(patient$followup_months[patient$dmr_event == 0])),
  "",
  "\\input{../06_Tables/table_01_baseline_characteristics.tex}",
  "",
  "\\section{Longitudinal Molecular Patterns}",
  sprintf("The longitudinal file contains %d complete molecular measurements. Baseline or first-observed LOG-MRD had median %s, while the last observed LOG-MRD per patient had median %s. Across all visits, %d observations (%s) met the DMR threshold and %d observations (%s) met the CMR threshold.",
          nrow(long), report_median_iqr(first_visits$log_mrd), report_median_iqr(last_visit_rows$log_mrd),
          n_dmr_visits, report_pct(n_dmr_visits, nrow(long)), n_cmr_visits, report_pct(n_cmr_visits, nrow(long))),
  "",
  "Figure~\\ref{fig:log-mrd-trajectories} shows strong early decline in the smoothed LOG-MRD trend, followed by a flatter low-level trajectory. Individual trajectories vary substantially, supporting the need for patient-specific random effects in downstream longitudinal modeling.",
  "",
  "\\input{../06_Tables/table_02_longitudinal_by_time.tex}",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.92\\linewidth]{../05_Figures/figure_01_log_mrd_trajectories.pdf}",
  "\\caption{Patient-level LOG-MRD trajectories with smoothed population trend and DMR/CMR thresholds.}",
  "\\label{fig:log-mrd-trajectories}",
  "\\end{figure}",
  "",
  "\\section{Time to Deep Molecular Response}",
  sprintf("DMR was observed for %d of %d patients (%s). Among event patients, the observed time to first DMR had median %s months. The interval-survival table contains %d at-risk intervals and encodes only the first interval ending in observed DMR for each patient.",
          n_events, nrow(patient), report_pct(n_events, nrow(patient)), report_median_iqr(event_times), nrow(intervals)),
  "",
  "The Kaplan--Meier curve in Figure~\\ref{fig:time-to-dmr} shows a rapid early decrease in the probability of remaining without observed DMR, followed by a long tail of delayed or unobserved responders. Because DMR is observed only at monitoring visits, interval-based survival modeling is more appropriate than treating event times as exact.",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.78\\linewidth]{../05_Figures/figure_02_time_to_dmr_km.pdf}",
  "\\caption{Kaplan--Meier curve for time until first observed DMR. Censoring marks indicate patients without observed DMR by last follow-up.}",
  "\\label{fig:time-to-dmr}",
  "\\end{figure}",
  "",
  "\\section{Cytogenetic Response}",
  sprintf("Cytogenetic measurements were most complete at baseline (n=%d) and became sparser over later landmarks (3 months n=%d, 6 months n=%d, 12 months n=%d, 18 months n=%d, 24 months n=%d). Median Philadelphia chromosome positivity was %s at baseline and %s at 12 months among available records.",
          ph_available[["ph_baseline_pct"]], ph_available[["ph_3m_pct"]], ph_available[["ph_6m_pct"]],
          ph_available[["ph_12m_pct"]], ph_available[["ph_18m_pct"]], ph_available[["ph_2y_pct"]],
          report_median_iqr(patient$ph_baseline_pct), report_median_iqr(patient$ph_12m_pct)),
  "",
  "Figure~\\ref{fig:cytogenetic-response} suggests a marked reduction in Philadelphia chromosome positivity after treatment initiation, though later landmark summaries should be interpreted cautiously because the available sample size decreases and follow-up schedules are irregular.",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.82\\linewidth]{../05_Figures/figure_03_cytogenetic_response.pdf}",
  "\\caption{Distribution of Philadelphia chromosome positive cells over follow-up among patients with available cytogenetic data.}",
  "\\label{fig:cytogenetic-response}",
  "\\end{figure}",
  "",
  "\\section{Visit Schedule and Missingness}",
  sprintf("The monitoring schedule was irregular. Visit times had median %s months from treatment start, and positive inter-visit gaps had median %s months. The most common missing or invalid raw longitudinal fields were %s.",
          report_median_iqr(long$t_months), report_median_iqr(long$gap_months[long$gap_months > 0]), top_missing_text),
  "",
  "Figure~\\ref{fig:visit-timing} shows that many measurements occurred early after treatment initiation, with progressively fewer late measurements. Figure~\\ref{fig:raw-missingness} shows that essential molecular fields were mostly complete, but missing LOG-MRD, ratio, laboratory date, specimen type, and ABL copy values drove most exclusions.",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.95\\linewidth]{../05_Figures/figure_04_visit_timing_and_gaps.pdf}",
  "\\caption{Monitoring time distribution and inter-visit gap distribution in the cleaned longitudinal data.}",
  "\\label{fig:visit-timing}",
  "\\end{figure}",
  "",
  "\\begin{figure}[!htbp]",
  "\\centering",
  "\\includegraphics[width=0.78\\linewidth]{../05_Figures/figure_05_raw_missingness.pdf}",
  "\\caption{Raw longitudinal missingness in essential modeling fields.}",
  "\\label{fig:raw-missingness}",
  "\\end{figure}",
  "",
  "\\section{Modeling Implications}",
  "\\begin{itemize}",
  "\\item The cleaned longitudinal data are suitable for modeling LOG-MRD trajectories with patient-specific intercepts and slopes.",
  "\\item The response endpoint should be handled as interval observed because DMR is detected at irregular monitoring visits rather than at a continuous exact event time.",
  "\\item Baseline covariate adjustment should be treated as secondary or sensitivity analysis because complete core covariates are available for only a subset of patients.",
  "\\item The public modeling datasets are privacy preserving; direct patient names are excluded from public outputs and retained only in the private linkage file.",
  "\\end{itemize}",
  "",
  "\\section{Generated Reproducible Outputs}",
  "The EDA pipeline generated the public model-ready CSV files under \\texttt{03\\_Data/Processed}, high-resolution 600 dpi PNG figures and vector PDF figures under \\texttt{05\\_Figures}, table CSV and LaTeX fragments under \\texttt{06\\_Tables}, and this compiled report under \\texttt{02\\_LaTeX}. The driver script \\texttt{04\\_Code/R/run\\_all.R} rebuilds the complete workflow and runs verification checks for completeness, privacy, nonempty artifacts, and interval-event consistency.",
  "",
  "\\end{document}"
)
writeLines(report_lines, report_file, useBytes = TRUE)

cat("Wrote EDA tables to:", TABLE_DIR, "\n")
cat("Wrote figures to:", FIGURE_DIR, "\n")
cat("Wrote LaTeX report to:", report_file, "\n")
