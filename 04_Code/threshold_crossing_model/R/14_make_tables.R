## =============================================================================
## 14_make_tables.R  —  build the LaTeX tables for the simulation section
##
## Reads the aggregated CSVs and writes booktabs tables into outputs/tables/,
## which was empty until now (every table in the manuscript was hand-made, so
## nothing was traceable to a specific run).
##
## OUTPUT (outputs/tables/)
##   tab_sim_central.tex     bias (MCse) and coverage (MCse) by parameter x n
##   tab_sim_misspec.tex     coverage (MCse) by scenario x parameter, n = 150
##   tab_sbc.tex             SBC rank-uniformity test
##   tab_sim_diagnostics.tex convergence / diagnostic rates per cell
##   tables_provenance.txt   which run produced them, and when
##
## EVERY table carries the dgm_tag it came from. The summary CSV can contain
## several mechanisms at once (legacy fits are still on disk), and silently
## pooling them would be the worst kind of error, so exactly one tag is
## selected and it is recorded in the caption.
##
## RUN (from the R/ folder):
##   Rscript 14_make_tables.R
##   $env:DGM_TAG='v1_legacy_legacy_hazard'; Rscript 14_make_tables.R   # legacy
## =============================================================================

.tab_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
OUT  <- normalizePath(file.path(.tab_dir, "..", "outputs"), mustWork = FALSE)
SIMD <- file.path(OUT, "sim_redesign")
TABD <- file.path(OUT, "tables")
dir.create(TABD, recursive = TRUE, showWarnings = FALSE)

## ---- what to report ---------------------------------------------------------
## Ordered; only these appear in the tables. Everything else stays in the CSVs.
PLAB <- c(
  alpha   = "$\\alpha_{\\mathrm{MRD}}$",
  beta1   = "$\\beta_{\\mathrm{time}}$",
  beta2   = "$\\beta_{\\mathrm{time}^2}$",
  gamma2  = "$\\gamma_{\\mathrm{gap}}$",
  sigma_y = "$\\sigma_y$",
  tau0    = "$\\tau_{b0}$",
  tau1    = "$\\tau_{b1}$")
MISSPEC_PARAMS <- c("alpha", "beta1", "sigma_y", "tau1")

## ---- helpers ----------------------------------------------------------------
fmt <- function(x, d = 3, dash = "---") {
  if (length(x) == 0 || is.na(x)) return(dash)
  formatC(as.numeric(x), format = "f", digits = d)
}
## value with its Monte Carlo SE in parentheses
fmt_mc <- function(v, se, d = 3, dash = "---") {
  if (length(v) == 0 || is.na(v)) return(dash)
  if (is.na(se)) return(fmt(v, d))
  sprintf("%s (%s)", fmt(v, d), fmt(se, d))
}
tex_write <- function(lines, file) {
  con <- file(file, "w"); writeLines(lines, con); close(con)
  message("  wrote ", basename(file))
}
need <- function(f) {
  if (!file.exists(f)) { message("SKIP: missing ", basename(f)); return(FALSE) }
  TRUE
}

PROV <- character(0)
note <- function(...) PROV <<- c(PROV, paste0(...))

## ---- select exactly one mechanism ------------------------------------------
sfile <- file.path(SIMD, "sim_redesign_summary.csv")
have_sim <- need(sfile)
TAG <- NA_character_
if (have_sim) {
  S <- utils::read.csv(sfile, stringsAsFactors = FALSE)
  ## Assert the schema up front. Below, columns are reached with `$`, which does
  ## PARTIAL matching on data.frames -- if 07 ever renames a column, `$` would
  ## silently resolve to a different one instead of returning NULL. Failing here
  ## with a named list is far better than producing a plausible wrong table.
  SIM_REQ <- c("dgm_tag","scen_id","label","tag","n","param","nrep","bias",
               "bias_mcse","coverage","cov_mcse","nfit","rate_bad",
               "rate_marginal","max_rhat","min_ebfmi","tot_div")
  if (length(setdiff(SIM_REQ, names(S))))
    stop("sim_redesign_summary.csv is missing: ",
         paste(setdiff(SIM_REQ, names(S)), collapse = ", "),
         "\n  Re-run 07_aggregate_sim.R.")
  tags <- sort(unique(S$dgm_tag))
  TAG <- Sys.getenv("DGM_TAG", "")
  if (!nzchar(TAG)) {
    ## default to the production mechanism if present, else the only one there
    TAG <- if ("v2_cohort_posterior_hazard" %in% tags) "v2_cohort_posterior_hazard"
           else tags[length(tags)]
  }
  if (!TAG %in% tags)
    stop("DGM_TAG '", TAG, "' not in the summary. Available: ", paste(tags, collapse = ", "))
  if (length(tags) > 1)
    message("NOTE: summary holds ", length(tags), " mechanisms (",
            paste(tags, collapse = ", "), "); reporting '", TAG, "' only.")
  S <- S[S$dgm_tag == TAG, ]
  note("dgm_tag          : ", TAG)
  note("other mechanisms : ", paste(setdiff(tags, TAG), collapse = ", "))
  note("source           : ", sfile)
  note("summary mtime    : ", format(file.mtime(sfile)))
}

## ---- 1. central cells: bias and coverage by parameter x n -------------------
if (have_sim) {
  C <- S[S$tag == "central", ]
  ns <- sort(unique(C$n))
  ncol_spec <- paste(rep("c", 2 * length(ns)), collapse = "")
  hdr_n <- paste(sprintf("\\multicolumn{2}{c}{$n=%d$}", ns), collapse = " & ")
  hdr_2 <- paste(rep("Bias (MCse) & Cov. (MCse)", length(ns)), collapse = " & ")
  cmid  <- paste(sprintf("\\cmidrule(lr){%d-%d}", 2 + 2*(seq_along(ns)-1), 3 + 2*(seq_along(ns)-1)),
                 collapse = " ")
  body <- character(0)
  for (p in names(PLAB)) {
    cells <- character(0)
    for (nn in ns) {
      r <- C[C$param == p & C$n == nn, ]
      if (nrow(r) == 0) { cells <- c(cells, "---", "---"); next }
      cells <- c(cells,
                 fmt_mc(r$bias[1],     r$bias_mcse[1], 3),
                 fmt_mc(r$coverage[1], r$cov_mcse[1],  3))
    }
    body <- c(body, paste0(PLAB[[p]], " & ", paste(cells, collapse = " & "), " \\\\"))
  }
  nrep_txt <- paste(sprintf("$n=%d$: %d", ns,
    vapply(ns, function(nn) suppressWarnings(max(C$nrep[C$n == nn], na.rm = TRUE)), numeric(1))),
    collapse = ", ")
  tex_write(c(
    "% generated by 14_make_tables.R -- do not edit by hand",
    sprintf("%% dgm_tag: %s", TAG),
    "\\begin{table}[!htbp]", "\\centering",
    "\\caption{Operating characteristics of the Bayesian joint interval-hazard",
    " estimator under correct specification. Bias is the mean posterior mean minus",
    " the data-generating value; coverage is of the equal-tailed $95\\%$ credible",
    " interval. Monte Carlo standard errors in parentheses. Replicates per cell:",
    paste0(" ", nrep_txt, ". Fits failing convergence diagnostics ($\\hat{R}>1.05$,"),
    " E-BFMI $<0.2$, or any divergent transition) are excluded; see",
    " Table~\\ref{tab:sim-diagnostics}.}",
    "\\label{tab:sim-central}", "\\small",
    sprintf("\\begin{tabular}{l%s}", ncol_spec), "\\toprule",
    paste0("Parameter & ", hdr_n, " \\\\"), cmid,
    paste0(" & ", hdr_2, " \\\\"), "\\midrule",
    body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(TABD, "tab_sim_central.tex"))
}

## ---- 2. misspecification cells ---------------------------------------------
if (have_sim) {
  M <- S[S$tag == "misspec", ]
  if (nrow(M)) {
    scen <- unique(M[order(M$scen_id), c("scen_id", "label")])
    body <- character(0); any_dash <- FALSE
    for (i in seq_len(nrow(scen))) {
      cells <- character(0)
      for (p in MISSPEC_PARAMS) {
        r <- M[M$scen_id == scen$scen_id[i] & M$param == p, ]
        if (nrow(r) == 0) { cells <- c(cells, "---"); any_dash <- TRUE; next }
        cells <- c(cells, fmt_mc(r$coverage[1], r$cov_mcse[1], 3))
      }
      d <- M[M$scen_id == scen$scen_id[i], ]
      conv <- fmt(1 - suppressWarnings(max(d$rate_bad, na.rm = TRUE)), 2)
      body <- c(body, paste0(gsub("&", "\\\\&", scen$label[i]), " & ",
                             paste(cells, collapse = " & "), " & ", conv, " \\\\"))
    }
    dash_note <- if (any_dash) paste0(
      " A dash marks a parameter with no counterpart in the data-generating",
      " mechanism for that departure, for which bias and coverage are therefore",
      " undefined rather than poor: under a monotone trajectory there is no",
      " quadratic coefficient, and under an exponential trajectory the time slope",
      " and random-slope SD are not on a comparable scale. Reporting them would",
      " describe the absence of a parameter, not the performance of the estimator.")
      else ""
    tex_write(c(
      "% generated by 14_make_tables.R -- do not edit by hand",
      sprintf("%% dgm_tag: %s", TAG),
      "\\begin{table}[!htbp]", "\\centering",
      "\\caption{Robustness to misspecification: coverage of the $95\\%$ credible",
      " interval under one departure at a time ($n=150$), with Monte Carlo standard",
      paste0(" errors in parentheses.", dash_note,
             " The final column is the proportion of fits meeting all convergence",
             " diagnostics.}"),
      "\\label{tab:sim-misspec}", "\\small",
      sprintf("\\begin{tabular}{l%sc}", paste(rep("c", length(MISSPEC_PARAMS)), collapse="")),
      "\\toprule",
      paste0("Departure & ", paste(PLAB[MISSPEC_PARAMS], collapse = " & "),
             " & Converged \\\\"),
      "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
      file.path(TABD, "tab_sim_misspec.tex"))
  }
}

## ---- 3. diagnostics --------------------------------------------------------
if (have_sim) {
  D <- unique(S[, c("scen_id","label","n","nfit","rate_bad","rate_marginal",
                    "max_rhat","min_ebfmi","tot_div")])
  D <- D[order(D$scen_id), ]
  body <- apply(D, 1, function(r) sprintf("%s & %s & %s & %s & %s & %s & %s \\\\",
    gsub("&", "\\\\&", r[["label"]]), r[["nfit"]],
    fmt(as.numeric(r[["rate_bad"]]), 3), fmt(as.numeric(r[["rate_marginal"]]), 3),
    fmt(as.numeric(r[["max_rhat"]]), 3), fmt(as.numeric(r[["min_ebfmi"]]), 3),
    r[["tot_div"]]))
  tex_write(c(
    "% generated by 14_make_tables.R -- do not edit by hand",
    sprintf("%% dgm_tag: %s", TAG),
    "\\begin{table}[!htbp]", "\\centering",
    "\\caption{Sampler diagnostics per simulation cell. ``Excluded'' is the",
    " proportion of fits with $\\hat{R}>1.05$, E-BFMI $<0.2$ or any divergent",
    " transition, which are dropped from the performance measures; ``marginal'' is",
    " the proportion with $1.01<\\hat{R}\\le1.05$, which are retained. The two are",
    " reported separately because a single combined threshold at $\\hat{R}>1.01$",
    " conflates benign autocorrelation at $4\\times500$ draws with genuine",
    " non-convergence.}",
    "\\label{tab:sim-diagnostics}", "\\small",
    "\\begin{tabular}{lcccccc}", "\\toprule",
    "Cell & Fits & Excluded & Marginal & $\\max\\hat{R}$ & $\\min$ E-BFMI & Divergences \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(TABD, "tab_sim_diagnostics.tex"))
}

## ---- 4. SBC ----------------------------------------------------------------
## Pick by SCHEMA, not by existence. sbc_summary.csv may be a stale file from
## before the binning fix, carrying only (param, nsbc, chisq) -- and computed with
## equal expected counts, so its chi-square values are inflated and must not be
## reported. Take the first candidate that has every column the table needs.
SBC_REQ <- c("param", "nsbc", "chisq", "df", "p")
sbc_candidates <- file.path(OUT, "sbc",
                            c("sbc_summary.csv", "sbc_summary_corrected.csv"))
bfile <- NA_character_
for (cand in sbc_candidates) {
  if (!file.exists(cand)) next
  hdr <- names(utils::read.csv(cand, nrows = 1, stringsAsFactors = FALSE))
  if (all(SBC_REQ %in% hdr)) { bfile <- cand; break }
  message("SKIP: ", basename(cand), " lacks ",
          paste(setdiff(SBC_REQ, hdr), collapse = ", "),
          " -- pre-fix schema, not reportable.")
}
if (is.na(bfile)) {
  message("SKIP: no SBC summary with a complete schema. Re-run 10_sbc.R.")
} else {
  B <- utils::read.csv(bfile, stringsAsFactors = FALSE)
  ## Use [[ ]] throughout, never $. `$` on a data.frame does PARTIAL matching, so
  ## a missing column can silently resolve to a different one (B$p -> B$param),
  ## which is how this block failed before: as.integer(B$df[i]) returned
  ## integer(0) and sprintf() then returned character(0).
  col <- function(d, nm) if (nm %in% names(d)) d[[nm]] else rep(NA, nrow(d))
  if ("dgm_tag" %in% names(B) && !is.na(TAG) && TAG %in% B[["dgm_tag"]])
    B <- B[B[["dgm_tag"]] == TAG, ]
  B <- B[order(as.numeric(col(B, "p"))), , drop = FALSE]
  stopifnot(nrow(B) > 0)
  ph <- as.numeric(col(B, "p_holm"))
  body <- vapply(seq_len(nrow(B)), function(i) sprintf("%s & %s & %s & %s & %s & %s \\\\",
    if (B[["param"]][i] %in% names(PLAB)) PLAB[[B[["param"]][i]]]
      else paste0("\\texttt{", gsub("_", "\\\\_", B[["param"]][i]), "}"),
    fmt(col(B, "nsbc")[i], 0), fmt(col(B, "chisq")[i], 2), fmt(col(B, "df")[i], 0),
    fmt(col(B, "p")[i], 3), fmt(ph[i], 3)), character(1))
  tex_write(c(
    "% generated by 14_make_tables.R -- do not edit by hand",
    "\\begin{table}[!htbp]", "\\centering",
    "\\caption{Simulation-based calibration. Ranks of the true value among",
    sprintf(" $L=%s$ thinned posterior draws are Uniform$\\{0,\\dots,L\\}$ under correct",
            if ("L" %in% names(B)) B[["L"]][1] else "99"),
    " calibration. The $\\chi^2$ statistic uses exact bin probabilities: the",
    " $L+1$ integer ranks divide into equiprobable bins only when the number of",
    " bins divides $L+1$, and assuming equal expected counts otherwise inflates",
    " the statistic. $p_{\\mathrm{Holm}}$ adjusts for testing all parameters",
    " simultaneously.}",
    "\\label{tab:sbc}", "\\small",
    "\\begin{tabular}{lccccc}", "\\toprule",
    "Parameter & SBC draws & $\\chi^2$ & df & $p$ & $p_{\\mathrm{Holm}}$ \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular}", "\\end{table}"),
    file.path(TABD, "tab_sbc.tex"))
  note("sbc source       : ", bfile)
}

## ---- provenance ------------------------------------------------------------
note("generated        : ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
note("R                : ", R.version.string)
tex_write(PROV, file.path(TABD, "tables_provenance.txt"))
message("\nTables -> ", TABD)
message("\\input them from the manuscript; do not edit the .tex files by hand.")
