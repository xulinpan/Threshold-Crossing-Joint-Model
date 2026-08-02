## =============================================================================
## 07_aggregate_sim.R  —  robust aggregator for 06_simulation_redesign.R
## Tolerates error-rows (different columns) and partial checkpoint sets.
##
## FIXED 2026-07-29
##  (1) bias_mcse was NA in every row. It was computed as
##          summarise(bias = mean(bias), bias_mcse = sd(bias)/sqrt(n()))
##      and dplyr evaluates a summarise() sequentially, so `bias` was already
##      the scalar mean by the time bias_mcse ran -> sd(scalar) = NA. The MCse
##      is now computed FIRST, from the replicate-level values.
##  (2) Scenario labels are carried through, so the summary no longer has to be
##      decoded against the grid construction inside 06.
##  (3) The single comp_fail flag is split. R-hat > 1.05 / E-BFMI < 0.2 /
##      divergences (fit_bad) are genuinely non-converged and are EXCLUDED from
##      bias and coverage; 1.01 < R-hat <= 1.05 (fit_marginal) is benign at
##      4 x 500 draws and is only reported. Previously a handful of fits with
##      R-hat up to 1.78 were being averaged into the headline results while
##      100% of the n=300 cell was flagged for R-hat 1.02.
##  (4) Parameters with no generative counterpart under a departure
##      (truth_defined = FALSE) are dropped instead of reporting coverage
##      against a value that does not exist.
## =============================================================================
suppressMessages(library(dplyr))

.agg_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
outd <- file.path(.agg_dir, "..", "outputs", "sim_redesign")

## Pool every checkpoint directory (one per DGM tag) but keep the tag as a column
## so mechanisms are never silently mixed. Legacy runs used a bare "checkpoints".
ckpt_dirs <- list.dirs(outd, recursive = FALSE)
ckpt_dirs <- ckpt_dirs[grepl("^checkpoints", basename(ckpt_dirs))]
if (!length(ckpt_dirs)) { message("No checkpoint directories under ", outd); }
fs <- unlist(lapply(ckpt_dirs, function(d)
  list.files(d, pattern = "\\.rds$", full.names = TRUE)), use.names = FALSE)

if (length(fs) == 0) {
  message("No checkpoints yet.")
} else {
  rows <- lapply(fs, function(f) {
    x <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(x)) return(NULL)
    if (!"dgm_tag" %in% names(x)) {
      tg <- sub("^checkpoints_?", "", basename(dirname(f)))
      x$dgm_tag <- if (nzchar(tg)) tg else "v1_legacy_legacy_hazard"
    }
    x
  })
  res <- dplyr::bind_rows(rows)                  # fills missing cols with NA

  ## ---- back-compatibility with pre-fix checkpoints -------------------------
  if (!"truth_defined" %in% names(res)) res$truth_defined <- TRUE
  res$truth_defined[is.na(res$truth_defined)] <- TRUE
  ## Old schema stored only comp_fail = (ndiv>0 | ebfmi<0.2 | rhat>1.01). Re-derive
  ## the two-level flag wherever it is absent, so a pool of old and new
  ## checkpoints is handled row by row rather than all-or-nothing.
  if (!"fit_bad" %in% names(res))      res$fit_bad      <- NA_integer_
  if (!"fit_marginal" %in% names(res)) res$fit_marginal <- NA_integer_
  need <- is.na(res$fit_bad)
  if (any(need)) {
    res$fit_bad[need] <- as.integer(
      (!is.na(res$ndiv[need])     & res$ndiv[need]  > 0) |
      (!is.na(res$ebfmi[need])    & res$ebfmi[need] < 0.2) |
      (!is.na(res$rhat_max[need]) & res$rhat_max[need] > 1.05))
    res$fit_marginal[need] <- as.integer(!is.na(res$rhat_max[need]) &
                                 res$rhat_max[need] > 1.01 & res$rhat_max[need] <= 1.05)
  }

  ## Apply the undefined-truth rule RETROACTIVELY to checkpoints written before
  ## truth_for_scenario() existed. Those fits stored a bias against a parameter
  ## the DGM never had (e.g. beta2 under a monotone trajectory), so without this
  ## the legacy summary keeps publishing coverage 0.000 for it. Safe because
  ## scen_id ordering and the monotone/exponential forms are unchanged.
  gfile0 <- file.path(outd, "scenario_grid.csv")
  if (file.exists(gfile0)) {
    g0 <- utils::read.csv(gfile0, stringsAsFactors = FALSE)
    if ("undefined_params" %in% names(g0)) {
      for (i in seq_len(nrow(g0))) {
        up <- trimws(strsplit(as.character(g0$undefined_params[i]), ";")[[1]])
        up <- up[nzchar(up)]
        if (length(up))
          res$truth_defined[res$scen_id == g0$scen_id[i] & res$param %in% up] <- FALSE
      }
    }
  }

  ok <- res[!is.na(res$param) & !is.na(res$bias) & res$truth_defined, ]
  nfit <- length(unique(res$task_id))
  nerr <- length(unique(res$task_id[is.na(res$bias) & is.na(res$param)]))

  ## ---- fit-level diagnostic rates (per scenario, one row per fit) ----------
  fitlvl <- res %>%
    distinct(dgm_tag, scen_id, tag, n, task_id, .keep_all = TRUE) %>%
    group_by(dgm_tag, scen_id, tag, n) %>%
    summarise(nfit = n(),
              rate_bad      = mean(fit_bad,      na.rm = TRUE),
              rate_marginal = mean(fit_marginal, na.rm = TRUE),
              max_rhat      = suppressWarnings(max(rhat_max, na.rm = TRUE)),
              min_ebfmi     = suppressWarnings(min(ebfmi,    na.rm = TRUE)),
              tot_div       = sum(ndiv, na.rm = TRUE), .groups = "drop")

  ## ---- performance measures, on converged fits only ------------------------
  good <- ok[is.na(ok$fit_bad) | ok$fit_bad == 0, ]
  agg <- good %>%
    group_by(dgm_tag, scen_id, tag, n, param) %>%
    summarise(nrep      = n(),
              ## MCse FIRST: `bias` must still be the vector of replicate biases
              bias_mcse = sd(bias) / sqrt(n()),
              bias      = mean(bias),
              coverage  = mean(covered),
              cov_mcse  = sqrt(mean(covered) * (1 - mean(covered)) / n()),
              ci_width  = mean(ci_width),
              .groups   = "drop") %>%
    relocate(bias, .before = bias_mcse) %>%
    left_join(fitlvl, by = c("dgm_tag","scen_id","tag","n"))

  ## ---- attach scenario labels if 06 wrote the grid -------------------------
  gfile <- file.path(outd, "scenario_grid.csv")
  if (file.exists(gfile)) {
    g <- utils::read.csv(gfile, stringsAsFactors = FALSE)
    keep <- intersect(c("scen_id","label","trajectory","floor_kind","err_family",
                        "re_family","informative"), names(g))
    agg <- left_join(agg, g[, keep, drop = FALSE], by = "scen_id") %>%
      relocate(any_of("label"), .after = scen_id)
  } else {
    message("NOTE: ", gfile, " not found -- summary has scen_id but no labels. ",
            "Re-run 06 to write it.")
  }

  ## ---- parameters dropped because they have no true value ------------------
  dropped <- res[!is.na(res$param) & !res$truth_defined, ]
  if (nrow(dropped)) {
    dd <- dropped %>% count(scen_id, param, name = "nrep")
    write.csv(dd, file.path(outd, "sim_redesign_undefined_params.csv"), row.names = FALSE)
    message(sprintf("%d (scenario, parameter) cells have no generative truth and were dropped -> sim_redesign_undefined_params.csv",
                    nrow(dd)))
  }

  write.csv(good, file.path(outd, "sim_redesign_raw.csv"),     row.names = FALSE)
  write.csv(agg,  file.path(outd, "sim_redesign_summary.csv"), row.names = FALSE)
  write.csv(fitlvl, file.path(outd, "sim_redesign_diagnostics.csv"), row.names = FALSE)
  nbad <- sum(fitlvl$rate_bad * fitlvl$nfit)
  message(sprintf("Aggregated %d fits (%d errored, %d non-converged and excluded) across %d scenarios.",
                  nfit, nerr, round(nbad), length(unique(agg$scen_id))))
  message("-> sim_redesign_summary.csv | sim_redesign_diagnostics.csv")
}
