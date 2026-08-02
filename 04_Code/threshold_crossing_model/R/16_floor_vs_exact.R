## =============================================================================
## 16_floor_vs_exact.R  —  the cost of treating assay-floor values as exact
##
## REGENERATES TABLE 10 (tab:sim-bias). The published Table 10 reports a
## floor-censored bias in tau_b1 of -0.10 / -0.14 / -0.08 at n = 87 / 150 / 300.
## The only simulation output in this repository, outputs/simulation_summary.csv,
## reports +9.01 / +11.70 / +16.11 for the same estimator and cells -- opposite
## sign, two orders of magnitude apart -- and no script here reproduces the
## published values. That output came from 04_simulation_study.R, which estimates
## by $optimize(), i.e. the mode of the joint posterior over hyperparameters AND
## the latent random effects; the joint mode is not a consistent estimator of
## variance components, and its tau_b1 bias correspondingly grows with n. Both
## sets of numbers are therefore unusable and the claim that floor censoring is
## "a necessity rather than a refinement" currently has no supporting evidence.
##
## This script re-establishes the contrast properly: same datasets, same HMC
## estimator, fitted twice, differing only in whether floor-coded observations
## enter as left-censored (floor_ind as generated) or as exact values at cF.
## Pairing the two fits on the same data removes between-replicate variation
## from the comparison.
##
## OUTPUT
##   outputs/floor_vs_exact/floor_vs_exact_raw.csv
##   outputs/floor_vs_exact/floor_vs_exact_summary.csv
##   06_Tables/table_10_sim_bias.tex        (written by 14_make_tables.R)
##
## RUN (from the R/ folder). Sharding mirrors 06; 2 fits per replicate.
##   $env:QUICK_TEST=1; Rscript 16_floor_vs_exact.R          # smoke test
##   $env:N_SHARDS=4; $env:SHARD_ID=0; Rscript 16_floor_vs_exact.R
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior) })

.fx_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
source(file.path(.fx_dir, "dgm_common.R"))
dgm_describe()

N_SHARDS <- as.integer(Sys.getenv("N_SHARDS", "1"))
SHARD_ID <- as.integer(Sys.getenv("SHARD_ID", "0"))
CFG <- list(
  TIME_BUDGET_HOURS = 10,
  n_rep = as.integer(Sys.getenv("N_REP_FLOOR", "100")),
  sample_sizes = c(87L, 150L, 300L),
  chains = 4L, parallel_chains = 4L,
  iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,
  seed = 20260731L,
  out_dir = file.path(.fx_dir, "..", "outputs", "floor_vs_exact"))
CFG$ckpt_dir <- file.path(CFG$out_dir, paste0("checkpoints_", dgm_tag()))
dir.create(CFG$ckpt_dir, recursive = TRUE, showWarnings = FALSE)
## QUICK_TEST forces a single shard. N_SHARDS/SHARD_ID persist in a PowerShell
## session once a sharded run has been launched, so a smoke test invoked from
## that same window would otherwise run only its own slice and -- if SHARD_ID
## were not 0 -- skip aggregation, which is precisely the path being tested.
if (Sys.getenv("QUICK_TEST", "0") == "1") {
  CFG$n_rep <- 2L
  if (N_SHARDS != 1L || SHARD_ID != 0L)
    message(sprintf("QUICK_TEST: ignoring inherited N_SHARDS=%d SHARD_ID=%d; running single-shard.",
                    N_SHARDS, SHARD_ID))
  N_SHARDS <- 1L; SHARD_ID <- 0L
  message("QUICK_TEST: 2 reps.")
}

STAN <- normalizePath(file.path(.fx_dir, "..", "stan", "interval_hazard_joint_pp.stan"))
mod  <- cmdstan_model(STAN)
CENTRAL <- list(trajectory = "quadratic", floor_kind = "fixed",
                err_family = "gaussian", re_family = "gaussian", informative = FALSE)
## the three parameters the floor contrast is about
REPORT <- c("beta1", "tau1", "sigma_y")

## ---- queue ------------------------------------------------------------------
queue <- do.call(rbind, lapply(CFG$sample_sizes, function(nn)
  data.frame(n = nn, rep = seq_len(CFG$n_rep))))
queue$task_id <- sprintf("n%03d_r%03d", queue$n, queue$rep)
queue <- queue[order(queue$n, queue$rep), ]
queue$shard <- (seq_len(nrow(queue)) - 1L) %% N_SHARDS
done <- sub("\\.rds$", "", list.files(CFG$ckpt_dir, pattern = "\\.rds$"))
todo <- queue[queue$shard == SHARD_ID & !(queue$task_id %in% done), ]
message(sprintf("Queue %d | shard %d has %d | %d remaining.",
                nrow(queue), SHARD_ID, sum(queue$shard == SHARD_ID), nrow(todo)))

## ---- main loop --------------------------------------------------------------
start <- Sys.time(); ndone <- 0L
for (i in seq_len(nrow(todo))) {
  if (as.numeric(difftime(Sys.time(), start, units = "hours")) > CFG$TIME_BUDGET_HOURS) {
    message("Time budget reached; re-run to resume."); break }
  tk <- as.list(todo[i, ]); f <- file.path(CFG$ckpt_dir, paste0(tk$task_id, ".rds"))
  if (file.exists(f)) next
  set.seed(CFG$seed + 1000L * tk$n + tk$rep)
  out <- tryCatch({
    ds <- simulate_dataset(tk$n, CENTRAL)
    sd_cens <- make_stan_data(ds)
    ## The ONLY difference: floor-coded observations declared exact rather than
    ## left-censored. The y values themselves are identical (both are cF), so
    ## this isolates the likelihood contribution, not the data.
    sd_exct <- sd_cens; sd_exct$floor_ind <- rep(0L, sd_cens$Nobs)

    one <- function(sdat, tag) {
      ft <- fit_hmc(mod, c(sdat, stan_priors()), CFG)
      do.call(rbind, lapply(REPORT, function(p) {
        x <- as.numeric(ft$draws[[p]]); qi <- quantile(x, c(0.025, 0.975))
        data.frame(estimator = tag, param = p, post_mean = mean(x),
                   truth = TRUTH[[p]], bias = mean(x) - TRUTH[[p]],
                   covered = as.integer(TRUTH[[p]] >= qi[1] & TRUTH[[p]] <= qi[2]),
                   fit_bad = ft$fit_bad, rhat_max = ft$rhat_max, ndiv = ft$ndiv)
      }))
    }
    m <- rbind(one(sd_cens, "censored"), one(sd_exct, "exact"))
    cbind(task_id = tk$task_id, n = tk$n, rep = tk$rep, dgm_tag = dgm_tag(), m)
  }, error = function(e) data.frame(
      task_id = tk$task_id, n = tk$n, rep = tk$rep, dgm_tag = dgm_tag(),
      estimator = NA, param = NA, fit_bad = 1L, error_msg = conditionMessage(e)))
  saveRDS(out, f); ndone <- ndone + 1L
  if (ndone %% 5 == 0) message(sprintf("  ...%d replicates (%.1f h)", ndone,
      as.numeric(difftime(Sys.time(), start, units = "hours"))))
}
message(sprintf("Shard %d: %d new replicates.", SHARD_ID, ndone))

## ---- aggregate (shard 0) ----------------------------------------------------
if (SHARD_ID == 0L) {
  fs <- list.files(CFG$ckpt_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(fs)) {
    R <- do.call(rbind, lapply(fs, readRDS))
    R <- R[!is.na(R$param) & (is.na(R$fit_bad) | R$fit_bad == 0), ]
    write.csv(R, file.path(CFG$out_dir, "floor_vs_exact_raw.csv"), row.names = FALSE)
    ag <- do.call(rbind, lapply(split(R, list(R$n, R$estimator, R$param), drop = TRUE),
      function(d) {
        k <- nrow(d)
        data.frame(n = d$n[1], estimator = d$estimator[1], param = d$param[1],
                   nrep = k, truth = d$truth[1],
                   bias = mean(d$bias), bias_mcse = sd(d$bias)/sqrt(k),
                   rmse = sqrt(mean(d$bias^2)),
                   coverage = mean(d$covered), row.names = NULL) }))
    ag <- ag[order(ag$n, ag$param, ag$estimator), ]
    write.csv(ag, file.path(CFG$out_dir, "floor_vs_exact_summary.csv"), row.names = FALSE)
    print(ag, row.names = FALSE, digits = 3)
    message("\n-> floor_vs_exact_summary.csv. Run 14_make_tables.R to build Table 10.")
  }
}
