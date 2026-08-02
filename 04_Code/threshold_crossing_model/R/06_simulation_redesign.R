## =============================================================================
## 06_simulation_redesign.R  —  DESKTOP / OVERNIGHT edition (SHARDED)
## ADEMP simulation study for Reviewer point 4 (Statistics in Medicine).
## Evaluates the PROPOSED BAYESIAN estimator with FULL HMC (cmdstanr).
##
## SPEEDUP: the queue is split into N_SHARDS disjoint pieces so you can run
## several R processes at once and use all your cores. Each task has a FIXED
## shard assignment, so shards never collide and never repeat work. Checkpoints
## are shared, so any process resumes the study after a stop.
##
## HOW TO RUN (recommended: use the launcher, 4 shards x 4 chains = 16 cores):
##   PowerShell, from this R/ folder:
##     ./run_shards.ps1
##   or manually launch 4 windows, each:
##     $env:N_SHARDS=4; $env:SHARD_ID=0; Rscript 06_simulation_redesign.R   # 0,1,2,3
##   Single-process fallback (no sharding): just Rscript 06_simulation_redesign.R
##
## Aggregate current results anytime with: Rscript 07_aggregate_sim.R
## REQUIREMENTS: R>=4.1, cmdstanr + CmdStan>=2.30, posterior, dplyr.
##
## CHANGED 2026-07-29
##  (1) The data-generating mechanism moved to dgm_common.R, which is now the
##      only copy (06/10/11 previously each had their own and they had drifted).
##  (2) The "informative monitoring" cell did not implement informative
##      monitoring. It duplicated ~30% of visit times INDEPENDENTLY of the
##      latent trajectory, which also produced zero-length at-risk intervals
##      rescued by pmax(., 1e-6). It is now an NHPP whose intensity depends on
##      the patient's own latent MRD value, so the cell tests what it claims.
##  (3) Bias/coverage are no longer reported for parameters with no generative
##      counterpart under a departure (see truth_for_scenario in dgm_common.R).
##      Previously beta2 coverage was published as 0.000 under the monotone
##      trajectory purely because the truth had no quadratic term, and beta1's
##      bias was measured against -3.56 when the generative slope was -4.06.
##  (4) Checkpoints live in checkpoints_<dgm_tag>/ so fits made under different
##      mechanisms can never be pooled. Legacy fits stay in checkpoints/.
##  (5) The scenario grid is written to outputs/sim_redesign/scenario_grid.csv so
##      the summary can be read without decoding scen_id against this script.
##  (6) Replicate counts are explicit and reproduce what is reported: the old
##      CFG said n_rep_central = 67 while 200 (n=87) and 131 (n=150) checkpoints
##      existed, so the committed config did not regenerate the summary.
## =============================================================================

suppressMessages({ library(cmdstanr); library(posterior); library(dplyr) })

.sim_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
source(file.path(.sim_dir, "dgm_common.R"))

## ---- 0. Configuration ------------------------------------------------------
N_SHARDS <- as.integer(Sys.getenv("N_SHARDS", "1"))
SHARD_ID <- as.integer(Sys.getenv("SHARD_ID", "0"))   # 0 .. N_SHARDS-1
intenv <- function(k, d) as.integer(Sys.getenv(k, as.character(d)))
CFG <- list(
  TIME_BUDGET_HOURS = 10,
  ## Monte Carlo SE on a coverage near 0.95 is sqrt(.95*.05/nrep): 0.015 at 200
  ## reps, 0.022 at 100, 0.040 at 30. 30 reps cannot support a coverage claim to
  ## two decimals, so the misspecification cells default to 100.
  n_rep_central = intenv("N_REP_CENTRAL", 200L),
  n_rep_n300    = intenv("N_REP_N300",    100L),
  n_rep_misspec = intenv("N_REP_MISSPEC", 100L),
  chains = 4L, parallel_chains = 4L,     # 4 cores per fit; 4 shards -> 16 cores
  iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,
  seed = 20260711L,
  out_dir = file.path(.sim_dir, "..", "outputs", "sim_redesign"))
CFG$ckpt_dir <- file.path(CFG$out_dir, paste0("checkpoints_", dgm_tag()))
dir.create(CFG$ckpt_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(CFG$seed)
## QUICK_TEST=1 : tiny run (a handful of fits) to verify the harness end-to-end.
## QUICK_TEST forces a single shard: N_SHARDS/SHARD_ID persist in a PowerShell
## session once a sharded run has been launched, so a smoke test from that same
## window would otherwise run only its own slice and, unless SHARD_ID were 0,
## skip the aggregation step that the smoke test exists to verify.
if (Sys.getenv("QUICK_TEST", "0") == "1") {
  CFG$n_rep_central <- 2L; CFG$n_rep_n300 <- 1L; CFG$n_rep_misspec <- 1L
  if (N_SHARDS != 1L || SHARD_ID != 0L)
    message(sprintf("QUICK_TEST: ignoring inherited N_SHARDS=%d SHARD_ID=%d; running single-shard.",
                    N_SHARDS, SHARD_ID))
  N_SHARDS <- 1L; SHARD_ID <- 0L
  message("QUICK_TEST mode: reps reduced to a handful for a smoke test.") }
## Priors-as-data version of the model. Passing PRIOR_SETS$current reproduces
## interval_hazard_joint.stan exactly, so PRIOR_SET=current is not a change of
## model -- it is the same likelihood and the same priors, expressed as data.
STAN_INTERVAL <- normalizePath(file.path(.sim_dir, "..", "stan",
                                         "interval_hazard_joint_pp.stan"))
dgm_describe()
message(sprintf("Shard %d of %d | detected %d cores | %d parallel chains/fit.",
                SHARD_ID, N_SHARDS, parallel::detectCores(), CFG$parallel_chains))
message("Checkpoints: ", CFG$ckpt_dir)

## ---- 1. Scenario grid ------------------------------------------------------
## Execution order (scen_id) is preserved from the original script so scen_id
## keeps meaning across runs: n=87 and n=150 first, then the misspecification
## cells, then the slow n=300 last.
cell <- function(label, trajectory, floor_kind, err_family, re_family,
                 informative, n, tag, nrep)
  data.frame(label = label, trajectory = trajectory, floor_kind = floor_kind,
             err_family = err_family, re_family = re_family,
             informative = informative, n = n, tag = tag, nrep = nrep,
             stringsAsFactors = FALSE)

grid <- rbind(
  cell("Correct specification, n=87",  "quadratic",  "fixed",         "gaussian","gaussian",FALSE, 87L,"central",CFG$n_rep_central),
  cell("Correct specification, n=150", "quadratic",  "fixed",         "gaussian","gaussian",FALSE,150L,"central",CFG$n_rep_central),
  cell("Trajectory: monotone",         "monotone",   "fixed",         "gaussian","gaussian",FALSE,150L,"misspec",CFG$n_rep_misspec),
  cell("Trajectory: exponential",      "exponential","fixed",         "gaussian","gaussian",FALSE,150L,"misspec",CFG$n_rep_misspec),
  cell("Heterogeneous assay floor",    "quadratic",  "heterogeneous", "gaussian","gaussian",FALSE,150L,"misspec",CFG$n_rep_misspec),
  cell("Heavy-tailed (t3) errors/RE",  "quadratic",  "fixed",         "t3",      "t3",      FALSE,150L,"misspec",CFG$n_rep_misspec),
  cell("Informative monitoring",       "quadratic",  "fixed",         "gaussian","gaussian",TRUE, 150L,"misspec",CFG$n_rep_misspec),
  cell("Correct specification, n=300", "quadratic",  "fixed",         "gaussian","gaussian",FALSE,300L,"central",CFG$n_rep_n300)
)
grid$scen_id <- seq_len(nrow(grid))
## record which parameters have no generative truth in each cell
grid$undefined_params <- vapply(seq_len(nrow(grid)), function(i)
  paste(attr(truth_for_scenario(grid[i, ]), "undefined"), collapse = ";"),
  character(1))
if (SHARD_ID == 0L) {
  dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
  write.csv(grid, file.path(CFG$out_dir, "scenario_grid.csv"), row.names = FALSE)
}

queue <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s) {
  g <- grid[s, ]
  cbind(g[rep(1, g$nrep), setdiff(names(g), "nrep")], rep = seq_len(g$nrep),
        row.names = NULL) }))
queue$task_id <- sprintf("s%02d_%s_n%d_r%03d", queue$scen_id, queue$tag, queue$n, queue$rep)
queue <- queue[order(queue$scen_id, queue$rep), ]                 # deterministic order
queue$shard <- (seq_len(nrow(queue)) - 1L) %% N_SHARDS            # FIXED assignment

## ---- 2. Select this shard's outstanding tasks ------------------------------
done_ids <- sub("\\.rds$", "", list.files(CFG$ckpt_dir, pattern = "\\.rds$"))
todo <- queue[queue$shard == SHARD_ID & !(queue$task_id %in% done_ids), ]
message(sprintf("Queue total %d | this shard %d tasks | %d remaining after resume.",
                nrow(queue), sum(queue$shard == SHARD_ID), nrow(todo)))
mod <- cmdstan_model(STAN_INTERVAL)

## ---- 3. Main loop: checkpointed + time-budgeted ----------------------------
start <- Sys.time(); n_done_now <- 0L
for (r in seq_len(nrow(todo))) {
  if (as.numeric(difftime(Sys.time(), start, units = "hours")) > CFG$TIME_BUDGET_HOURS) {
    message("Time budget reached; stopping cleanly. Re-run to resume."); break }
  task <- as.list(todo[r, ]); f <- file.path(CFG$ckpt_dir, paste0(task$task_id, ".rds"))
  if (file.exists(f)) next
  set.seed(CFG$seed + task$scen_id * 10000L + task$rep)
  out <- tryCatch({
    ds <- simulate_dataset(task$n, task)
    ft <- fit_hmc(mod, c(make_stan_data(ds), stan_priors()), CFG)
    tr <- truth_for_scenario(task)
    m  <- do.call(rbind, lapply(PARAMS, function(p) metrics_one(ft$draws, p, tr[[p]])))
    m$ndiv <- ft$ndiv; m$ebfmi <- ft$ebfmi; m$rhat_max <- ft$rhat_max
    m$ess_min <- ft$ess_min; m$ess_tail_min <- ft$ess_tail_min
    m$rhat_max_all <- ft$rhat_max_all; m$ess_min_all <- ft$ess_min_all
    m$prior_set <- PRIOR_SET
    m$fit_bad <- ft$fit_bad; m$fit_marginal <- ft$fit_marginal
    cbind(task_id = task$task_id, scen_id = task$scen_id, label = task$label,
          tag = task$tag, n = task$n, rep = task$rep, dgm_tag = dgm_tag(), m)
  }, error = function(e) data.frame(
      task_id = task$task_id, scen_id = task$scen_id, label = task$label,
      tag = task$tag, n = task$n, rep = task$rep, dgm_tag = dgm_tag(),
      param = NA, fit_bad = 1L, error_msg = conditionMessage(e)))
  saveRDS(out, f); n_done_now <- n_done_now + 1L
  if (n_done_now %% 10 == 0) message(sprintf("  [shard %d] ...%d fits this session (%.1f h)",
      SHARD_ID, n_done_now, as.numeric(difftime(Sys.time(), start, units = "hours")))) }
message(sprintf("Shard %d finished this session: %d new fits.", SHARD_ID, n_done_now))

## ---- 4. Aggregate (only shard 0 writes, to avoid concurrent-write clashes) --
if (SHARD_ID == 0L) source(file.path(.sim_dir, "07_aggregate_sim.R"), local = TRUE)
