# Launches N sharded worker processes of 06_simulation_redesign.R.
# Run from the R/ folder:  ./run_shards.ps1
#
# The DGM switches are set EXPLICITLY here rather than relying on defaults, so
# the run is self-documenting and every shard generates data the same way.
# See dgm_common.R for what each switch means.
#   DGM_VERSION     v2_cohort  visit gaps matching the cohort (median 6 months)
#   TRUTH_SOURCE    posterior  truth = fitted posterior means, read from
#                              outputs/posterior_summary_interval.csv
#   EVENT_MECHANISM hazard     events from the likelihood the Stan model fits
# Checkpoints land in outputs/sim_redesign/checkpoints_<DGM_VERSION>_<TRUTH_SOURCE>_<EVENT_MECHANISM>/
# so they can never be pooled with the pre-2026-07-29 legacy fits.

$ErrorActionPreference = "Stop"

$N = 4                      # number of shards (4 x 4 chains = 16 cores; use 5 for 20)
$DGM_VERSION     = "v2_cohort"
$TRUTH_SOURCE    = "posterior"
$EVENT_MECHANISM = "hazard"
# current | wide2 | wide4 | vague  (see PRIOR_SETS in dgm_common.R)
# "current" reproduces the published priors and keeps the existing checkpoint
# directory; anything else appends the set name to the tag and starts a fresh
# directory, so prior specifications can never be pooled into one summary.
$PRIOR_SET       = "current"

$STAN = Join-Path (Get-Location) "..\stan\interval_hazard_joint.stan"

# ---------------------------------------------------------------------------
# Pre-flight. Each check below has cost 4 wasted processes or a silently
# invalid study if it is skipped.
# ---------------------------------------------------------------------------

# 1. QUICK_TEST leaks into children. Start-Process inherits the parent
#    environment, so a QUICK_TEST=1 left over from a smoke test would run the
#    whole study at 2 replicates per cell and produce output indistinguishable
#    from a real run. Same for the replicate overrides. Clear them explicitly.
foreach ($v in @("QUICK_TEST","N_REP_CENTRAL","N_REP_N300","N_REP_MISSPEC")) {
  if (Test-Path "env:$v") {
    Write-Host "  clearing inherited `$env:$v = $((Get-Item "env:$v").Value)" -ForegroundColor Yellow
    Remove-Item "env:$v"
  }
}

# 2. TRUTH_SOURCE=posterior reads the fitted posterior summary. If it is absent,
#    every shard dies on start with a stop() from dgm_common.R.
if ($TRUTH_SOURCE -eq "posterior") {
  $psum = Join-Path (Get-Location) "..\outputs\posterior_summary_interval.csv"
  if (-not (Test-Path $psum)) {
    Write-Host "ABORT: TRUTH_SOURCE=posterior needs $psum" -ForegroundColor Red
    Write-Host "       Run 03_numeric_results.R first, or set TRUTH_SOURCE='legacy'." -ForegroundColor Red
    exit 1
  }
}

# 3. Compile the Stan model ONCE, before any shard starts. Otherwise shard 0
#    begins a ~60 s compile and shards 1..3 launch into the same output path and
#    race on the same executable.
Write-Host "Compiling Stan model (once, before launching shards)..." -ForegroundColor Cyan
Rscript -e "invisible(cmdstanr::cmdstan_model('$($STAN -replace '\\','/')'))"
if ($LASTEXITCODE -ne 0) {
  Write-Host "ABORT: Stan model failed to compile. Nothing launched." -ForegroundColor Red
  exit 1
}
Write-Host "Compile OK." -ForegroundColor Green

# 4. Warn if oversubscribed. N shards x 4 chains each.
$cores = [Environment]::ProcessorCount
if ($N * 4 -gt $cores) {
  Write-Host "WARNING: $N shards x 4 chains = $($N*4) processes on $cores logical cores." -ForegroundColor Yellow
  Write-Host "         Oversubscription mostly costs wall time, not correctness." -ForegroundColor Yellow
}

# 5. Confirm the mechanism looks like the cohort before spending hours on it.
Write-Host ""
Write-Host "Reminder: sanity-check the mechanism first if you have not already:" -ForegroundColor Cyan
Write-Host "  Rscript 13_dgm_sanity_check.R" -ForegroundColor Cyan
Write-Host ""

# ---------------------------------------------------------------------------
# Launch
# ---------------------------------------------------------------------------
for ($i = 0; $i -lt $N; $i++) {
  Start-Process powershell -ArgumentList @(
    "-NoExit","-Command",
    "`$env:DGM_VERSION='$DGM_VERSION'; `$env:TRUTH_SOURCE='$TRUTH_SOURCE'; `$env:EVENT_MECHANISM='$EVENT_MECHANISM'; `$env:PRIOR_SET='$PRIOR_SET'; `$env:N_SHARDS=$N; `$env:SHARD_ID=$i; Remove-Item env:QUICK_TEST -ErrorAction SilentlyContinue; Rscript 06_simulation_redesign.R"
  ) -WorkingDirectory (Get-Location)
  Start-Sleep -Seconds 2
}
Write-Host "Launched $N shards ($DGM_VERSION / $TRUTH_SOURCE / $EVENT_MECHANISM / priors=$PRIOR_SET)."
Write-Host "NOTE: shard 0 aggregates when ITS queue empties, so the summary it"
Write-Host "      writes may be partial. Re-run 07 after all shards finish:"
Write-Host "  Rscript 07_aggregate_sim.R"
