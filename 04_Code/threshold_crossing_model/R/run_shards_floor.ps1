# Launches N sharded workers of 16_floor_vs_exact.R (regenerates Table 10).
# Run from the R/ folder:  ./run_shards_floor.ps1
#
# Fits every simulated dataset TWICE -- identical data, differing only in whether
# assay-floor observations enter the likelihood left-censored or as exact values.
# 100 replicates x 3 sample sizes x 2 fits = 600 fits.

$ErrorActionPreference = "Stop"

$N = 4
$DGM_VERSION     = "v2_cohort"
$TRUTH_SOURCE    = "posterior"
$EVENT_MECHANISM = "hazard"
$PRIOR_SET       = "current"

$STAN = Join-Path (Get-Location) "..\stan\interval_hazard_joint_pp.stan"

# Clear anything inherited from a previous run in this window, which would
# otherwise silently shrink the study or misassign shards.
foreach ($v in @("QUICK_TEST","N_SHARDS","SHARD_ID","N_REP_FLOOR")) {
  if (Test-Path "env:$v") {
    Write-Host "  clearing inherited `$env:$v = $((Get-Item "env:$v").Value)" -ForegroundColor Yellow
    Remove-Item "env:$v"
  }
}

if ($TRUTH_SOURCE -eq "posterior") {
  $psum = Join-Path (Get-Location) "..\outputs\posterior_summary_interval.csv"
  if (-not (Test-Path $psum)) {
    Write-Host "ABORT: TRUTH_SOURCE=posterior needs $psum" -ForegroundColor Red; exit 1
  }
}

Write-Host "Compiling Stan model once before launching..." -ForegroundColor Cyan
Rscript -e "invisible(cmdstanr::cmdstan_model('$($STAN -replace '\\','/')'))"
if ($LASTEXITCODE -ne 0) { Write-Host "ABORT: compile failed." -ForegroundColor Red; exit 1 }
Write-Host "Compile OK." -ForegroundColor Green

for ($i = 0; $i -lt $N; $i++) {
  Start-Process powershell -ArgumentList @(
    "-NoExit","-Command",
    "`$env:DGM_VERSION='$DGM_VERSION'; `$env:TRUTH_SOURCE='$TRUTH_SOURCE'; `$env:EVENT_MECHANISM='$EVENT_MECHANISM'; `$env:PRIOR_SET='$PRIOR_SET'; `$env:N_SHARDS=$N; `$env:SHARD_ID=$i; Remove-Item env:QUICK_TEST -ErrorAction SilentlyContinue; Rscript 16_floor_vs_exact.R"
  ) -WorkingDirectory (Get-Location)
  Start-Sleep -Seconds 2
}
Write-Host "Launched $N shards. 600 fits; expect roughly 3-4 hours."
Write-Host "When all four windows finish, aggregate with:"
Write-Host "  `$env:SHARD_ID=0; `$env:N_SHARDS=1; Rscript 16_floor_vs_exact.R"
Write-Host "(it skips every completed fit and only re-aggregates)"
