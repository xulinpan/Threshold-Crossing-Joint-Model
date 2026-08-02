# Launches N sharded workers of 17_calibration_study.R.
# Regenerates Table 12 and provides the first simulation of the
# latent threshold-crossing model.  Run from the R/ folder: ./run_shards_calib.ps1
#
# 5 cells x 100 replicates = 500 fits. The threshold arm is the slower one
# (soft-min over a G=60 grid per patient) and is also the cell of interest for
# divergent transitions.

$ErrorActionPreference = "Stop"
$N = 4
$DGM_VERSION     = "v2_cohort"
$TRUTH_SOURCE    = "posterior"
$EVENT_MECHANISM = "hazard"
$PRIOR_SET       = "current"

foreach ($v in @("QUICK_TEST","N_SHARDS","SHARD_ID","N_REP_CALIB")) {
  if (Test-Path "env:$v") {
    Write-Host "  clearing inherited `$env:$v = $((Get-Item "env:$v").Value)" -ForegroundColor Yellow
    Remove-Item "env:$v"
  }
}

$psum = Join-Path (Get-Location) "..\outputs\posterior_summary_interval.csv"
if (-not (Test-Path $psum)) { Write-Host "ABORT: missing $psum" -ForegroundColor Red; exit 1 }

Write-Host "Compiling both models once before launching..." -ForegroundColor Cyan
foreach ($s in @("interval_hazard_joint_pp.stan","threshold_crossing_joint.stan")) {
  $p = (Join-Path (Get-Location) "..\stan\$s") -replace '\\','/'
  Rscript -e "invisible(cmdstanr::cmdstan_model('$p'))"
  if ($LASTEXITCODE -ne 0) { Write-Host "ABORT: $s failed to compile." -ForegroundColor Red; exit 1 }
}
Write-Host "Compile OK." -ForegroundColor Green

for ($i = 0; $i -lt $N; $i++) {
  Start-Process powershell -ArgumentList @(
    "-NoExit","-Command",
    "`$env:DGM_VERSION='$DGM_VERSION'; `$env:TRUTH_SOURCE='$TRUTH_SOURCE'; `$env:EVENT_MECHANISM='$EVENT_MECHANISM'; `$env:PRIOR_SET='$PRIOR_SET'; `$env:N_SHARDS=$N; `$env:SHARD_ID=$i; Remove-Item env:QUICK_TEST -ErrorAction SilentlyContinue; Rscript 17_calibration_study.R"
  ) -WorkingDirectory (Get-Location)
  Start-Sleep -Seconds 2
}
Write-Host "Launched $N shards. 500 fits."
Write-Host "Aggregate when finished:  `$env:SHARD_ID=0; `$env:N_SHARDS=1; Rscript 17_calibration_study.R"
