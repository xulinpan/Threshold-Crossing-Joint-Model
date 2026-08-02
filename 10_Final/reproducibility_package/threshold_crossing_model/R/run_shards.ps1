# Launches N sharded worker processes of 06_simulation_redesign.R.
# Run from the R/ folder:  ./run_shards.ps1
$N = 4                      # number of shards (4 x 4 chains = 16 cores; use 5 for 20)
for ($i = 0; $i -lt $N; $i++) {
  Start-Process powershell -ArgumentList @(
    "-NoExit","-Command",
    "`$env:N_SHARDS=$N; `$env:SHARD_ID=$i; Rscript 06_simulation_redesign.R"
  ) -WorkingDirectory (Get-Location)
  Start-Sleep -Seconds 4    # stagger starts so the compiled model is reused, not rebuilt
}
Write-Host "Launched $N shards. Aggregate anytime with: Rscript 07_aggregate_sim.R"
