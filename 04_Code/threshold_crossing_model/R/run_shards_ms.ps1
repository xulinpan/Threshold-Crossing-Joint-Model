$N = 4                      # shards (4 x 4 chains = 16 cores; set 5 for all 20)
for ($i = 0; $i -lt $N; $i++) {
  Start-Process powershell -ArgumentList @("-NoExit","-Command",
    "`$env:N_SHARDS=$N; `$env:SHARD_ID=$i; Rscript 08_simulation_multistate.R") -WorkingDirectory (Get-Location)
  Start-Sleep -Seconds 4 }
Write-Host "Launched $N MS shards. Aggregate anytime: Rscript 09_aggregate_ms.R"
