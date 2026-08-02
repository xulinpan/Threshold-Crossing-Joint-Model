## Robust aggregator: tolerates error-rows (different columns) and partial files.
suppressMessages(library(dplyr))
ckpt <- file.path("..","outputs","sim_redesign","checkpoints")
outd <- file.path("..","outputs","sim_redesign")
fs <- list.files(ckpt, pattern="\\.rds$", full.names=TRUE)
if (length(fs)==0){ message("No checkpoints yet."); } else {
  rows <- lapply(fs, function(f) tryCatch(readRDS(f), error=function(e) NULL))
  res  <- dplyr::bind_rows(rows)                          # fills missing cols with NA
  ok   <- res[!is.na(res$param) & !is.na(res$bias), ]     # keep successful metric rows
  nfit <- length(unique(res$task_id))
  nerr <- length(unique(res$task_id[is.na(res$bias)]))
  agg <- ok %>% group_by(scen_id, tag, n, param) %>%
    summarise(nrep=n(), bias=mean(bias), bias_mcse=sd(bias)/sqrt(n()),
              coverage=mean(covered), cov_mcse=sqrt(mean(covered)*(1-mean(covered))/n()),
              ci_width=mean(ci_width), comp_fail=mean(comp_fail), .groups="drop")
  write.csv(ok,  file.path(outd,"sim_redesign_raw.csv"), row.names=FALSE)
  write.csv(agg, file.path(outd,"sim_redesign_summary.csv"), row.names=FALSE)
  message(sprintf("Aggregated %d fits (%d errored) across %d scenarios -> sim_redesign_summary.csv",
                  nfit, nerr, length(unique(ok$scen_id)))) }
