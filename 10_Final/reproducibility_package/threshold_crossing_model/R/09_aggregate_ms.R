suppressMessages(library(dplyr))
ck<-file.path("..","outputs","sim_ms","checkpoints"); od<-file.path("..","outputs","sim_ms")
fs<-list.files(ck,pattern="\\.rds$",full.names=TRUE)
if (length(fs)==0){message("No MS checkpoints yet.")} else {
  res<-dplyr::bind_rows(lapply(fs,function(f) tryCatch(readRDS(f),error=function(e) NULL))); res<-res[!is.na(res$param)&!is.na(res$bias),]
  agg<-res%>%group_by(scen_id,tag,n,relapse,cD,param)%>%
    summarise(nrep=n(),bias=mean(bias),bias_mcse=sd(bias)/sqrt(n()),
              coverage=mean(covered),cov_mcse=sqrt(mean(covered)*(1-mean(covered))/n()),
              ci_width=mean(ci_width),onset_rate=mean(onset_rate),
              relapse_rate=mean(relapse_rate,na.rm=TRUE),comp_fail=mean(comp_fail),.groups="drop")
  write.csv(res,file.path(od,"sim_ms_raw.csv"),row.names=FALSE)
  write.csv(agg,file.path(od,"sim_ms_summary.csv"),row.names=FALSE)
  message(sprintf("Aggregated %d MS fits -> sim_ms_summary.csv",length(fs))) }
