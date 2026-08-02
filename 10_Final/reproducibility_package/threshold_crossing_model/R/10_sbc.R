## =============================================================================
## 10_sbc.R  —  Simulation-Based Calibration (SBC) for the joint interval model
## (Reviewer point 4, required item 3). Verifies that the HMC posterior is
## correctly calibrated: draw theta ~ prior, simulate data, fit, and record the
## rank of each true theta among (thinned, ~independent) posterior draws. Ranks
## should be Uniform{0..L}; departures reveal mis-calibration.
##
## Self-contained; reuses interval_hazard_joint.stan. Checkpointed like 06.
## RUN:  $env:QUICK_TEST=1; Rscript 10_sbc.R   (smoke test)  then  Rscript 10_sbc.R
## Aggregate/plot:  Rscript -e "source('10_sbc.R', echo=FALSE)"  (auto at end)
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior) })

CFG <- list(
  N_SBC = 150L,          # number of SBC iterations (>=100 recommended)
  n = 87L,               # cohort size per SBC draw (n=87 keeps fits fast)
  L = 100L,              # thinned posterior draws used for the rank (0..L)
  chains = 4L, parallel_chains = 4L, iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,
  TIME_BUDGET_HOURS = 10, seed = 20260713L,
  out_dir  = file.path("..","outputs","sbc"),
  ckpt_dir = file.path("..","outputs","sbc","checkpoints"))
dir.create(CFG$ckpt_dir, recursive=TRUE, showWarnings=FALSE); set.seed(CFG$seed)
if (Sys.getenv("QUICK_TEST","0")=="1"){ CFG$N_SBC<-3L; message("QUICK_TEST: 3 SBC draws.") }
STAN <- normalizePath(file.path("..","stan","interval_hazard_joint.stan"))
FLOOR <- -5.0
PARAMS <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","gamma0","gamma1","gamma2","alpha")

## ---- draw theta from the EXACT priors of interval_hazard_joint.stan ---------
draw_from_prior <- function() list(
  beta0=rnorm(1,-2.5,2), beta1=rnorm(1,-1,1), beta2=rnorm(1,0,0.5), beta_bm=rnorm(1,0,1),
  sigma_y=rexp(1,1), tau0=rexp(1,1), tau1=rexp(1,1),
  gamma0=rnorm(1,-2,2), gamma1=rnorm(1,0,1), gamma2=rnorm(1,0,1), alpha=rnorm(1,-0.5,0.75))

## ---- simulate one dataset from given parameters (matches the fitted model) --
simulate_from_params <- function(th, n){
  b0<-rnorm(n,0,th$tau0); b1<-rnorm(n,0,th$tau1)
  visit<-function(){ t<-0;v<-numeric(0); while(t<5){t<-t+rexp(1,4); if(t<5)v<-c(v,t)}; sort(unique(round(c(0.25,v),3))) }
  long<-list(); iv<-list()
  for (i in 1:n){
    vt<-visit(); if (length(vt)<2) vt<-c(0.25,0.75); ell<-log1p(vt)
    m<-th$beta0+(th$beta1+b1[i])*ell+th$beta2*ell^2+b0[i]
    bm<-rbinom(length(vt),1,0.9); y<-m+th$beta_bm*bm+rnorm(length(vt),0,th$sigma_y)
    fl<-as.integer(y<=FLOOR); y[fl==1]<-FLOOR
    long[[i]]<-data.frame(pid=i,ell=ell,y=y,floor_ind=fl,bm=bm)
    ts<-head(vt,-1); te<-tail(vt,-1); midl<-log1p(.5*(ts+te)); gapl<-log1p(te-ts); dl<-pmax(te-ts,1e-6)
    mm<-th$beta0+(th$beta1+b1[i])*midl+th$beta2*midl^2+b0[i]
    logh<-th$gamma0+th$gamma1*midl+th$gamma2*gapl+th$alpha*mm
    ev<-rbinom(length(dl),1,1-exp(-exp(logh)*dl)); first<-which(ev==1)[1]
    keep<-if(is.na(first)) seq_along(ts) else seq_len(first)
    e<-rep(0L,length(keep)); if(!is.na(first)) e[first]<-1L
    iv[[i]]<-data.frame(pid=i,t_start=ts[keep],t_end=te[keep],event=e)
  }
  L<-do.call(rbind,long); I<-do.call(rbind,iv)
  list(N=n,Nobs=nrow(L),pid=L$pid,ell=L$ell,y=L$y,floor_ind=L$floor_ind,bm=L$bm,cF=FLOOR,
       Nint=nrow(I),pid_int=I$pid,midlog=log1p(.5*(I$t_start+I$t_end)),
       gaplog=log1p(I$t_end-I$t_start),delta_len=pmax(I$t_end-I$t_start,1e-6),event=I$event)
}

## ---- rank of the true value among L thinned posterior draws ----------------
posterior_rank <- function(draws, th, L){
  out<-integer(length(PARAMS)); names(out)<-PARAMS
  for (p in PARAMS){ x<-as.numeric(draws[[p]]); idx<-round(seq(1,length(x),length.out=L))
    out[p]<-sum(x[idx] < th[[p]]) }                       # rank in 0..L
  out }

mod<-cmdstan_model(STAN); start<-Sys.time()
done<-sub("\\.rds$","",list.files(CFG$ckpt_dir,pattern="\\.rds$"))
for (s in seq_len(CFG$N_SBC)){
  id<-sprintf("sbc_%04d",s); f<-file.path(CFG$ckpt_dir,paste0(id,".rds")); if (id %in% done) next
  if (as.numeric(difftime(Sys.time(),start,units="hours"))>CFG$TIME_BUDGET_HOURS){message("budget reached");break}
  set.seed(CFG$seed+s)
  out<-tryCatch({ th<-draw_from_prior(); sd<-simulate_from_params(th,CFG$n)
    fit<-mod$sample(data=sd,chains=CFG$chains,parallel_chains=CFG$parallel_chains,
      iter_warmup=CFG$iter_warmup,iter_sampling=CFG$iter_sampling,adapt_delta=CFG$adapt_delta,
      max_treedepth=CFG$max_treedepth,refresh=0,show_messages=FALSE,show_exceptions=FALSE)
    dr<-fit$draws(format="draws_df"); r<-posterior_rank(dr,th,CFG$L)
    dg<-fit$diagnostic_summary(quiet=TRUE)
    data.frame(sbc=s,param=names(r),rank=as.integer(r),L=CFG$L,
               ndiv=sum(dg$num_divergent),comp_fail=as.integer(sum(dg$num_divergent)>0))
  }, error=function(e) data.frame(sbc=s,param=NA,rank=NA,L=CFG$L))
  saveRDS(out,f); if (s%%10==0) message(sprintf("  SBC %d/%d",s,CFG$N_SBC))
}

## ---- aggregate: uniformity test per parameter ------------------------------
fs<-list.files(CFG$ckpt_dir,pattern="\\.rds$",full.names=TRUE)
if (length(fs)>0){
  res<-do.call(rbind,lapply(fs,readRDS)); res<-res[!is.na(res$rank),]
  chisq_unif<-function(rk,L){ b<-cut(rk,breaks=seq(-.5,L+.5,length.out=min(20,L+1)+1))
    o<-table(b); e<-rep(mean(o),length(o)); sum((o-e)^2/e) }
  summ<-do.call(rbind,lapply(split(res,res$param),function(d)
    data.frame(param=d$param[1],nsbc=nrow(d),chisq=round(chisq_unif(d$rank,d$L[1]),1))))
  write.csv(res, file.path(CFG$out_dir,"sbc_ranks.csv"), row.names=FALSE)
  write.csv(summ,file.path(CFG$out_dir,"sbc_summary.csv"),row.names=FALSE)
  message(sprintf("SBC: %d draws aggregated. Ranks -> sbc_ranks.csv (histogram per param should be flat/Uniform).",
                  length(unique(res$sbc))))
}
