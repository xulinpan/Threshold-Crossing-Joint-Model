## =============================================================================
## 08_simulation_multistate.R  —  companion to 06 (Reviewer point 4.4)
## ADEMP simulation for the MULTI-STATE spline threshold-crossing model, fit
## with FULL HMC (multistate_spline.stan). Adds the axes that pertain to this
## model and that 06 cannot cover: THRESHOLD at cD in {-4.5,-5.0} and RELAPSE
## FREQUENCY {low,moderate,high} (controlled by the rebound SD sigma_zeta).
## Same sharding / checkpoint / time-budget design as 06.
##
## RUN (PowerShell, from this R/ folder):
##   $env:QUICK_TEST=1; Rscript 08_simulation_multistate.R    # smoke test first
##   ./run_shards_ms.ps1                                      # real run (sharded)
##   Rscript 09_aggregate_ms.R                                # roll up anytime
## REQUIREMENTS: R>=4.1, cmdstanr+CmdStan, posterior, dplyr.
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior); library(dplyr) })

N_SHARDS <- as.integer(Sys.getenv("N_SHARDS","1")); SHARD_ID <- as.integer(Sys.getenv("SHARD_ID","0"))
CFG <- list(
  TIME_BUDGET_HOURS = 10,
  n_rep_central = 200L, n_rep_n300 = 100L, n_rep_scenario = 100L,
  sample_sizes = c(87L,150L,300L),
  chains = 4L, parallel_chains = 4L, iter_warmup = 800L, iter_sampling = 800L,
  adapt_delta = 0.95, max_treedepth = 12L,
  KAPPA = log1p(2), W = log1p(1),                # spline knot (2y), durability window (1y)
  seed = 20260712L,
  out_dir  = file.path("..","outputs","sim_ms"),
  ckpt_dir = file.path("..","outputs","sim_ms","checkpoints"))
dir.create(CFG$ckpt_dir, recursive=TRUE, showWarnings=FALSE); set.seed(CFG$seed)
if (Sys.getenv("QUICK_TEST","0")=="1"){
  CFG$n_rep_central<-2L; CFG$n_rep_n300<-1L; CFG$n_rep_scenario<-1L
  message("QUICK_TEST mode: reps reduced for a smoke test.") }
STAN_MS <- normalizePath(file.path("..","multistate","stan","multistate_spline.stan"))
message(sprintf("Shard %d of %d | %d cores | %d chains/fit.",
                SHARD_ID,N_SHARDS,parallel::detectCores(),CFG$parallel_chains))

## ---- Truth + generative trajectory -----------------------------------------
## Convex spline: m(ell)=beta0 + beta1*ell + beta2*ell^2 + zeta*(ell-kappa)_+  (+ RE)
TRUTH <- list(beta0=-2.0, beta1=-3.6, beta2=1.0, beta_bm=0.6,
              sigma_y=0.8, tau0=0.55, tau1=0.75, sigma_thr=0.15)
## relapse frequency is set by sigma_zeta (larger rebound -> more relapse)
SIGMA_ZETA <- c(low=3.0, moderate=6.0, high=10.0)
mfull <- function(a,g,beta2,ze,kappa,ell) a + g*ell + beta2*ell^2 + ze*pmax(ell-kappa,0)
rre <- function(fam,n,sd) if (fam=="gaussian") rnorm(n,0,sd) else sd*rt(n,3)/sqrt(3)
floor_draw <- function(kind,n,cF) if (kind=="fixed") rep(cF,n) else cF+rnorm(n,0,0.3)
simulate_visits <- function(n,max_t=5,rate=4) lapply(seq_len(n),function(i){
  t<-0; tm<-numeric(0); while(t<max_t){t<-t+rexp(1,rate); if(t<max_t) tm<-c(tm,t)}
  sort(unique(round(c(0.25,tm),3))) })

simulate_dataset <- function(n, scen){
  cF <- -5.0; cD <- scen$cD; kap <- CFG$KAPPA
  sz <- SIGMA_ZETA[[scen$relapse]]
  b0<-rre(scen$re_family,n,TRUTH$tau0); b1<-rre(scen$re_family,n,TRUTH$tau1)
  zeta<-rre("gaussian",n,sz)
  if (scen$trajectory=="quadratic") zeta<-rep(0,n)   # misspec: no rebound in truth
  visits<-simulate_visits(n); floors<-floor_draw(scen$floor_kind,n,cF)
  long<-list()
  onset<-integer(n); onL<-onR<-ellC<-numeric(n); relapse<-integer(n); rlL<-rlR<-numeric(n)
  for (i in seq_len(n)){
    vt<-visits[[i]]; if (length(vt)<3) vt<-c(0.25,1,2.5); ell<-log1p(vt)
    a<-TRUTH$beta0+b0[i]; g<-TRUTH$beta1+b1[i]; ze<-zeta[i]
    m<-mfull(a,g,TRUTH$beta2,ze,kap,ell)                 # latent trajectory at visits
    Ci<-cD+rnorm(1,0,TRUTH$sigma_thr)                    # per-patient effective threshold
    bm_i<-rbinom(length(vt),1,0.9)                       # 90% bone marrow (identifies beta_bm)
    mu<-m+TRUTH$beta_bm*bm_i
    y<-mu+rre(scen$err_family,length(mu),TRUTH$sigma_y)
    fl<-as.integer(y<=floors[i]); y[fl==1]<-floors[i]
    long[[i]]<-data.frame(pid=i, ell=ell, y=y, floor_ind=fl, bm=bm_i)
    ellC[i]<-tail(ell,1)
    k<-which(m<=Ci)[1]                                   # first documented DMR (onset)
    if (is.na(k)){ onset[i]<-0L; onL[i]<-onR[i]<-0; relapse[i]<-0L; rlL[i]<-rlR[i]<-0
    } else {
      onset[i]<-1L; onR[i]<-ell[k]; onL[i]<-if (k>1) ell[k-1] else 0
      j<-k+which(m[(k+1):length(m)]>Ci)[1]               # first upward crossing after onset
      if (!is.na(j)){ relapse[i]<-1L; rlR[i]<-ell[j]; rlL[i]<-ell[j-1]
      } else { relapse[i]<-0L; rlL[i]<-rlR[i]<-0 }
    }
  }
  list(long=do.call(rbind,long), N=n, cF=cF, cD=cD,
       onset=onset,onL=onL,onR=onR,ellC=ellC,relapse=relapse,rlL=rlL,rlR=rlR,
       sigma_zeta_true=sz) }

make_stan_data <- function(ds){ L<-ds$long; list(
  N=ds$N, Nobs=nrow(L), pid=as.integer(L$pid), ell=L$ell, y=L$y,
  floor_ind=as.integer(L$floor_ind), bm=as.integer(L$bm),
  cF=ds$cF, cD=ds$cD, W=CFG$W, kappa=CFG$KAPPA,
  onset=ds$onset, onL=ds$onL, onR=ds$onR, ellC=ds$ellC,
  relapse=ds$relapse, rlL=ds$rlL, rlR=ds$rlR) }

## ---- HMC + metrics ---------------------------------------------------------
fit_hmc <- function(mod,sd){ fit<-mod$sample(data=sd,chains=CFG$chains,
  parallel_chains=CFG$parallel_chains,iter_warmup=CFG$iter_warmup,iter_sampling=CFG$iter_sampling,
  adapt_delta=CFG$adapt_delta,max_treedepth=CFG$max_treedepth,refresh=0,
  show_messages=FALSE,show_exceptions=FALSE)
  dr<-fit$draws(format="draws_df"); dg<-fit$diagnostic_summary(quiet=TRUE)
  list(draws=dr,ndiv=sum(dg$num_divergent),ebfmi=min(dg$ebfmi),
       rhat_max=max(summarise_draws(dr,"rhat")$rhat,na.rm=TRUE)) }
PARAMS <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_zeta","sigma_thr")
metrics_one <- function(draws,param,truth){ x<-as.numeric(draws[[param]]); qi<-quantile(x,c(.025,.975))
  data.frame(param=param,post_mean=mean(x),truth=truth,bias=mean(x)-truth,
             covered=as.integer(truth>=qi[1]&truth<=qi[2]),ci_width=diff(qi)) }

## ---- Queue: central (correctly specified) + reviewer axes ------------------
central <- do.call(rbind, lapply(CFG$sample_sizes, function(nn)
  data.frame(trajectory="spline", cD=-4.5, floor_kind="fixed", err_family="gaussian",
    re_family="gaussian", relapse="moderate", n=nn, tag="central",
    nrep=if (nn==300L) CFG$n_rep_n300 else CFG$n_rep_central, stringsAsFactors=FALSE)))
scen <- transform(rbind(
  data.frame(trajectory="spline",    cD=-5.0, floor_kind="fixed",         err_family="gaussian", re_family="gaussian", relapse="moderate"),
  data.frame(trajectory="spline",    cD=-4.5, floor_kind="fixed",         err_family="gaussian", re_family="gaussian", relapse="low"),
  data.frame(trajectory="spline",    cD=-4.5, floor_kind="fixed",         err_family="gaussian", re_family="gaussian", relapse="high"),
  data.frame(trajectory="quadratic", cD=-4.5, floor_kind="fixed",         err_family="gaussian", re_family="gaussian", relapse="moderate"), # no-rebound truth, spline fit (misspec)
  data.frame(trajectory="spline",    cD=-4.5, floor_kind="heterogeneous", err_family="gaussian", re_family="gaussian", relapse="moderate"),
  data.frame(trajectory="spline",    cD=-4.5, floor_kind="fixed",         err_family="t3",       re_family="t3",       relapse="moderate")),
  n=150L, tag="scenario", nrep=CFG$n_rep_scenario, stringsAsFactors=FALSE)
grid <- rbind(central, scen); grid$scen_id <- seq_len(nrow(grid))
queue <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s){ g<-grid[s,]
  cbind(g[rep(1,g$nrep), setdiff(names(g),"nrep")], rep=seq_len(g$nrep), row.names=NULL) }))
queue$task_id <- sprintf("ms%02d_%s_n%d_r%03d", queue$scen_id, queue$tag, queue$n, queue$rep)
queue <- queue[order(queue$scen_id,queue$rep),]; queue$shard <- (seq_len(nrow(queue))-1L)%%N_SHARDS

done <- sub("\\.rds$","",list.files(CFG$ckpt_dir,pattern="\\.rds$"))
todo <- queue[queue$shard==SHARD_ID & !(queue$task_id %in% done),]
message(sprintf("Queue %d | shard %d has %d | %d remaining.",nrow(queue),SHARD_ID,sum(queue$shard==SHARD_ID),nrow(todo)))
mod <- cmdstan_model(STAN_MS)

## ---- Main loop -------------------------------------------------------------
start<-Sys.time(); nd<-0L
for (r in seq_len(nrow(todo))){
  if (as.numeric(difftime(Sys.time(),start,units="hours"))>CFG$TIME_BUDGET_HOURS){
    message("Time budget reached; stopping. Re-run to resume."); break }
  tk<-as.list(todo[r,]); f<-file.path(CFG$ckpt_dir,paste0(tk$task_id,".rds")); if (file.exists(f)) next
  set.seed(CFG$seed+tk$scen_id*10000L+tk$rep)
  out<-tryCatch({ ds<-simulate_dataset(tk$n,tk); ft<-fit_hmc(mod,make_stan_data(ds))
    truth<-modifyList(TRUTH, list(sigma_zeta=ds$sigma_zeta_true))
    m<-do.call(rbind,lapply(PARAMS,function(p) metrics_one(ft$draws,p,truth[[p]])))
    m$onset_rate<-mean(ds$onset); m$relapse_rate<-mean(ds$relapse[ds$onset==1])
    m$ndiv<-ft$ndiv; m$ebfmi<-ft$ebfmi; m$rhat_max<-ft$rhat_max
    m$comp_fail<-as.integer(ft$ndiv>0|ft$ebfmi<0.2|ft$rhat_max>1.01)
    cbind(task_id=tk$task_id,scen_id=tk$scen_id,tag=tk$tag,n=tk$n,rep=tk$rep,relapse=tk$relapse,cD=tk$cD,m)
  }, error=function(e) data.frame(task_id=tk$task_id,scen_id=tk$scen_id,tag=tk$tag,n=tk$n,rep=tk$rep,param=NA,comp_fail=1L))
  saveRDS(out,f); nd<-nd+1L
  if (nd%%10==0) message(sprintf("  [ms shard %d] ...%d fits (%.1f h)",SHARD_ID,nd,
      as.numeric(difftime(Sys.time(),start,units="hours")))) }
message(sprintf("MS shard %d done this session: %d fits.",SHARD_ID,nd))
if (SHARD_ID==0L) try(source("09_aggregate_ms.R",local=TRUE), silent=TRUE)
