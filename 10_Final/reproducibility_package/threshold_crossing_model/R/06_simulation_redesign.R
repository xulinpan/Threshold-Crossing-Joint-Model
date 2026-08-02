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
## =============================================================================

suppressMessages({ library(cmdstanr); library(posterior); library(dplyr) })

## ---- 0. Configuration ------------------------------------------------------
N_SHARDS <- as.integer(Sys.getenv("N_SHARDS", "1"))
SHARD_ID <- as.integer(Sys.getenv("SHARD_ID", "0"))   # 0 .. N_SHARDS-1
CFG <- list(
  TIME_BUDGET_HOURS = 10,
  n_rep_central     = 67L,     # FROZEN: n=87 (200) & n=150 (67) already done & strong
  n_rep_n300        = 30L,     # trimmed further (slow tail)
  n_rep_misspec     = 30L,     # reduced: compute-limited robustness demonstration
  sample_sizes      = c(87L, 150L, 300L),
  chains = 4L, parallel_chains = 4L,     # 4 cores per fit; 4 shards -> 16 cores
  iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,   # 0.95 safe now geometry is healthy (was 0.99)
  DO_SBC = FALSE, DO_COMPARATOR = FALSE, n_sbc = 256L,
  seed = 20260711L,
  out_dir  = file.path("..", "outputs", "sim_redesign"),
  ckpt_dir = file.path("..", "outputs", "sim_redesign", "checkpoints"))
dir.create(CFG$ckpt_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(CFG$seed)
## QUICK_TEST=1 : tiny run (a handful of fits) to verify the harness end-to-end.
if (Sys.getenv("QUICK_TEST","0")=="1"){
  CFG$n_rep_central<-2L; CFG$n_rep_n300<-1L; CFG$n_rep_misspec<-1L
  message("QUICK_TEST mode: reps reduced to a handful for a smoke test.") }
STAN_INTERVAL <- normalizePath(file.path("..", "stan", "interval_hazard_joint.stan"))
message(sprintf("Shard %d of %d | detected %d cores | %d parallel chains/fit.",
                SHARD_ID, N_SHARDS, parallel::detectCores(), CFG$parallel_chains))

## ---- 1. Truth + data-generating mechanisms ---------------------------------
TRUTH <- list(beta0=-2.07, beta1=-3.56, beta2=0.50, beta_bm=0.61,
              sigma_y=1.81, tau0=0.97, tau1=2.20,
              gamma0=-4.09, gamma1=-0.23, gamma2=-1.83, alpha=-1.27)
FLOOR <- -5.0
## latent trajectory m_i(ell); ell = log(1+t). Random slope b1 multiplies ell.
traj_fun <- function(kind, ell, b0, b1, th){
  base <- th$beta0 + b0
  if (kind=="quadratic")        base + (th$beta1+b1)*ell + th$beta2*ell^2
  else if (kind=="monotone")    base + (th$beta1+b1)*ell - 0.5*ell            # misspecified
  else if (kind=="exponential") base + (th$beta1+b1)*(1-exp(-2*ell))          # misspecified
  else stop("unknown trajectory") }
rre <- function(family,n,sd) if (family=="gaussian") rnorm(n,0,sd) else sd*rt(n,3)/sqrt(3)
floor_draw <- function(kind,n) if (kind=="fixed") rep(FLOOR,n) else FLOOR + rnorm(n,0,0.3)
simulate_visits <- function(n, delta1, max_t=5, base_rate=4) lapply(seq_len(n), function(i){
  t<-0; times<-numeric(0)
  while (t<max_t){ t<-t+rexp(1, base_rate); if (t<max_t) times<-c(times,t) }
  if (delta1>0 && length(times)>2) times<-sort(c(times, times[runif(length(times))<delta1*0.3]))
  sort(unique(round(c(0.25,times),3))) })
## Generate one dataset. Longitudinal mean = m_i + beta_bm*bm; events are drawn
## from the SAME interval-hazard likelihood the Stan model fits, so the central
## scenario is correctly specified. bm is a 90/10 BM/PB mix so that beta0 and
## beta_bm are identified (all-BM would confound them and wreck HMC geometry).
simulate_dataset <- function(n, scen){
  b0<-rre(scen$re_family,n,TRUTH$tau0); b1<-rre(scen$re_family,n,TRUTH$tau1)
  visits<-simulate_visits(n, scen$vis_delta1); floors<-floor_draw(scen$floor_kind,n)
  long<-list(); intervals<-list()
  for (i in seq_len(n)){
    vt<-visits[[i]]; if (length(vt)<2) vt<-c(0.25,0.75); ell<-log1p(vt)
    m  <- traj_fun(scen$trajectory, ell, b0[i], b1[i], TRUTH)     # latent trajectory
    bm <- rbinom(length(vt),1,0.9)                                # 90% bone marrow
    mu <- m + TRUTH$beta_bm*bm                                    # observation mean
    y  <- mu + rre(scen$err_family, length(mu), TRUTH$sigma_y)
    fl <- as.integer(y<=floors[i]); y[fl==1]<-floors[i]
    long[[i]] <- data.frame(pid=i, t=vt, ell=ell, y=y, floor_ind=fl, bm=bm)
    ## at-risk intervals between consecutive visits
    ts<-head(vt,-1); te<-tail(vt,-1)
    if (length(ts)>=1){
      midl<-log1p(0.5*(ts+te)); gapl<-log1p(te-ts); dl<-pmax(te-ts,1e-6)
      m_mid<-traj_fun(scen$trajectory, midl, b0[i], b1[i], TRUTH) # latent at midpoint
      logh<-TRUTH$gamma0+TRUTH$gamma1*midl+TRUTH$gamma2*gapl+TRUTH$alpha*m_mid
      p<-1-exp(-exp(logh)*dl)
      ev<-rbinom(length(p),1,p); first<-which(ev==1)[1]
      keep<-if (is.na(first)) seq_along(ts) else seq_len(first)   # exit at first DMR
      evc<-rep(0L,length(keep)); if (!is.na(first)) evc[first]<-1L
      intervals[[i]]<-data.frame(pid=i, t_start=ts[keep], t_end=te[keep], event=evc)
    }
  }
  list(long=do.call(rbind,long), intervals=do.call(rbind,intervals)) }
make_stan_data <- function(ds){ L<-ds$long; I<-ds$intervals; list(
  N=length(unique(L$pid)), Nobs=nrow(L), pid=L$pid, ell=L$ell, y=L$y,
  floor_ind=L$floor_ind, bm=L$bm, cF=FLOOR,
  Nint=nrow(I), pid_int=I$pid, midlog=log1p(0.5*(I$t_start+I$t_end)),
  gaplog=log1p(I$t_end-I$t_start), delta_len=pmax(I$t_end-I$t_start,1e-6), event=I$event) }

## ---- 2. HMC fit + diagnostics ----------------------------------------------
fit_hmc <- function(mod, standata){
  fit<-mod$sample(data=standata, chains=CFG$chains, parallel_chains=CFG$parallel_chains,
                  iter_warmup=CFG$iter_warmup, iter_sampling=CFG$iter_sampling,
                  adapt_delta=CFG$adapt_delta, max_treedepth=CFG$max_treedepth,
                  refresh=0, show_messages=FALSE, show_exceptions=FALSE)
  dr<-fit$draws(format="draws_df"); dg<-fit$diagnostic_summary(quiet=TRUE)
  list(draws=dr, ndiv=sum(dg$num_divergent), ebfmi=min(dg$ebfmi),
       rhat_max=max(summarise_draws(dr,"rhat")$rhat, na.rm=TRUE)) }

## ---- 3. Performance measures -----------------------------------------------
PARAMS <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","gamma0","gamma1","gamma2","alpha")
metrics_one <- function(draws, param, truth){
  x<-as.numeric(draws[[param]]); qi<-quantile(x,c(0.025,0.975))
  data.frame(param=param, post_mean=mean(x), truth=truth, bias=mean(x)-truth,
             covered=as.integer(truth>=qi[1] & truth<=qi[2]), ci_width=diff(qi)) }

## ---- 4. Build FULL queue with FIXED shard assignment -----------------------
## Execution order (scen_id): n=87 and n=150 FIRST (keep their finished
## checkpoints: n=87=scen1, n=150=scen2), then misspecification cells, then the
## slow n=300 LAST. Reordering does not disturb scen1/scen2 task_ids.
central_fast <- do.call(rbind, lapply(c(87L,150L), function(nn)
  data.frame(trajectory="quadratic", floor_kind="fixed",
             err_family="gaussian", re_family="gaussian", vis_delta1=0, n=nn,
             tag="central", nrep=CFG$n_rep_central, stringsAsFactors=FALSE)))
## Misspecification scenarios (one departure at a time), n=150.
## NB: thresholds -4.5/-5.0 and relapse frequency pertain to the MULTI-STATE
## model and are evaluated in the companion script (08), not here.
misspec <- transform(rbind(
  data.frame(trajectory="monotone",    floor_kind="fixed",         err_family="gaussian", re_family="gaussian", vis_delta1=0),
  data.frame(trajectory="exponential", floor_kind="fixed",         err_family="gaussian", re_family="gaussian", vis_delta1=0),
  data.frame(trajectory="quadratic",   floor_kind="heterogeneous", err_family="gaussian", re_family="gaussian", vis_delta1=0),
  data.frame(trajectory="quadratic",   floor_kind="fixed",         err_family="t3",       re_family="t3",       vis_delta1=0),
  data.frame(trajectory="quadratic",   floor_kind="fixed",         err_family="gaussian", re_family="gaussian", vis_delta1=1.0)),
  n=150L, tag="misspec", nrep=CFG$n_rep_misspec, stringsAsFactors=FALSE)
central_n300 <- data.frame(trajectory="quadratic", floor_kind="fixed",
  err_family="gaussian", re_family="gaussian", vis_delta1=0, n=300L,
  tag="central", nrep=CFG$n_rep_n300, stringsAsFactors=FALSE)
grid <- rbind(central_fast, misspec, central_n300); grid$scen_id <- seq_len(nrow(grid))
queue <- do.call(rbind, lapply(seq_len(nrow(grid)), function(s){
  g<-grid[s,]; cbind(g[rep(1,g$nrep), setdiff(names(g),"nrep")], rep=seq_len(g$nrep), row.names=NULL)}))
queue$task_id <- sprintf("s%02d_%s_n%d_r%03d", queue$scen_id, queue$tag, queue$n, queue$rep)
queue <- queue[order(queue$scen_id, queue$rep), ]                 # deterministic order
queue$shard <- (seq_len(nrow(queue)) - 1L) %% N_SHARDS            # FIXED assignment

## ---- 5. Select this shard's outstanding tasks ------------------------------
done_ids <- sub("\\.rds$","", list.files(CFG$ckpt_dir, pattern="\\.rds$"))
todo <- queue[queue$shard==SHARD_ID & !(queue$task_id %in% done_ids), ]
message(sprintf("Queue total %d | this shard %d tasks | %d remaining after resume.",
                nrow(queue), sum(queue$shard==SHARD_ID), nrow(todo)))
mod <- cmdstan_model(STAN_INTERVAL)

## ---- 6. Main loop: checkpointed + time-budgeted ----------------------------
start<-Sys.time(); n_done_now<-0L
for (r in seq_len(nrow(todo))){
  if (as.numeric(difftime(Sys.time(),start,units="hours")) > CFG$TIME_BUDGET_HOURS){
    message("Time budget reached; stopping cleanly. Re-run to resume."); break }
  task<-as.list(todo[r,]); f<-file.path(CFG$ckpt_dir, paste0(task$task_id,".rds"))
  if (file.exists(f)) next
  set.seed(CFG$seed + task$scen_id*10000L + task$rep)
  out<-tryCatch({
    ds<-simulate_dataset(task$n, task); ft<-fit_hmc(mod, make_stan_data(ds))
    m<-do.call(rbind, lapply(PARAMS, function(p) metrics_one(ft$draws, p, TRUTH[[p]])))
    m$ndiv<-ft$ndiv; m$ebfmi<-ft$ebfmi; m$rhat_max<-ft$rhat_max
    m$comp_fail<-as.integer(ft$ndiv>0 | ft$ebfmi<0.2 | ft$rhat_max>1.01)
    cbind(task_id=task$task_id, scen_id=task$scen_id, tag=task$tag, n=task$n, rep=task$rep, m)
  }, error=function(e) data.frame(task_id=task$task_id, scen_id=task$scen_id,
       tag=task$tag, n=task$n, rep=task$rep, param=NA, comp_fail=1L))
  saveRDS(out, f); n_done_now<-n_done_now+1L
  if (n_done_now %% 10 == 0) message(sprintf("  [shard %d] ...%d fits this session (%.1f h)",
      SHARD_ID, n_done_now, as.numeric(difftime(Sys.time(),start,units="hours")))) }
message(sprintf("Shard %d finished this session: %d new fits.", SHARD_ID, n_done_now))

## ---- 7. Aggregate (only shard 0 writes, to avoid concurrent-write clashes) --
if (SHARD_ID == 0L) source("07_aggregate_sim.R", local=TRUE)
