## =============================================================================
## 11_comparators.R  —  comparison with established alternatives
## (Reviewer point 4, required item 10). On the central simulated datasets,
## compares the PROPOSED Bayesian joint model against (a) a standard
## shared-parameter joint model (JMbayes2) that treats documented DMR as an
## exact/right-censored time and the floor as exact, and (b) an interval-
## censored survival model (icenReg) with no longitudinal coupling. Reports
## bias in the biomarker--event association (where defined) and a common
## interval-level Brier score.
##
## NOTE: the JMbayes2 and icenReg calls below are TEMPLATES — adapt them to the
## versions you have installed (both packages' interfaces change between
## releases). The proposed-model fit and the Brier computation are complete.
##
## RUN:  Rscript 11_comparators.R    (set R_REP small first to test)
## REQUIREMENTS: cmdstanr, posterior; optionally JMbayes2, nlme, survival, icenReg
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior) })

CFG <- list(R_REP = 100L, n = 150L, seed = 20260714L,
  chains=4L, parallel_chains=4L, iter_warmup=600L, iter_sampling=500L,
  adapt_delta=0.95, max_treedepth=11L,
  out_dir=file.path("..","outputs","comparators"))
dir.create(CFG$out_dir, recursive=TRUE, showWarnings=FALSE); set.seed(CFG$seed)
if (Sys.getenv("QUICK_TEST","0")=="1") CFG$R_REP<-2L
FLOOR<--5.0
TRUTH<-list(beta0=-2.07,beta1=-3.56,beta2=0.50,beta_bm=0.61,sigma_y=1.81,
            tau0=0.97,tau1=2.20,gamma0=-4.09,gamma1=-0.23,gamma2=-1.83,alpha=-1.27)

## ---- shared data generator (interval-hazard truth; matches proposed model) --
simulate <- function(th,n){
  b0<-rnorm(n,0,th$tau0); b1<-rnorm(n,0,th$tau1); long<-list(); iv<-list()
  for (i in 1:n){
    t<-0;vt<-numeric(0); while(t<5){t<-t+rexp(1,4); if(t<5)vt<-c(vt,t)}; vt<-sort(unique(round(c(0.25,vt),3)))
    ell<-log1p(vt); m<-th$beta0+(th$beta1+b1[i])*ell+th$beta2*ell^2+b0[i]
    bm<-rbinom(length(vt),1,0.9); y<-m+th$beta_bm*bm+rnorm(length(vt),0,th$sigma_y)
    fl<-as.integer(y<=FLOOR); y[fl==1]<-FLOOR
    long[[i]]<-data.frame(pid=i,t=vt,ell=ell,y=y,floor_ind=fl,bm=bm)
    ts<-head(vt,-1);te<-tail(vt,-1);midl<-log1p(.5*(ts+te));gapl<-log1p(te-ts);dl<-pmax(te-ts,1e-6)
    mm<-th$beta0+(th$beta1+b1[i])*midl+th$beta2*midl^2+b0[i]
    logh<-th$gamma0+th$gamma1*midl+th$gamma2*gapl+th$alpha*mm
    ev<-rbinom(length(dl),1,1-exp(-exp(logh)*dl)); first<-which(ev==1)[1]
    keep<-if(is.na(first)) seq_along(ts) else seq_len(first); e<-rep(0L,length(keep)); if(!is.na(first)) e[first]<-1L
    iv[[i]]<-data.frame(pid=i,t_start=ts[keep],t_end=te[keep],event=e,
                        Lyr=ts[keep],Ryr=ifelse(e==1,te[keep],Inf))
  }
  list(long=do.call(rbind,long), iv=do.call(rbind,iv)) }

brier <- function(p,y) mean((p-y)^2)

## ---- (a) proposed model ----------------------------------------------------
mod<-cmdstan_model(normalizePath(file.path("..","stan","interval_hazard_joint.stan")))
fit_proposed <- function(d){
  L<-d$long; I<-d$iv
  sd<-list(N=length(unique(L$pid)),Nobs=nrow(L),pid=L$pid,ell=L$ell,y=L$y,
    floor_ind=L$floor_ind,bm=L$bm,cF=FLOOR,Nint=nrow(I),pid_int=I$pid,
    midlog=log1p(.5*(I$t_start+I$t_end)),gaplog=log1p(I$t_end-I$t_start),
    delta_len=pmax(I$t_end-I$t_start,1e-6),event=I$event)
  f<-mod$sample(data=sd,chains=CFG$chains,parallel_chains=CFG$parallel_chains,
    iter_warmup=CFG$iter_warmup,iter_sampling=CFG$iter_sampling,adapt_delta=CFG$adapt_delta,
    max_treedepth=CFG$max_treedepth,refresh=0,show_messages=FALSE,show_exceptions=FALSE)
  dr<-f$draws(format="draws_df")
  ph<-1-exp(-exp(mean(dr$gamma0)+mean(dr$gamma1)*sd$midlog+mean(dr$gamma2)*sd$gaplog+
              mean(dr$alpha)*(mean(dr$beta0)+mean(dr$beta1)*sd$midlog+mean(dr$beta2)*sd$midlog^2))*sd$delta_len)
  list(alpha_hat=mean(dr$alpha), brier=brier(ph,I$event)) }

## ---- (b) standard joint model (JMbayes2) -- TEMPLATE, adapt to your version -
fit_jmbayes2 <- function(d){
  if (!requireNamespace("JMbayes2",quietly=TRUE)) return(list(alpha_hat=NA,brier=NA))
  ## L<-d$long; pat<-...first documented DMR time as EXACT/right-censored (floor as exact)
  ## lme_fit  <- nlme::lme(y ~ ell + I(ell^2), random = ~ ell | pid, data = L)
  ## cox_fit  <- survival::coxph(Surv(time, status) ~ 1, data = pat, model = TRUE)
  ## jm_fit   <- JMbayes2::jm(cox_fit, lme_fit, time_var = "ell")
  ## extract association coefficient -> alpha_hat; predict -> brier
  list(alpha_hat=NA, brier=NA)                              # TODO: fill for your JMbayes2
}

## ---- (c) interval-censored survival (icenReg) -- TEMPLATE -------------------
fit_icenreg <- function(d){
  if (!requireNamespace("icenReg",quietly=TRUE)) return(list(brier=NA))
  ## pat <- one row per patient with L=Lyr, R=Ryr (Inf if censored)
  ## fit <- icenReg::ic_par(cbind(L, R) ~ 1, data = pat, model = "ph", dist = "weibull")
  ## interval-level predicted probs -> brier   (no longitudinal coupling)
  list(brier=NA)                                            # TODO: fill for your icenReg
}

## ---- run ------------------------------------------------------------------
rows<-list()
for (r in seq_len(CFG$R_REP)){
  set.seed(CFG$seed+r); d<-simulate(TRUTH,CFG$n)
  pr<-tryCatch(fit_proposed(d),error=function(e) list(alpha_hat=NA,brier=NA))
  jm<-tryCatch(fit_jmbayes2(d), error=function(e) list(alpha_hat=NA,brier=NA))
  ic<-tryCatch(fit_icenreg(d),  error=function(e) list(brier=NA))
  rows[[r]]<-data.frame(rep=r,
    proposed_alpha_bias=pr$alpha_hat-TRUTH$alpha, proposed_brier=pr$brier,
    jmbayes2_alpha_bias=jm$alpha_hat-TRUTH$alpha, jmbayes2_brier=jm$brier,
    icenreg_brier=ic$brier)
  if (r%%10==0) message(sprintf("  comparators %d/%d",r,CFG$R_REP))
}
res<-do.call(rbind,rows); write.csv(res,file.path(CFG$out_dir,"comparators_raw.csv"),row.names=FALSE)
summ<-data.frame(
  method=c("Proposed (joint interval+floor)","JMbayes2 (standard joint)","icenReg (interval-only)"),
  alpha_bias=c(mean(res$proposed_alpha_bias,na.rm=TRUE),mean(res$jmbayes2_alpha_bias,na.rm=TRUE),NA),
  brier=c(mean(res$proposed_brier,na.rm=TRUE),mean(res$jmbayes2_brier,na.rm=TRUE),mean(res$icenreg_brier,na.rm=TRUE)))
write.csv(summ,file.path(CFG$out_dir,"comparators_summary.csv"),row.names=FALSE)
print(summ); message("Comparator summary -> comparators_summary.csv")
