## =============================================================================
## 11_comparators.R  —  comparison with established alternatives
## (Reviewer point 4, required item 10). On the central simulated datasets,
## compares the PROPOSED Bayesian joint model against (a) a standard
## shared-parameter joint model (JMbayes2) that treats documented DMR as an
## exact/right-censored time and the floor as exact, and (b) an interval-
## censored survival model (icenReg) with no longitudinal coupling.
##
## STATUS. The JMbayes2 and icenReg arms are still TEMPLATES (both packages
## change interfaces between releases). The script now REFUSES to write a
## summary in which a comparator arm is entirely NA, instead of silently
## producing a two-row table of NAs that looks like a result.
##
## CHANGED 2026-07-29
##  (1) The private copy of the data-generating mechanism (a fourth copy) is
##      replaced by dgm_common.R.
##  (2) The Brier score was wrong. It plugged POSTERIOR MEANS into the hazard and
##      omitted b0/b1 entirely, so it scored a population-average predicted
##      probability against subject-level events. It is now the posterior mean of
##      the predicted probability, computed per draw and INCLUDING the
##      subject-specific random effects.
##      Still an APPARENT (in-sample) score: every method is scored on the data
##      it was fitted to. That is a common yardstick but an optimistic one; a
##      within-patient temporal split would be needed for an honest comparison.
##
## RUN:  Rscript 11_comparators.R    (set R_REP small first to test)
## REQUIREMENTS: cmdstanr, posterior; optionally JMbayes2, nlme, survival, icenReg
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior) })

.cmp_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
source(file.path(.cmp_dir, "dgm_common.R"))
dgm_describe()

CFG <- list(R_REP = 100L, n = 150L, seed = 20260714L,
  chains = 4L, parallel_chains = 4L, iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,
  out_dir = file.path(.cmp_dir, "..", "outputs", "comparators"))
dir.create(CFG$out_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(CFG$seed)
if (Sys.getenv("QUICK_TEST", "0") == "1") CFG$R_REP <- 2L

CENTRAL <- list(trajectory = "quadratic", floor_kind = "fixed",
                err_family = "gaussian", re_family = "gaussian", informative = FALSE)

brier <- function(p, y) mean((p - y)^2)

## ---- (a) proposed model ----------------------------------------------------
mod <- cmdstan_model(normalizePath(file.path(.cmp_dir, "..", "stan",
                                             "interval_hazard_joint.stan")))
fit_proposed <- function(d) {
  sdat <- make_stan_data(d)
  f <- mod$sample(data = sdat, chains = CFG$chains, parallel_chains = CFG$parallel_chains,
    iter_warmup = CFG$iter_warmup, iter_sampling = CFG$iter_sampling,
    adapt_delta = CFG$adapt_delta, max_treedepth = CFG$max_treedepth,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  dr <- f$draws(format = "draws_df")

  ## posterior mean of the interval-level event probability, per draw, with the
  ## subject-specific random effects included (b0/b1 are transformed parameters)
  mid <- sdat$midlog; gap <- sdat$gaplog; dl <- sdat$delta_len; pid <- sdat$pid_int
  B0 <- as.matrix(dr[, sprintf("b0[%d]", seq_len(sdat$N)), drop = FALSE])
  B1 <- as.matrix(dr[, sprintf("b1[%d]", seq_len(sdat$N)), drop = FALSE])
  acc <- numeric(length(mid))
  D <- nrow(dr)
  for (j in seq_len(D)) {
    m    <- dr$beta0[j] + dr$beta1[j]*mid + dr$beta2[j]*mid^2 +
            B0[j, pid] + B1[j, pid]*mid
    logh <- dr$gamma0[j] + dr$gamma1[j]*mid + dr$gamma2[j]*gap + dr$alpha[j]*m
    acc  <- acc + (1 - exp(-exp(logh) * dl))
  }
  ph <- acc / D
  list(alpha_hat = mean(dr$alpha), brier = brier(ph, sdat$event))
}

## ---- (b) standard joint model (JMbayes2) -- TEMPLATE, adapt to your version -
fit_jmbayes2 <- function(d) {
  if (!requireNamespace("JMbayes2", quietly = TRUE)) return(list(alpha_hat=NA, brier=NA))
  ## L<-d$long; pat<-...first documented DMR time as EXACT/right-censored (floor as exact)
  ## lme_fit  <- nlme::lme(y ~ ell + I(ell^2), random = ~ ell | pid, data = L)
  ## cox_fit  <- survival::coxph(Surv(time, status) ~ 1, data = pat, model = TRUE)
  ## jm_fit   <- JMbayes2::jm(cox_fit, lme_fit, time_var = "ell")
  ## extract association coefficient -> alpha_hat; predict -> brier
  list(alpha_hat = NA, brier = NA)                          # TODO: fill for your JMbayes2
}

## ---- (c) interval-censored survival (icenReg) -- TEMPLATE -------------------
fit_icenreg <- function(d) {
  if (!requireNamespace("icenReg", quietly = TRUE)) return(list(brier = NA))
  ## pat <- one row per patient with L = last event-free t_end, R = t_end at DMR (Inf if censored)
  ## fit <- icenReg::ic_par(cbind(L, R) ~ 1, data = pat, model = "ph", dist = "weibull")
  ## interval-level predicted probs -> brier   (no longitudinal coupling)
  list(brier = NA)                                          # TODO: fill for your icenReg
}

## ---- run ------------------------------------------------------------------
rows <- list()
for (r in seq_len(CFG$R_REP)) {
  set.seed(CFG$seed + r)
  d  <- simulate_dataset(CFG$n, CENTRAL)
  pr <- tryCatch(fit_proposed(d), error = function(e) list(alpha_hat=NA, brier=NA))
  jm <- tryCatch(fit_jmbayes2(d), error = function(e) list(alpha_hat=NA, brier=NA))
  ic <- tryCatch(fit_icenreg(d),  error = function(e) list(brier=NA))
  rows[[r]] <- data.frame(rep = r,
    proposed_alpha_bias = pr$alpha_hat - TRUTH$alpha, proposed_brier = pr$brier,
    jmbayes2_alpha_bias = jm$alpha_hat - TRUTH$alpha, jmbayes2_brier = jm$brier,
    icenreg_brier = ic$brier)
  if (r %% 10 == 0) message(sprintf("  comparators %d/%d", r, CFG$R_REP))
}
res <- do.call(rbind, rows)
write.csv(res, file.path(CFG$out_dir, "comparators_raw.csv"), row.names = FALSE)

mse <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
summ <- data.frame(
  method = c("Proposed (joint interval+floor)", "JMbayes2 (standard joint)",
             "icenReg (interval-only)"),
  nrep = c(sum(!is.na(res$proposed_brier)), sum(!is.na(res$jmbayes2_brier)),
           sum(!is.na(res$icenreg_brier))),
  alpha_bias = c(mse(res$proposed_alpha_bias), mse(res$jmbayes2_alpha_bias), NA),
  brier = c(mse(res$proposed_brier), mse(res$jmbayes2_brier), mse(res$icenreg_brier)))
write.csv(summ, file.path(CFG$out_dir, "comparators_summary.csv"), row.names = FALSE)
print(summ, row.names = FALSE)

empty <- summ$method[summ$nrep == 0]
if (length(empty))
  warning("These comparator arms produced no results and MUST NOT be reported as a ",
          "comparison until implemented: ", paste(empty, collapse = "; "),
          call. = FALSE, immediate. = TRUE)
message("Comparator summary -> comparators_summary.csv")
