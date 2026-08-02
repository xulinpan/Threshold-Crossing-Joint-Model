## =============================================================================
## 10_sbc.R  —  Simulation-Based Calibration (SBC) for the joint interval model
## (Reviewer point 4, required item 3). Draw theta ~ prior, simulate data, fit,
## record the rank of each true theta among L thinned posterior draws. Ranks are
## Uniform{0..L} iff the posterior is correctly calibrated.
##
## FIXED 2026-07-29
##  (1) The chi-square bins were not equiprobable. Ranks are integers on
##      0..L, i.e. (L+1) atoms; splitting them into B bins gives equal null
##      probabilities only when B divides (L+1). With L = 100 and B = 20 one bin
##      held 6 atoms and nineteen held 5, yet the test used equal expected
##      counts -- which inflates the statistic and makes the test
##      anti-conservative. Expected counts are now exact (proportional to the
##      number of atoms per bin), and L defaults to 99 so that 20 bins of 5
##      atoms divide evenly.
##  (2) Only a bare chi-square value was written, with no df, so it could not be
##      interpreted. df, a p-value and a Holm-adjusted p-value (11 parameters
##      are tested simultaneously) are now reported, together with the minimum
##      expected cell count so low-power runs are visible.
##  (3) The data-generating mechanism was a fourth private copy. It now comes
##      from dgm_common.R, so SBC and the simulation study cannot drift apart.
##
## RUN:  $env:QUICK_TEST=1; Rscript 10_sbc.R   (smoke test)  then  Rscript 10_sbc.R
## =============================================================================
suppressMessages({ library(cmdstanr); library(posterior) })

.sbc_dir <- local({
  a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1])))
  else tryCatch(dirname(normalizePath(sys.frame(1)$ofile)), error = function(e) getwd())
})
source(file.path(.sbc_dir, "dgm_common.R"))

## SBC validates the sampler against the model it fits, so the generative event
## mechanism must be the one in the Stan likelihood.
if (EVENT_MECHANISM != "hazard")
  stop("SBC requires EVENT_MECHANISM=hazard (the mechanism the Stan model encodes).")
dgm_describe()

CFG <- list(
  N_SBC = 150L,          # number of SBC iterations (>=100 recommended)
  n = 87L,               # cohort size per SBC draw (n=87 keeps fits fast)
  L = 99L,               # thinned draws -> ranks on 0..99 = 100 atoms = 20 x 5
  NBIN = 20L,
  chains = 4L, parallel_chains = 4L, iter_warmup = 600L, iter_sampling = 500L,
  adapt_delta = 0.95, max_treedepth = 11L,
  TIME_BUDGET_HOURS = 10, seed = 20260713L,
  out_dir  = file.path(.sbc_dir, "..", "outputs", "sbc"))
CFG$ckpt_dir <- file.path(CFG$out_dir, paste0("checkpoints_", dgm_tag()))
dir.create(CFG$ckpt_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(CFG$seed)
if (Sys.getenv("QUICK_TEST", "0") == "1") { CFG$N_SBC <- 3L; message("QUICK_TEST: 3 SBC draws.") }
STAN <- normalizePath(file.path(.sbc_dir, "..", "stan", "interval_hazard_joint_pp.stan"))

## ---- priors: derived from the SAME object the model is given -----------------
## Previously these were transcribed by hand from the .stan file, so a prior
## edited in one place and not the other would silently invalidate SBC. They are
## now generated from PRIOR_SETS[[PRIOR_SET]] in dgm_common.R, which is also what
## is passed to the model as data -- the two cannot disagree. This also means SBC
## follows PRIOR_SET: validating the wide2 fit requires drawing from wide2.
draw_from_prior <- function(set = PRIOR_SET) {
  p <- PRIOR_SETS[[set]]
  list(
    beta0   = rnorm(1, p$pr_beta0[1],   p$pr_beta0[2]),
    beta1   = rnorm(1, p$pr_beta1[1],   p$pr_beta1[2]),
    beta2   = rnorm(1, p$pr_beta2[1],   p$pr_beta2[2]),
    beta_bm = rnorm(1, p$pr_beta_bm[1], p$pr_beta_bm[2]),
    sigma_y = rexp(1, p$pr_sigma_y),
    tau0    = rexp(1, p$pr_tau0),
    tau1    = rexp(1, p$pr_tau1),
    gamma0  = rnorm(1, p$pr_gamma0[1],  p$pr_gamma0[2]),
    gamma1  = rnorm(1, p$pr_gamma1[1],  p$pr_gamma1[2]),
    gamma2  = rnorm(1, p$pr_gamma2[1],  p$pr_gamma2[2]),
    alpha   = rnorm(1, p$pr_alpha[1],   p$pr_alpha[2]))
}

SBC_SCEN <- list(trajectory = "quadratic", floor_kind = "fixed",
                 err_family = "gaussian", re_family = "gaussian", informative = FALSE)

## ---- rank of the true value among L thinned posterior draws ------------------
posterior_rank <- function(draws, th, L) {
  out <- integer(length(PARAMS)); names(out) <- PARAMS
  for (p in PARAMS) {
    x   <- as.numeric(draws[[p]])
    idx <- round(seq(1, length(x), length.out = L))
    out[p] <- sum(x[idx] < th[[p]])          # rank in 0..L
  }
  out
}

## ---- exact chi-square test of rank uniformity -------------------------------
## Atom a in 0..L is assigned to bin floor(a * B / (L + 1)); bins differ in size
## by at most one atom and the null probability of a bin is (atoms in it)/(L+1).
sbc_chisq <- function(rk, L, B = 20L) {
  rk <- rk[!is.na(rk)]
  if (!length(rk)) return(list(chisq = NA, df = NA, p = NA, min_exp = NA))
  bin_atom <- ((0:L) * B) %/% (L + 1L)                       # NB: %/% binds tighter than *
  w <- as.numeric(table(factor(bin_atom, levels = 0:(B - 1L))))
  p0 <- w / (L + 1L)
  bin_obs <- (rk * B) %/% (L + 1L)
  o <- as.numeric(table(factor(bin_obs, levels = 0:(B - 1L))))
  e <- length(rk) * p0
  stat <- sum((o - e)^2 / e); df <- B - 1L
  list(chisq = stat, df = df, p = pchisq(stat, df, lower.tail = FALSE),
       min_exp = min(e))
}

## ---- main loop --------------------------------------------------------------
mod <- cmdstan_model(STAN); start <- Sys.time()
done <- sub("\\.rds$", "", list.files(CFG$ckpt_dir, pattern = "\\.rds$"))
for (s in seq_len(CFG$N_SBC)) {
  id <- sprintf("sbc_%04d", s); f <- file.path(CFG$ckpt_dir, paste0(id, ".rds"))
  if (id %in% done) next
  if (as.numeric(difftime(Sys.time(), start, units = "hours")) > CFG$TIME_BUDGET_HOURS) {
    message("budget reached"); break }
  set.seed(CFG$seed + s)
  out <- tryCatch({
    th <- draw_from_prior()
    ds <- simulate_dataset(CFG$n, SBC_SCEN, th = th)
    ft <- fit_hmc(mod, c(make_stan_data(ds), stan_priors()), CFG)
    r  <- posterior_rank(ft$draws, th, CFG$L)
    data.frame(sbc = s, param = names(r), rank = as.integer(r), L = CFG$L,
               ndiv = ft$ndiv, ebfmi = ft$ebfmi, rhat_max = ft$rhat_max,
               fit_bad = ft$fit_bad, fit_marginal = ft$fit_marginal,
               dgm_tag = dgm_tag())
  }, error = function(e) data.frame(sbc = s, param = NA, rank = NA, L = CFG$L,
                                    dgm_tag = dgm_tag()))
  saveRDS(out, f)
  if (s %% 10 == 0) message(sprintf("  SBC %d/%d", s, CFG$N_SBC))
}

## ---- aggregate --------------------------------------------------------------
## Pool every checkpoint dir but keep dgm_tag, so mechanisms are never mixed.
ck <- list.dirs(CFG$out_dir, recursive = FALSE)
ck <- ck[grepl("^checkpoints", basename(ck))]
fs <- unlist(lapply(ck, function(d) list.files(d, pattern = "\\.rds$", full.names = TRUE)))
if (length(fs) > 0) {
  res <- do.call(rbind, lapply(fs, function(f) {
    x <- readRDS(f)
    if (!"dgm_tag" %in% names(x)) {
      tg <- sub("^checkpoints_?", "", basename(dirname(f)))
      x$dgm_tag <- if (nzchar(tg)) tg else "v1_legacy_legacy_hazard"
    }
    if (!"fit_bad" %in% names(x)) x$fit_bad <- NA_integer_
    x[, c("sbc","param","rank","L","fit_bad","dgm_tag")]
  }))
  res <- res[!is.na(res$rank), ]
  ## SBC on non-converged fits is not informative about calibration
  nbad <- sum(!is.na(res$fit_bad) & res$fit_bad == 1) / length(PARAMS)
  res  <- res[is.na(res$fit_bad) | res$fit_bad == 0, ]

  summ <- do.call(rbind, lapply(split(res, list(res$dgm_tag, res$param), drop = TRUE),
    function(d) {
      tt <- sbc_chisq(d$rank, d$L[1], CFG$NBIN)
      data.frame(dgm_tag = d$dgm_tag[1], param = d$param[1], nsbc = nrow(d),
                 L = d$L[1], nbin = CFG$NBIN, min_expected = round(tt$min_exp, 2),
                 chisq = round(tt$chisq, 2), df = tt$df, p = signif(tt$p, 4))
    }))
  summ$p_holm <- NA_real_
  for (tg in unique(summ$dgm_tag)) {
    i <- summ$dgm_tag == tg
    summ$p_holm[i] <- signif(p.adjust(summ$p[i], method = "holm"), 4)
  }
  summ <- summ[order(summ$dgm_tag, summ$p), ]
  if (any(summ$min_expected < 5, na.rm = TRUE))
    message("WARNING: some bins have expected count < 5 -- increase N_SBC or reduce NBIN.")
  if (nbad > 0)
    message(sprintf("NOTE: %d non-converged SBC fits excluded.", round(nbad)))
  write.csv(res,  file.path(CFG$out_dir, "sbc_ranks.csv"),   row.names = FALSE)
  write.csv(summ, file.path(CFG$out_dir, "sbc_summary.csv"), row.names = FALSE)
  print(summ, row.names = FALSE)
  message(sprintf("SBC: %d draws aggregated -> sbc_summary.csv (rank histograms should be flat).",
                  length(unique(res$sbc))))
}
