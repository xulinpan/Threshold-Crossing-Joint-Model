## 01_simulate_ms.R --------------------------------------------------------
## Generative model + Stan-data builder for the bidirectional multi-state
## latent threshold-crossing model. Validated against a Python marginal-ML
## refit (parameters, including sigma_thr, recovered accurately).
## -------------------------------------------------------------------------
CF <- -5.0; CD <- -4.5; WIN <- 1.0     # floor, DMR threshold, durability window (yr)
TRUTH_MS <- list(beta0 = -2.0, beta1 = -3.6, beta2 = 1.30, beta_bm = 0.30,
                 sigma_y = 0.80, tau0 = 0.55, tau1 = 0.75, sigma_thr = 0.15)
ellf <- function(t) log1p(t)
mquad <- function(a, g, b2, l) a + g * l + b2 * l^2

simulate_ms <- function(n, seed = 1, p = TRUTH_MS, p_bm = 0.4) {
  set.seed(seed)
  long <- list(); pat <- vector("list", n)
  for (i in seq_len(n)) {
    u0 <- rnorm(1, 0, p$tau0); u1 <- rnorm(1, 0, p$tau1)
    a <- p$beta0 + u0; g <- p$beta1 + u1; b2 <- p$beta2
    Ci <- runif(1, 3, 10)
    Cthr <- rnorm(1, CD, p$sigma_thr)                 # per-patient effective threshold
    ## visits
    t <- numeric(0); cur <- 0
    repeat { cur <- cur + rlnorm(1, log(0.55), 0.5); if (cur > Ci) break; t <- c(t, cur) }
    if (length(t) == 0 || min(t) > 0.3) t <- c(runif(1, .03, .12), t)
    t <- sort(t)
    bm <- rbinom(length(t), 1, p_bm)
    y  <- mquad(a, g, b2, ellf(t)) + p$beta_bm * bm + rnorm(length(t), 0, p$sigma_y)
    fl <- as.integer(y <= CF); yobs <- ifelse(fl == 1, CF, y)
    ## true onset / relapse from quadratic roots vs Cthr
    disc <- g^2 - 4 * b2 * (a - Cthr)
    t_on <- Inf; t_off <- Inf
    if (disc > 0) {
      lon <- (-g - sqrt(disc)) / (2 * b2); loff <- (-g + sqrt(disc)) / (2 * b2)
      if (lon >= 0) { t_on <- expm1(lon); t_off <- expm1(loff) }
      else if (mquad(a, g, b2, 0) <= Cthr) { t_on <- 0; t_off <- expm1(loff) }
    }
    onset <- as.integer(t_on <= Ci)
    relapse <- as.integer(onset == 1 && t_off <= Ci)
    ## bracket a time by surrounding visits -> (ell_L, ell_R)
    edges <- c(0, t)
    br <- function(tt) { k <- findInterval(tt, edges); L <- edges[max(k, 1)]
                         R <- if (k < length(edges)) edges[k + 1] else Ci; c(ellf(L), ellf(R)) }
    onb <- if (onset == 1) br(t_on) else c(ellf(Ci), ellf(Ci))
    rlb <- if (relapse == 1) br(t_off) else c(ellf(Ci), ellf(Ci))
    for (j in seq_along(t)) long[[length(long) + 1]] <-
      data.frame(pid = i, t = t[j], ell = ellf(t[j]), y = yobs[j], floor = fl[j], bm = bm[j])
    pat[[i]] <- data.frame(pid = i, Cend = Ci, ellC = ellf(Ci),
                           onset = onset, onL = onb[1], onR = onb[2],
                           relapse = relapse, rlL = rlb[1], rlR = rlb[2],
                           sustained = as.integer(onset == 1 && (t_off - t_on >= WIN)),
                           t_on = t_on, t_off = t_off)
  }
  list(long = do.call(rbind, long), pat = do.call(rbind, pat), n = n)
}

build_stan_ms <- function(coh) {
  P <- coh$pat
  list(N = coh$n, Nobs = nrow(coh$long),
       pid = coh$long$pid, ell = coh$long$ell, y = coh$long$y,
       floor_ind = coh$long$floor, bm = coh$long$bm,
       cF = CF, cD = CD, W = WIN,
       onset = P$onset, onL = P$onL, onR = P$onR, ellC = P$ellC,
       relapse = P$relapse, rlL = P$rlL, rlR = P$rlR)
}

if (sys.nframe() == 0) {
  d <- simulate_ms(400, seed = 3); P <- d$pat
  cat(sprintf("n=%d onset=%.3f relapse|onset=%.3f sustained|onset=%.3f floor=%.3f visits/pt=%.1f\n",
      d$n, mean(P$onset), mean(P$relapse[P$onset == 1]),
      mean(P$sustained[P$onset == 1]), mean(d$long$floor), nrow(d$long) / d$n))
}
