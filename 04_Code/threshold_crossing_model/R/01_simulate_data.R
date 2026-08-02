## 01_simulate_data.R ------------------------------------------------------
## Generative model (matches bayesian_hierarchical_threshold_crossing_model.tex)
## and Stan-data builders for both the threshold-crossing and interval-hazard
## models. Sourced by 02/03/04.
## -------------------------------------------------------------------------

ell <- function(t) log1p(t)   # l(t) = log(1 + t)

## latent trajectory value at time t (years) for given random effects
m_latent <- function(t, u0, u1, p = TRUTH) {
  L <- ell(t)
  p$beta0 + p$beta1 * L + p$beta2 * L^2 + u0 + u1 * L
}

## ---- simulate one cohort -------------------------------------------------
## informative = TRUE  -> visit intensity depends on latent MRD (NHPP thinning)
## dense        = TRUE  -> shorter target visit gap
simulate_cohort <- function(n, informative = FALSE, dense = FALSE,
                            seed = 1, p = TRUTH,
                            p_bm = 0.90) {
  set.seed(seed)
  gap_med <- if (dense) 0.33 else 0.55
  gap_sd  <- 0.5
  d0 <- -0.7; d1 <- 0.35; d2 <- 0.2; wsd <- 0.5   # informative-visit params

  long <- list(); pat <- vector("list", n)
  for (i in seq_len(n)) {
    u0 <- rnorm(1, 0, p$tau0); u1 <- rnorm(1, 0, p$tau1)
    Ci <- runif(1, 1.5, 7.0)                       # follow-up horizon (years)

    ## ---- visit times ----
    if (!informative) {
      t <- numeric(0); cur <- 0
      repeat {
        cur <- cur + rlnorm(1, log(gap_med), gap_sd)
        if (cur > Ci) break
        t <- c(t, cur)
      }
    } else {
      w <- rnorm(1, 0, wsd)
      grid <- seq(1e-3, Ci, length.out = 50)
      lam_grid <- exp(d0 + d1 * (m_latent(grid, u0, u1, p) + 4.5) + d2 * ell(grid) + w)
      lam_max <- max(lam_grid) * 1.2 + 1e-6
      t <- numeric(0); cur <- 0
      repeat {
        cur <- cur + rexp(1, lam_max)
        if (cur > Ci) break
        lam <- exp(d0 + d1 * (m_latent(cur, u0, u1, p) + 4.5) + d2 * ell(cur) + w)
        if (runif(1) < lam / lam_max) t <- c(t, cur)
      }
    }
    if (length(t) == 0 || min(t) > 0.25) t <- c(runif(1, 0.02, 0.12), t)
    t <- sort(t)

    bm  <- rbinom(length(t), 1, p_bm)
    mu  <- m_latent(t, u0, u1, p) + p$beta_bm * bm
    ys  <- rnorm(length(t), mu, p$sigma_y)
    fl  <- as.integer(ys <= C_F)
    yobs <- ifelse(fl == 1, C_F, ys)

    ## true latent onset (first crossing of c_D)
    fg <- seq(1e-3, Ci, length.out = 400)
    mfg <- m_latent(fg, u0, u1, p)
    onset <- if (any(mfg <= C_D)) fg[which(mfg <= C_D)[1]] else Inf

    ## documented DMR = first visit with observed y <= c_D
    dv <- which(yobs <= C_D)
    delta <- as.integer(length(dv) > 0)
    edges <- c(0, t)
    if (delta == 1) { k <- dv[1]; L <- edges[k]; R <- edges[k + 1] }
    else            { L <- NA;    R <- NA }

    long[[i]] <- data.frame(pid = i, t = t, ell = ell(t), y = yobs,
                            floor = fl, bm = bm)
    pat[[i]]  <- data.frame(pid = i, delta = delta, L = L, R = R,
                            Cend = Ci, true_onset = onset)
  }
  list(long = do.call(rbind, long), pat = do.call(rbind, pat), n = n)
}

## ---- build Stan data: threshold-crossing model ---------------------------
build_stan_threshold <- function(coh, G = GRID, kappa = KAPPA) {
  n <- coh$n; pat <- coh$pat
  ellA <- matrix(0, n, G); ellB <- matrix(0, n, G); ellC <- matrix(0, n, G)
  for (i in seq_len(n)) {
    Ci <- pat$Cend[i]
    if (pat$delta[i] == 1) {
      ellA[i, ] <- ell(seq(0, pat$L[i], length.out = G))
      ellB[i, ] <- ell(seq(0, pat$R[i], length.out = G))
    } else {
      ellA[i, ] <- ell(seq(0, Ci, length.out = G))
      ellB[i, ] <- ellA[i, ]
    }
    ellC[i, ] <- ell(seq(0, Ci, length.out = G))
  }
  list(N = n, Nobs = nrow(coh$long),
       pid = coh$long$pid, ell = coh$long$ell, y = coh$long$y,
       floor_ind = coh$long$floor, bm = coh$long$bm,
       cF = C_F, cD = C_D,
       delta = pat$delta, G = G, ellA = ellA, ellB = ellB, ellC = ellC,
       kappa = kappa)
}

## ---- build at-risk intervals + Stan data: interval-hazard model ----------
build_stan_interval <- function(coh) {
  n <- coh$n; pat <- coh$pat
  rows <- list()
  for (i in seq_len(n)) {
    ti <- coh$long$t[coh$long$pid == i]
    yi <- coh$long$y[coh$long$pid == i]
    edges <- c(0, ti)
    dv <- which(yi <= C_D)
    kmax <- if (length(dv) > 0) dv[1] else length(ti)
    for (k in seq_len(kmax)) {
      L <- edges[k]; R <- edges[k + 1]
      ev <- as.integer(length(dv) > 0 && k == dv[1])
      rows[[length(rows) + 1]] <- data.frame(
        pid = i, midlog = log1p((L + R) / 2), gaplog = log1p(R - L),
        delta_len = R - L, event = ev)
    }
  }
  iv <- do.call(rbind, rows)
  list(N = n, Nobs = nrow(coh$long),
       pid = coh$long$pid, ell = coh$long$ell, y = coh$long$y,
       floor_ind = coh$long$floor, bm = coh$long$bm, cF = C_F,
       Nint = nrow(iv), pid_int = iv$pid, midlog = iv$midlog,
       gaplog = iv$gaplog, delta_len = iv$delta_len, event = iv$event)
}

if (sys.nframe() == 0) {
  coh <- simulate_cohort(120, informative = FALSE, seed = 3)
  cat(sprintf("n=%d  obs=%d  floor_rate=%.3f  documented_DMR=%.3f  median_visits=%.1f\n",
              coh$n, nrow(coh$long), mean(coh$long$floor), mean(coh$pat$delta),
              median(table(coh$long$pid))))
}
