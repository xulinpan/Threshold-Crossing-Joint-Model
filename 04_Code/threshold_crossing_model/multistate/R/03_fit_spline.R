## 03_fit_spline.R ---------------------------------------------------------
## Fit the SPLINE (hinge) multi-state threshold-crossing model to the real CML
## cohort and compare relapse calibration against the convex-quadratic model.
## The linear-spline term (patient-specific slope change zeta_i after knot
## kappa) lets individual trajectories rebound, resolving relapse
## under-prediction. Everything stays closed-form (piecewise-quadratic).
## -------------------------------------------------------------------------
this_dir <- local({ a <- commandArgs(FALSE); m <- grep("^--file=", a, value = TRUE)
  if (length(m)) dirname(normalizePath(sub("^--file=", "", m[1]))) else getwd() })
suppressMessages({ library(cmdstanr); library(posterior) })
cD <- -4.5; cF <- -5.0; KAPPA <- log1p(2.0)     # knot at 2 years

## locate the real longitudinal file
find_data <- function() {
  cand <- c(file.path(this_dir, "../../../../03_Data/Processed/real_longitudinal_analysis.csv"),
            file.path(this_dir, "../../../03_Data/Processed/real_longitudinal_analysis.csv"))
  for (p in cand) if (file.exists(p)) return(normalizePath(p))
  stop("real_longitudinal_analysis.csv not found; set the path manually.")
}
d <- read.csv(find_data()); d$t <- d$t_months / 12; d <- d[order(d$patient_num, d$t), ]

## build per-patient inputs (confirmed relapse = two consecutive visits above cD)
ids <- sort(unique(d$patient_num)); N <- length(ids)
pid <- ell <- y <- fl <- bm <- integer(0); ll <- list()
onset <- relapse <- integer(N); onL <- onR <- ellC <- rlL <- rlR <- numeric(N)
for (ii in seq_along(ids)) {
  g <- d[d$patient_num == ids[ii], ]; t <- g$t; yy <- g$log_mrd; b <- g$sample_bm
  edges <- c(0, t); Ci <- max(t); dv <- which(yy <= cD); on <- as.integer(length(dv) > 0)
  oL <- Ci; oR <- Ci; rel <- 0L; rL <- Ci; rR <- Ci
  if (on == 1) {
    k <- dv[1]; oL <- edges[k]; oR <- t[k]
    if (k < length(t)) for (r in (k + 1):(length(t) - 1)) if (r + 1 <= length(t)) if (yy[r] > cD && yy[r + 1] > cD) { rel <- 1L; rL <- t[r - 1]; rR <- t[r]; break }
  }
  onset[ii] <- on; relapse[ii] <- rel; ellC[ii] <- log1p(Ci)
  onL[ii] <- log1p(oL); onR[ii] <- log1p(oR); rlL[ii] <- log1p(rL); rlR[ii] <- log1p(rR)
  for (j in seq_along(t)) { pid <- c(pid, ii); ell <- c(ell, log1p(t[j]))
    y <- c(y, ifelse(yy[j] <= cF, cF, yy[j])); fl <- c(fl, as.integer(yy[j] <= cF)); bm <- c(bm, b[j]) }
}
sdata <- list(N = N, Nobs = length(pid), pid = pid, ell = ell, y = y, floor_ind = fl, bm = bm,
              cF = cF, cD = cD, W = 1.0, kappa = KAPPA,
              onset = onset, onL = onL, onR = onR, ellC = ellC,
              relapse = relapse, rlL = rlL, rlR = rlR)
cat(sprintf("N=%d Nobs=%d onset=%d confirmed_relapse=%d\n", N, length(pid), sum(onset), sum(relapse)))

init_sp <- function() list(beta0 = -1.8, beta1 = -3.6, beta2 = 0.5, beta_bm = 0.5,
                           sigma_y = 1.6, tau0 = 1.0, tau1 = 1.8, sigma_zeta = 1.0,
                           sigma_thr = 0.3, z0 = rep(0, N), z1 = rep(0, N), zz = rep(0, N))
mod <- cmdstan_model(file.path(this_dir, "..", "stan", "multistate_spline.stan"))
## Longer run + higher adapt_delta: sigma_thr is pushed toward 0 (the flexible
## spline explains crossings), which stiffens the geometry and lowers ESS.
fit <- mod$sample(data = sdata, chains = 4, parallel_chains = 4,
                  iter_warmup = 2000, iter_sampling = 2000, adapt_delta = 0.999,
                  max_treedepth = 14, init = init_sp, seed = 20260709, refresh = 200)

dr <- as_draws_df(fit$draws()); OUT <- file.path(this_dir, "..", "outputs")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
saveRDS(dr, file.path(OUT, "spline_real_draws.rds"))   # full draws for trace plots
keys <- c("beta0","beta1","beta2","beta_bm","sigma_y","tau0","tau1","sigma_zeta","sigma_thr")
sm <- subset(summarise_draws(dr, "mean","sd", ~quantile2(.x, probs = c(0.05,0.95)), "rhat","ess_bulk"),
             variable %in% keys)
write.csv(sm, file.path(OUT, "spline_real_posterior.csv"), row.names = FALSE)
## relapse calibration
drdf <- as.data.frame(dr); pr <- colMeans(drdf[grep("^p_relapse_byC\\[", names(drdf))])
obs <- ifelse(onset == 1 & relapse == 1, 1, 0)
cat(sprintf("relapse: pred=%.3f obs=%.3f Brier=%.3f\n", mean(pr), mean(obs), mean((pr - obs)^2)))
cat("divergences:", sum(fit$diagnostic_summary()$num_divergent), "\n")
print(sm[, c("variable","mean","q5","q95","rhat")], digits = 3)
message("Spline fit complete -> ", OUT)
