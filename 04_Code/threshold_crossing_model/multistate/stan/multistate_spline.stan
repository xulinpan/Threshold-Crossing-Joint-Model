//==============================================================================
// Multi-state latent threshold-crossing model with a SPLINE (hinge) trajectory.
//
// The convex-quadratic trajectory of multistate_threshold.stan cannot rebound
// enough to reproduce molecular relapse. Here the latent trajectory gains a
// patient-specific LINEAR-SPLINE term: a slope change zeta_i after a knot kappa,
//
//   m_i(ell) = beta0 + beta1 ell + beta2 ell^2 + b0_i + b1_i ell
//              + zeta_i * fmax(0, ell - kappa),   zeta_i ~ N(0, sigma_zeta^2).
//
// This is piecewise-quadratic (two convex pieces joined at kappa), so the
// running minimum (onset) and post-knot values (relapse) remain CLOSED FORM
// -- no grid, no soft-min -- while zeta_i lets individual patients rebound
// after kappa and thus produce a documented relapse.
//==============================================================================
functions {
  // min of A l^2 + B l + C over [lo, hi], A > 0 (convex)
  real qmin(real A, real B, real C, real lo, real hi) {
    real lv = -B / (2 * A);
    real lc = fmin(fmax(lv, lo), hi);
    return A * lc^2 + B * lc + C;
  }
  // running minimum of the piecewise-quadratic trajectory over [0, lam]
  real runmin_pw(real a, real g, real b2, real ze, real kn, real lam) {
    real m1 = qmin(b2, g, a, 0, fmin(lam, kn));
    real m2 = (lam > kn) ? qmin(b2, g + ze, a - ze * kn, kn, lam)
                         : positive_infinity();
    return fmin(m1, m2);
  }
  // full trajectory value (handles the hinge)
  real mfull(real a, real g, real b2, real ze, real kn, real l) {
    return a + g * l + b2 * l^2 + ze * fmax(0, l - kn);
  }
}
data {
  int<lower=1> N;
  int<lower=1> Nobs;
  array[Nobs] int<lower=1,upper=N> pid;
  vector[Nobs] ell;
  vector[Nobs] y;
  array[Nobs] int<lower=0,upper=1> floor_ind;
  array[Nobs] int<lower=0,upper=1> bm;
  real cF;
  real cD;
  real<lower=0> W;
  real<lower=0> kappa;                     // spline knot on the ell scale (e.g. log(1+2))
  array[N] int<lower=0,upper=1> onset;
  vector[N] onL;
  vector[N] onR;
  vector[N] ellC;
  array[N] int<lower=0,upper=1> relapse;
  vector[N] rlL;
  vector[N] rlR;
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> beta2;
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=0> tau0;
  real<lower=0> tau1;
  real<lower=0> sigma_zeta;                // post-knot slope (rebound) SD
  real<lower=0.02> sigma_thr;              // small floor: spline explains crossings,
                                           // so sigma_thr is weakly identified and
                                           // otherwise funnels toward 0 (poor mixing).
  vector[N] z0;
  vector[N] z1;
  vector[N] zz;                            // non-centred rebound effects
}
transformed parameters {
  vector[N] b0   = z0 * tau0;
  vector[N] b1   = z1 * tau1;
  vector[N] zeta = zz * sigma_zeta;
}
model {
  beta0 ~ normal(-2.5, 2);
  beta1 ~ normal(-1, 1);
  beta2 ~ normal(0, 1);
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(1);
  tau0 ~ exponential(1);
  tau1 ~ exponential(1);
  sigma_zeta ~ normal(0, 2);               // half-normal; allows large rebound heterogeneity
  sigma_thr ~ normal(0, 0.5);
  z0 ~ std_normal();
  z1 ~ std_normal();
  zz ~ std_normal();

  // longitudinal assay model with left-censoring
  for (o in 1:Nobs) {
    real mu = mfull(beta0 + b0[pid[o]], beta1 + b1[pid[o]], beta2, zeta[pid[o]], kappa, ell[o])
              + beta_bm * bm[o];
    if (floor_ind[o] == 1)
      target += normal_lcdf(cF | mu, sigma_y);
    else
      target += normal_lpdf(y[o] | mu, sigma_y);
  }

  // transitions
  for (i in 1:N) {
    real a = beta0 + b0[i];
    real g = beta1 + b1[i];
    real ze = zeta[i];
    if (onset[i] == 1) {
      real MA = runmin_pw(a, g, beta2, ze, kappa, onL[i]);
      real MB = runmin_pw(a, g, beta2, ze, kappa, onR[i]);
      // clamped CDF difference (robust: never NaN/-inf if the two terms are
      // equal, saturated, or -- for relapse below -- out of order at init).
      target += log(fmax(Phi((MA - cD) / sigma_thr) - Phi((MB - cD) / sigma_thr), 1e-12));
      if (relapse[i] == 1) {
        real mR = mfull(a, g, beta2, ze, kappa, rlR[i]);
        real mL = mfull(a, g, beta2, ze, kappa, rlL[i]);
        target += log(fmax(Phi((mR - cD) / sigma_thr) - Phi((mL - cD) / sigma_thr), 1e-12));
      } else {
        real mC = mfull(a, g, beta2, ze, kappa, ellC[i]);
        target += normal_lcdf((cD - mC) / sigma_thr | 0, 1);
      }
    } else {
      real MC = runmin_pw(a, g, beta2, ze, kappa, ellC[i]);
      target += normal_lcdf((MC - cD) / sigma_thr | 0, 1);
    }
  }
}
generated quantities {
  vector[N] p_onset_byC;
  vector[N] p_relapse_byC;
  vector[Nobs] log_lik_long;
  for (o in 1:Nobs) {
    real mu = mfull(beta0 + b0[pid[o]], beta1 + b1[pid[o]], beta2, zeta[pid[o]], kappa, ell[o])
              + beta_bm * bm[o];
    log_lik_long[o] = (floor_ind[o] == 1) ? normal_lcdf(cF | mu, sigma_y)
                                          : normal_lpdf(y[o] | mu, sigma_y);
  }
  for (i in 1:N) {
    real a = beta0 + b0[i];
    real g = beta1 + b1[i];
    real ze = zeta[i];
    p_onset_byC[i]   = 1 - Phi((runmin_pw(a, g, beta2, ze, kappa, ellC[i]) - cD) / sigma_thr);
    p_relapse_byC[i] = Phi((mfull(a, g, beta2, ze, kappa, ellC[i]) - cD) / sigma_thr) * p_onset_byC[i];
  }
}
