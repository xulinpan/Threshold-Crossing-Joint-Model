//==============================================================================
// Bidirectional multi-state latent threshold-crossing joint model.
//
// States for molecular response depth relative to the DMR threshold c_D:
//   0 = not yet in DMR   (latent m has not crossed c_D)
//   1 = in DMR           (m has crossed below c_D)
//   2 = molecular relapse (m has crossed back above c_D after onset)
//
// The latent log-MRD trajectory is quadratic in ell = log(1+t):
//   m_i(ell) = a_i + g_i * ell + beta2 * ell^2,  a_i = beta0 + b0_i, g_i = beta1 + b1_i.
// With beta2 > 0 (convex) the trajectory declines then rebounds, so it crosses
// c_D downward once (ONSET) and upward once (RELAPSE). Both crossings are
// interval-observed at irregular visits and represented probabilistically with
// threshold uncertainty sigma_thr (distinct from measurement error sigma_y, so
// the two are separately identified).
//
// Key device: the running minimum M_i(t) = min_{u<=t} m_i(u) has a CLOSED FORM
// for a convex quadratic (value at the clamped vertex) -- no soft-min grid, so
// the geometry is smooth and divergence-prone soft-min ridges are avoided.
// After the vertex the trajectory is monotone increasing, so the upward
// (relapse) crossing uses the point value m_i(t) directly, and durability
// (sustained DMR over a window w) is m_i(t_on + w) <= c_D -- a running-MAX
// condition that reduces to a point value.
//==============================================================================
functions {
  // running minimum of a_i + g*ell + b2*ell^2 over ell in [0, lam], b2 > 0
  real runmin(real a, real g, real b2, real lam) {
    real lv = -g / (2 * b2);            // vertex
    real lc = fmin(fmax(lv, 0), lam);   // clamp to [0, lam]
    return a + g * lc + b2 * lc^2;
  }
  real mquad(real a, real g, real b2, real l) { return a + g * l + b2 * l^2; }
}
data {
  // ---- longitudinal ----
  int<lower=1> N;
  int<lower=1> Nobs;
  array[Nobs] int<lower=1,upper=N> pid;
  vector[Nobs] ell;                       // log(1+t)
  vector[Nobs] y;
  array[Nobs] int<lower=0,upper=1> floor_ind;
  array[Nobs] int<lower=0,upper=1> bm;
  real cF;                                // assay floor (-5.0)
  real cD;                                // DMR threshold (-4.5)
  real<lower=0> W;                        // durability window (years) for reporting

  // ---- transition information per patient (ell scale) ----
  array[N] int<lower=0,upper=1> onset;    // 1 if DMR documented in (onL,onR]
  vector[N] onL;                          // ell at last visit before onset (or C_end if censored)
  vector[N] onR;                          // ell at documenting visit  (ignored if onset==0)
  vector[N] ellC;                         // ell at end of follow-up (C_end)
  array[N] int<lower=0,upper=1> relapse;  // 1 if relapse documented in (rlL,rlR] (only if onset)
  vector[N] rlL;                          // ell at last visit before relapse
  vector[N] rlR;                          // ell at documenting relapse visit
}
parameters {
  real beta0;
  real beta1;
  real<lower=0> beta2;                    // convex curvature (>0)
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=0> tau0;
  real<lower=0> tau1;
  real<lower=0> sigma_thr;
  real<lower=0> sigma_rho;   // relapse frailty SD (0 -> relapse width = sigma_thr)
  vector[N] z0;
  vector[N] z1;
}
transformed parameters {
  vector[N] b0 = z0 * tau0;
  vector[N] b1 = z1 * tau1;
}
model {
  // ---- priors ----
  beta0 ~ normal(-2.5, 2);
  beta1 ~ normal(-1, 1);
  beta2 ~ normal(0, 1);                   // half-normal via <lower=0>
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(1);
  tau0 ~ exponential(1);
  tau1 ~ exponential(1);
  sigma_thr ~ normal(0, 0.5);             // half-normal
  sigma_rho ~ normal(0, 0.5);             // relapse frailty (half-normal)
  z0 ~ std_normal();
  z1 ~ std_normal();

  // ---- longitudinal assay model with left-censoring ----
  for (o in 1:Nobs) {
    real mu = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
              + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    if (floor_ind[o] == 1)
      target += normal_lcdf(cF | mu, sigma_y);
    else
      target += normal_lpdf(y[o] | mu, sigma_y);
  }

  // ---- transition likelihoods ----
  // Gaussian relapse frailty marginalizes to a wider relapse-specific width:
  //   int Phi((x+rho)/sigma_thr) N(rho;0,sigma_rho) drho = Phi(x / sigma_rel).
  real sigma_rel = sqrt(square(sigma_thr) + square(sigma_rho));
  for (i in 1:N) {
    real a = beta0 + b0[i];
    real g = beta1 + b1[i];
    // ONSET (state 0 -> 1) via running minimum
    if (onset[i] == 1) {
      real MA = runmin(a, g, beta2, onL[i]);
      real MB = runmin(a, g, beta2, onR[i]);
      // clamped CDF difference (robust to equal/saturated/out-of-order terms)
      target += log(fmax(Phi((MA - cD) / sigma_thr) - Phi((MB - cD) / sigma_thr), 1e-12));
    } else {
      real MC = runmin(a, g, beta2, ellC[i]);
      target += normal_lcdf((MC - cD) / sigma_thr | 0, 1);   // no onset by C
    }
    // RELAPSE (state 1 -> 2) via post-vertex point value (monotone increasing)
    if (onset[i] == 1) {
      if (relapse[i] == 1) {
        real mR = mquad(a, g, beta2, rlR[i]);
        real mL = mquad(a, g, beta2, rlL[i]);
        target += log(fmax(Phi((mR - cD) / sigma_rel) - Phi((mL - cD) / sigma_rel), 1e-12));
      } else {
        real mC = mquad(a, g, beta2, ellC[i]);
        target += normal_lcdf((cD - mC) / sigma_rel | 0, 1); // no relapse by C
      }
    }
  }
}
generated quantities {
  vector[N] p_onset_byC;     // P(onset by end of follow-up)
  vector[N] p_relapse_byC;   // P(relapse by end of follow-up | pathwise)
  vector[N] p_sustained_W;   // P(sustained DMR through W years after onset)
  vector[Nobs] log_lik_long;
  real sigma_rel = sqrt(square(sigma_thr) + square(sigma_rho));
  for (o in 1:Nobs) {
    real mu = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
              + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    log_lik_long[o] = (floor_ind[o] == 1) ? normal_lcdf(cF | mu, sigma_y)
                                          : normal_lpdf(y[o] | mu, sigma_y);
  }
  for (i in 1:N) {
    real a = beta0 + b0[i];
    real g = beta1 + b1[i];
    real MC = runmin(a, g, beta2, ellC[i]);
    p_onset_byC[i]   = 1 - Phi((MC - cD) / sigma_thr);
    real mC = mquad(a, g, beta2, ellC[i]);
    p_relapse_byC[i] = Phi((mC - cD) / sigma_rel) * p_onset_byC[i];
    // sustained: latent value W years after the (approx) onset time still <= cD.
    // onset ell approximated by the smaller root of m(ell)=cD when it exists.
    real disc = g^2 - 4 * beta2 * (a - cD);
    if (disc > 0) {
      real lon = (-g - sqrt(disc)) / (2 * beta2);
      real ton = expm1(fmax(lon, 0));
      real lW  = log1p(ton + W);
      p_sustained_W[i] = Phi((cD - mquad(a, g, beta2, lW)) / sigma_rel) * p_onset_byC[i];
    } else {
      p_sustained_W[i] = 0;
    }
  }
}
