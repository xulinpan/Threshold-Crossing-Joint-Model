//==============================================================================
// [SUPERSEDED / DEVELOPMENT VERSION]
// Single-state latent threshold-crossing joint model using a SOFT-MIN grid.
// Retained only for the recovery study; on the CML data it produced ~219
// divergent transitions from the soft-min geometry.
//
// >>> The PRIMARY model of this project is the multi-state model with an
//     analytic (piecewise-quadratic / linear-spline) trajectory:
//         multistate/stan/multistate_spline.stan
//     which computes the running minimum AND maximum in closed form (onset,
//     durable response, relapse), needs no soft-min grid, and samples without
//     divergences. Prefer that model for all applied analysis.
//
// Bayesian hierarchical latent threshold-crossing joint model
// for irregular molecular monitoring in CML.
//
// Longitudinal latent log-MRD trajectory with assay-floor LEFT-CENSORING,
// plus a probabilistic latent threshold-crossing likelihood for interval-
// observed deep molecular response (DMR). DMR onset is defined as the first
// time the latent trajectory crosses c_D, represented through the running
// minimum M_i(t) approximated by a differentiable soft-min over a grid.
//
// Reference: bayesian_hierarchical_threshold_crossing_model.tex (Sections 1-2).
//==============================================================================
data {
  // ---- longitudinal ----
  int<lower=1> N;                       // number of patients
  int<lower=1> Nobs;                    // number of longitudinal observations
  array[Nobs] int<lower=1,upper=N> pid; // patient index per observation
  vector[Nobs] ell;                     // log(1 + t) at each observation
  vector[Nobs] y;                       // observed log-MRD (floor rows hold c_F)
  array[Nobs] int<lower=0,upper=1> floor_ind; // 1 = assay-floor (left-censored)
  array[Nobs] int<lower=0,upper=1> bm;  // 1 = bone-marrow sample source
  real cF;                              // assay floor, e.g. -5.0
  real cD;                              // DMR threshold, e.g. -4.5

  // ---- DMR interval / threshold-crossing ----
  array[N] int<lower=0,upper=1> delta;  // 1 = DMR documented in (L_i,R_i]; 0 = censored
  int<lower=1> G;                       // grid size for running-minimum soft-min
  matrix[N,G] ellA;   // log(1+u) grid on [0, L_i]      (event) or [0, C_i] (censored)
  matrix[N,G] ellB;   // log(1+u) grid on [0, R_i]      (event) ; ignored if censored
  matrix[N,G] ellC;   // log(1+u) grid on [0, C_end_i]  (for predicted P(DMR by end))
  real<lower=0> kappa;                  // soft-min sharpness (e.g. 50)
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=0> tau0;
  real<lower=0> tau1;
  real<lower=0> sigma_thr;
  vector[N] z0;          // non-centered random intercept
  vector[N] z1;          // non-centered random slope
}
transformed parameters {
  vector[N] b0 = z0 * tau0;
  vector[N] b1 = z1 * tau1;
}
model {
  // ---- priors (weakly informative; see spec Section 2.1) ----
  beta0     ~ normal(-2.5, 2);
  beta1     ~ normal(-1, 1);
  beta2     ~ normal(0, 0.5);
  beta_bm   ~ normal(0, 1);
  sigma_y   ~ exponential(1);
  tau0      ~ exponential(1);
  tau1      ~ exponential(1);
  sigma_thr ~ normal(0, 0.25);          // Half-Normal(0, 0.25^2) via <lower=0>
  z0 ~ std_normal();
  z1 ~ std_normal();

  // ---- longitudinal assay model with left-censoring ----
  {
    vector[Nobs] mu;
    for (o in 1:Nobs)
      mu[o] = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
              + b0[pid[o]] + b1[pid[o]] * ell[o]
              + beta_bm * bm[o];
    for (o in 1:Nobs) {
      if (floor_ind[o] == 1)
        target += normal_lcdf(cF | mu[o], sigma_y);   // P(y* <= c_F)
      else
        target += normal_lpdf(y[o] | mu[o], sigma_y);
    }
  }

  // ---- latent threshold-crossing DMR likelihood ----
  for (i in 1:N) {
    vector[G] mA;
    vector[G] mB;
    for (g in 1:G) {
      mA[g] = beta0 + beta1 * ellA[i,g] + beta2 * square(ellA[i,g])
              + b0[i] + b1[i] * ellA[i,g];
      mB[g] = beta0 + beta1 * ellB[i,g] + beta2 * square(ellB[i,g])
              + b0[i] + b1[i] * ellB[i,g];
    }
    real MA = -log_sum_exp(-kappa * mA) / kappa;   // soft running minimum
    real MB = -log_sum_exp(-kappa * mB) / kappa;
    real a = (MA - cD) / sigma_thr;
    real b = (MB - cD) / sigma_thr;
    if (delta[i] == 1) {
      // log[ Phi(a) - Phi(b) ], with a >= b since M is non-increasing in t
      target += log_diff_exp(normal_lcdf(a | 0, 1), normal_lcdf(b | 0, 1));
    } else {
      target += normal_lcdf(a | 0, 1);             // P(T^D > C_i) = Phi((M(C)-cD)/sthr)
    }
  }
}
generated quantities {
  vector[Nobs] mu_obs;
  vector[Nobs] y_rep;
  vector[Nobs] p_floor;
  vector[N]    pDMR_end;    // predicted P(DMR by end of follow-up)

  for (o in 1:Nobs) {
    mu_obs[o] = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
                + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    y_rep[o]  = normal_rng(mu_obs[o], sigma_y);
    p_floor[o] = Phi((cF - mu_obs[o]) / sigma_y);
  }
  for (i in 1:N) {
    vector[G] mC;
    for (g in 1:G)
      mC[g] = beta0 + beta1 * ellC[i,g] + beta2 * square(ellC[i,g])
              + b0[i] + b1[i] * ellC[i,g];
    real MC = -log_sum_exp(-kappa * mC) / kappa;
    pDMR_end[i] = 1 - Phi((MC - cD) / sigma_thr);
  }
}
