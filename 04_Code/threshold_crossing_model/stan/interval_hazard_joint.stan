//==============================================================================
// Interval-hazard joint model (pragmatic comparator).
// Targets FIRST DOCUMENTED DMR under the observed visit process, using a
// complementary log-log interval hazard driven by the shared latent
// log-MRD trajectory. Longitudinal part is identical to the threshold-
// crossing model (assay-floor left-censoring).
//
// Reference: bayesian_hierarchical_threshold_crossing_model.tex (Sec 1.6).
//==============================================================================
data {
  // ---- longitudinal ----
  int<lower=1> N;
  int<lower=1> Nobs;
  array[Nobs] int<lower=1,upper=N> pid;
  vector[Nobs] ell;                     // log(1+t)
  vector[Nobs] y;
  array[Nobs] int<lower=0,upper=1> floor_ind;
  array[Nobs] int<lower=0,upper=1> bm;
  real cF;

  // ---- at-risk intervals for documented DMR ----
  int<lower=1> Nint;
  array[Nint] int<lower=1,upper=N> pid_int;
  vector[Nint] midlog;                  // log(1 + t_mid)
  vector[Nint] gaplog;                  // log(1 + Delta)
  vector<lower=0>[Nint] delta_len;      // interval length Delta (years)
  array[Nint] int<lower=0,upper=1> event; // 1 = DMR documented in this interval
}
parameters {
  real beta0;
  real beta1;
  real beta2;
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=0> tau0;
  real<lower=0> tau1;
  real gamma0;
  real gamma1;
  real gamma2;
  real alpha;
  vector[N] z0;
  vector[N] z1;
}
transformed parameters {
  vector[N] b0 = z0 * tau0;
  vector[N] b1 = z1 * tau1;
}
model {
  // ---- priors ----
  beta0   ~ normal(-2.5, 2);
  beta1   ~ normal(-1, 1);
  beta2   ~ normal(0, 0.5);
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(1);
  tau0    ~ exponential(1);
  tau1    ~ exponential(1);
  gamma0  ~ normal(-2, 2);
  gamma1  ~ normal(0, 1);
  gamma2  ~ normal(0, 1);
  alpha   ~ normal(-0.5, 0.75);
  z0 ~ std_normal();
  z1 ~ std_normal();

  // ---- longitudinal ----
  for (o in 1:Nobs) {
    real mu = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
              + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    if (floor_ind[o] == 1)
      target += normal_lcdf(cF | mu, sigma_y);
    else
      target += normal_lpdf(y[o] | mu, sigma_y);
  }

  // ---- interval hazard for first documented DMR ----
  for (k in 1:Nint) {
    real m = beta0 + beta1 * midlog[k] + beta2 * square(midlog[k])
             + b0[pid_int[k]] + b1[pid_int[k]] * midlog[k];
    real logh = gamma0 + gamma1 * midlog[k] + gamma2 * gaplog[k] + alpha * m;
    real Hd = exp(logh) * delta_len[k];               // cumulative hazard over interval
    if (event[k] == 1)
      target += log1m_exp(-Hd);                       // log p_ik
    else
      target += -Hd;                                  // log(1 - p_ik)
  }
}
generated quantities {
  vector[Nobs] mu_obs;
  vector[Nobs] y_rep;
  vector[Nobs] p_floor;
  for (o in 1:Nobs) {
    mu_obs[o] = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
                + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    y_rep[o]  = normal_rng(mu_obs[o], sigma_y);
    p_floor[o] = Phi((cF - mu_obs[o]) / sigma_y);
  }
}
