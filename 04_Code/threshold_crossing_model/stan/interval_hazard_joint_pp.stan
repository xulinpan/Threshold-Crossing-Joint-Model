//==============================================================================
// interval_hazard_joint_pp.stan
// Identical likelihood to interval_hazard_joint.stan, with the PRIORS PASSED AS
// DATA so a sensitivity analysis can vary them without editing the model, and a
// prior_only switch for prior-predictive checks.
//
// WHY. The simulation study shows bias in exactly the parameters whose values
// sit far in the tails of the fixed priors, all of it directed toward the prior
// means: gamma2 (prior N(0,1)) and beta1 (prior N(-1,1)). For gamma2 the
// coverage does not recover as n grows. A reviewer will ask how much of the
// reported inference is prior rather than data, and that question cannot be
// answered while the priors are hard-coded.
//
// DEFAULTS. Passing the hyperparameters below reproduces
// interval_hazard_joint.stan exactly:
//   pr_beta0   = [-2.5, 2]     pr_gamma0 = [-2, 2]
//   pr_beta1   = [-1, 1]       pr_gamma1 = [0, 1]
//   pr_beta2   = [0, 0.5]      pr_gamma2 = [0, 1]
//   pr_beta_bm = [0, 1]        pr_alpha  = [-0.5, 0.75]
//   pr_sigma_y = 1  pr_tau0 = 1  pr_tau1 = 1     (exponential rates)
//   prior_only = 0
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
  array[Nint] int<lower=0,upper=1> event;

  // ---- priors as data: [location, scale] for normals, rate for exponentials --
  vector[2] pr_beta0;
  vector[2] pr_beta1;
  vector[2] pr_beta2;
  vector[2] pr_beta_bm;
  real<lower=0> pr_sigma_y;
  real<lower=0> pr_tau0;
  real<lower=0> pr_tau1;
  vector[2] pr_gamma0;
  vector[2] pr_gamma1;
  vector[2] pr_gamma2;
  vector[2] pr_alpha;

  // 1 = ignore the likelihood and sample the prior (prior-predictive check)
  int<lower=0,upper=1> prior_only;

  // 1 = emit the posterior-predictive quantities below; 0 = skip them.
  // These are per-observation and per-interval, so at n=150 they add roughly
  // 4,650 columns to every draw -- about 43 MB per chain, written to disk and
  // then parsed back into R. A study that only needs the 11 scalar parameters
  // (e.g. 15_prior_sensitivity.R) should pass gq = 0.
  int<lower=0,upper=1> gq;
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
  // ---- priors (identical in form to the fixed-prior model) ----
  beta0   ~ normal(pr_beta0[1], pr_beta0[2]);
  beta1   ~ normal(pr_beta1[1], pr_beta1[2]);
  beta2   ~ normal(pr_beta2[1], pr_beta2[2]);
  beta_bm ~ normal(pr_beta_bm[1], pr_beta_bm[2]);
  sigma_y ~ exponential(pr_sigma_y);
  tau0    ~ exponential(pr_tau0);
  tau1    ~ exponential(pr_tau1);
  gamma0  ~ normal(pr_gamma0[1], pr_gamma0[2]);
  gamma1  ~ normal(pr_gamma1[1], pr_gamma1[2]);
  gamma2  ~ normal(pr_gamma2[1], pr_gamma2[2]);
  alpha   ~ normal(pr_alpha[1], pr_alpha[2]);
  z0 ~ std_normal();
  z1 ~ std_normal();

  if (prior_only == 0) {
    // ---- longitudinal, with assay-floor left-censoring ----
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
      real Hd = exp(logh) * delta_len[k];
      if (event[k] == 1)
        target += log1m_exp(-Hd);
      else
        target += -Hd;
    }
  }
}
generated quantities {
  // Zero-length when gq = 0, so nothing is written to the sampler CSVs.
  vector[gq == 1 ? Nobs : 0] mu_obs;
  vector[gq == 1 ? Nobs : 0] y_rep;
  vector[gq == 1 ? Nobs : 0] p_floor;
  vector[gq == 1 ? Nint : 0] p_event;   // predicted P(DMR in interval k)
  if (gq == 1) {
    for (o in 1:Nobs) {
      mu_obs[o] = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
                  + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
      y_rep[o]  = normal_rng(mu_obs[o], sigma_y);
      p_floor[o] = Phi((cF - mu_obs[o]) / sigma_y);
    }
    for (k in 1:Nint) {
      real m = beta0 + beta1 * midlog[k] + beta2 * square(midlog[k])
               + b0[pid_int[k]] + b1[pid_int[k]] * midlog[k];
      real logh = gamma0 + gamma1 * midlog[k] + gamma2 * gaplog[k] + alpha * m;
      p_event[k] = 1 - exp(-exp(logh) * delta_len[k]);
    }
  }
}
