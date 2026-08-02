//==============================================================================
// interval_hazard_joint_t.stan
//
// SCOPE. This is the SIMULATION-SANDBOX Student-t model, for use with
// dgm_common.R. It is NOT the file that produced the manuscript's Student-t
// sensitivity analysis: the application uses a different data contract, and
// that fit was run with 04_Code/Stan/glw_joint_interval_dmr_independent_studentt.stan
// against the real cohort. Do not use this file for anything reported in the
// paper.
//
// The joint interval-hazard model with a STUDENT-t observation model and
// Student-t random effects, in place of Gaussian. Everything else -- the
// assay-floor left-censoring, the complementary log-log interval hazard, the
// priors-as-data interface -- is identical to interval_hazard_joint_pp.stan.
//
// WHY THIS EXISTS. Section "Sensitivity to the censoring point and to serial
// correlation" concludes that the large residual SD is intrinsic to the cohort
// and that "heavier-tailed measurement error or a more flexible individual
// trajectory are the more likely explanations". The simulation study shows that
// under exactly that departure -- t3 errors and random effects, rescaled to the
// same standard deviation, so only tail shape differs -- coverage of sigma_y
// falls to 0.27 and of tau_b1 to 0.63. The paper therefore names heavy tails as
// the probable mechanism and then reports sigma_y and tau_b1 as though
// unaffected. Fitting this model answers the question directly: if sigma_y and
// tau_b1 move materially when the tails are freed, the reported heterogeneity is
// conditional on Gaussianity and must be described as such.
//
// nu is estimated with a gamma(2, 0.1) prior (mean 20), the usual weakly
// informative choice that neither forces nor forbids heavy tails; nu above
// roughly 30 is indistinguishable from Gaussian in practice. Report the
// posterior for nu alongside sigma_y and tau_b1.
//
// Data contract is that of interval_hazard_joint_pp.stan, minus prior_only.
//==============================================================================
data {
  int<lower=1> N;
  int<lower=1> Nobs;
  array[Nobs] int<lower=1,upper=N> pid;
  vector[Nobs] ell;
  vector[Nobs] y;
  array[Nobs] int<lower=0,upper=1> floor_ind;
  array[Nobs] int<lower=0,upper=1> bm;
  real cF;

  int<lower=1> Nint;
  array[Nint] int<lower=1,upper=N> pid_int;
  vector[Nint] midlog;
  vector[Nint] gaplog;
  vector<lower=0>[Nint] delta_len;
  array[Nint] int<lower=0,upper=1> event;

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

  int<lower=0,upper=1> gq;
  // 1 = also give the random effects t tails; 0 = Gaussian random effects with
  // t observation error only. Fit both to see which layer carries the tails.
  int<lower=0,upper=1> t_random_effects;
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
  real<lower=2> nu;            // degrees of freedom; >2 so the variance exists
  vector[N] z0;
  vector[N] z1;
}
transformed parameters {
  vector[N] b0 = z0 * tau0;
  vector[N] b1 = z1 * tau1;
}
model {
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
  nu      ~ gamma(2, 0.1);

  if (t_random_effects == 1) {
    z0 ~ student_t(nu, 0, 1);
    z1 ~ student_t(nu, 0, 1);
  } else {
    z0 ~ std_normal();
    z1 ~ std_normal();
  }

  // ---- longitudinal, Student-t with assay-floor left-censoring ----
  for (o in 1:Nobs) {
    real mu = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
              + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
    if (floor_ind[o] == 1)
      target += student_t_lcdf(cF | nu, mu, sigma_y);
    else
      target += student_t_lpdf(y[o] | nu, mu, sigma_y);
  }

  // ---- interval hazard, unchanged ----
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
generated quantities {
  // sigma_y is the t SCALE, not the SD. Report both: the marginal observation SD
  // is sigma_y * sqrt(nu/(nu-2)), and that is the quantity comparable with the
  // Gaussian fit's sigma_y.
  real sd_marginal = sigma_y * sqrt(nu / (nu - 2));
  vector[gq == 1 ? Nobs : 0] y_rep;
  vector[gq == 1 ? Nobs : 0] p_floor;
  if (gq == 1) {
    for (o in 1:Nobs) {
      real mu = beta0 + beta1 * ell[o] + beta2 * square(ell[o])
                + b0[pid[o]] + b1[pid[o]] * ell[o] + beta_bm * bm[o];
      y_rep[o]   = student_t_rng(nu, mu, sigma_y);
      p_floor[o] = student_t_cdf(cF | nu, mu, sigma_y);
    }
  }
}
