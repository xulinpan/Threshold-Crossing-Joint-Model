// glw_joint_interval_dmr_studentt.stan
//
// Student-t variant of glw_joint_interval_dmr.stan. IDENTICAL data contract, so
// it drops into 06_fit_stan_joint_model.R by changing the stan_file path only.
// The sole change is the longitudinal observation model: Student-t in place of
// Gaussian, with the degrees of freedom nu estimated.
//
// WHY. sigma_y is reported at ~1.81 on the log10 scale, against a marginal SD of
// the observed log-MRD of 1.83 and an analytical assay CV of order 0.2-0.3. It
// is about a third of the biomarker's dynamic range. Because 48% of
// observations sit at the assay floor, the inferred latent spread is driven
// substantially by the ASSUMED SHAPE of the latent distribution below -5, where
// no data are observed -- so the Gaussian assumption is doing real work, not
// merely providing a convenient error term. Simulation shows that under t3
// errors rescaled to the same SD, coverage of sigma_y falls to 0.27 and of the
// random-slope SD to 0.63. This fit determines how much of the reported
// between-patient heterogeneity survives when the tails are freed.
//
// WHAT TO COMPARE. sigma_y here is the t SCALE, not a standard deviation. The
// comparable quantity is sd_marginal = sigma_y * sqrt(nu/(nu-2)), reported in
// generated quantities. Compare that, tau_b, and alpha_mrd against the Gaussian
// fit. nu above roughly 30 is effectively Gaussian; nu near its lower bound
// indicates the tails were carrying the variance.
//
// The prior nu ~ gamma(2, 0.1) has mean 20 and is the conventional weakly
// informative choice: it neither forces heavy tails nor excludes them.
data {
  int<lower=1> N_obs;
  int<lower=1> N_pat;
  int<lower=1> N_int;
  array[N_obs] int<lower=1, upper=N_pat> id_obs;
  array[N_int] int<lower=1, upper=N_pat> id_int;
  vector[N_obs] y;
  array[N_obs] int<lower=0, upper=1> is_floor;
  real floor_value;
  vector<lower=0>[N_obs] t_obs;
  vector<lower=0>[N_int] t_start;
  vector<lower=0>[N_int] t_end;
  vector<lower=0>[N_int] gap;
  vector[N_obs] sample_bm;
  array[N_int] int<lower=0, upper=1> event_interval;
}

parameters {
  real beta0;
  real beta_time;
  real beta_time2;
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=2> nu;                 // > 2 so the variance exists

  real gamma0;
  real gamma_time;
  real gamma_gap;
  real alpha_mrd;

  vector<lower=0>[2] tau_b;
  cholesky_factor_corr[2] L_b;
  matrix[2, N_pat] z_b;
}

transformed parameters {
  matrix[N_pat, 2] b;
  b = (diag_pre_multiply(tau_b, L_b) * z_b)';
}

model {
  beta0 ~ normal(-2.5, 2);
  beta_time ~ normal(-1, 1);
  beta_time2 ~ normal(0, 0.5);
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(1);
  nu ~ gamma(2, 0.1);

  gamma0 ~ normal(-2, 2);
  gamma_time ~ normal(0, 1);
  gamma_gap ~ normal(0, 1);
  alpha_mrd ~ normal(-0.5, 0.75);

  tau_b ~ exponential(1);
  L_b ~ lkj_corr_cholesky(2);
  to_vector(z_b) ~ normal(0, 1);

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real lt = log1p(t_obs[n]);
    real mu = beta0 + beta_time * lt + beta_time2 * square(lt) +
      beta_bm * sample_bm[n] + b[i, 1] + b[i, 2] * lt;

    if (is_floor[n] == 1) {
      target += student_t_lcdf(floor_value | nu, mu, sigma_y);
    } else {
      y[n] ~ student_t(nu, mu, sigma_y);
    }
  }

  // event sub-model unchanged
  for (m in 1:N_int) {
    int i = id_int[m];
    real t_mid = 0.5 * (t_start[m] + t_end[m]);
    real lt_mid = log1p(t_mid);
    real mu_mid = beta0 + beta_time * lt_mid + beta_time2 * square(lt_mid) +
      b[i, 1] + b[i, 2] * lt_mid;
    real eta = gamma0 + gamma_time * lt_mid + gamma_gap * log1p(gap[m]) +
      alpha_mrd * mu_mid;
    real cum_hazard = exp(eta) * fmax(gap[m], 1e-6);

    if (event_interval[m] == 1) {
      target += log1m_exp(-cum_hazard);
    } else {
      target += -cum_hazard;
    }
  }
}

generated quantities {
  corr_matrix[2] Omega_b;
  vector[N_obs] log_lik_y;
  vector[N_int] event_prob;
  // the quantity comparable with the Gaussian fit's sigma_y
  real sd_marginal = sigma_y * sqrt(nu / (nu - 2));

  Omega_b = multiply_lower_tri_self_transpose(L_b);

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real lt = log1p(t_obs[n]);
    real mu = beta0 + beta_time * lt + beta_time2 * square(lt) +
      beta_bm * sample_bm[n] + b[i, 1] + b[i, 2] * lt;

    if (is_floor[n] == 1) {
      log_lik_y[n] = student_t_lcdf(floor_value | nu, mu, sigma_y);
    } else {
      log_lik_y[n] = student_t_lpdf(y[n] | nu, mu, sigma_y);
    }
  }

  for (m in 1:N_int) {
    int i = id_int[m];
    real t_mid = 0.5 * (t_start[m] + t_end[m]);
    real lt_mid = log1p(t_mid);
    real mu_mid = beta0 + beta_time * lt_mid + beta_time2 * square(lt_mid) +
      b[i, 1] + b[i, 2] * lt_mid;
    real eta = gamma0 + gamma_time * lt_mid + gamma_gap * log1p(gap[m]) +
      alpha_mrd * mu_mid;
    real cum_hazard = exp(eta) * fmax(gap[m], 1e-6);

    event_prob[m] = 1 - exp(-cum_hazard);
  }
}
