// glw_joint_interval_dmr_independent_studentt.stan
//
// Student-t variant of glw_joint_interval_dmr_INDEPENDENT.stan -- the model that
// produces Table 5 and every headline number in the manuscript.
//
// An earlier variant was built from glw_joint_interval_dmr.stan, which uses
// CORRELATED random effects (lkj_corr_cholesky). That model is not the primary
// analysis: Table 5's caption specifies the independent-random-effects model,
// and its posterior means differ (alpha_mrd -1.275 vs -1.256, beta_time -3.561
// vs -3.539). A sensitivity analysis must be run against the model whose results
// are reported, so this file mirrors the independent specification exactly.
//
// Identical data contract, so it drops into 06_fit_stan_joint_model_independent.R
// by changing the stan_file path only. The sole change from that model is the
// longitudinal observation model: Student-t with nu estimated, in place of
// Gaussian. Priors, random-effect structure and event sub-model are untouched.
//
// WHAT TO COMPARE. sigma_y is a t SCALE here, not a standard deviation, and with
// nu near 2 the implied variance is barely defined -- do not compare sigma_y
// directly to the Gaussian sigma_y as if both were SDs. Compare the scale, nu,
// tau_b[1], tau_b[2], alpha_mrd and beta_time2 against the Gaussian fit.
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
  real<lower=2> nu;

  real gamma0;
  real gamma_time;
  real gamma_gap;
  real alpha_mrd;

  vector<lower=0>[2] tau_b;
  matrix[2, N_pat] z_b;
}

transformed parameters {
  matrix[N_pat, 2] b;
  for (i in 1:N_pat) {
    b[i, 1] = tau_b[1] * z_b[1, i];
    b[i, 2] = tau_b[2] * z_b[2, i];
  }
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
  vector[N_obs] log_lik_y;
  vector[N_int] event_prob;
  real sd_marginal = sigma_y * sqrt(nu / (nu - 2));

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
