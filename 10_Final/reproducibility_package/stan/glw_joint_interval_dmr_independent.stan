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
      target += normal_lcdf(floor_value | mu, sigma_y);
    } else {
      y[n] ~ normal(mu, sigma_y);
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

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real lt = log1p(t_obs[n]);
    real mu = beta0 + beta_time * lt + beta_time2 * square(lt) +
      beta_bm * sample_bm[n] + b[i, 1] + b[i, 2] * lt;

    if (is_floor[n] == 1) {
      log_lik_y[n] = normal_lcdf(floor_value | mu, sigma_y);
    } else {
      log_lik_y[n] = normal_lpdf(y[n] | mu, sigma_y);
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
