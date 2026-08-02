functions {
  real soft_min_segment(vector x, real temperature) {
    return -temperature * log_sum_exp(-x / temperature);
  }

  real threshold_event_lpmf(
    int event_observed,
    real min_to_left,
    real min_to_right,
    real threshold,
    real sigma_threshold,
    real deterministic_eps,
    int threshold_model
  ) {
    if (threshold_model == 1) {
      if (event_observed == 1) {
        return log_inv_logit((min_to_left - threshold) / deterministic_eps) +
          log_inv_logit((threshold - min_to_right) / deterministic_eps);
      } else {
        return log_inv_logit((min_to_right - threshold) / deterministic_eps);
      }
    } else {
      if (event_observed == 1) {
        real log_F_left = normal_lcdf(min_to_left | threshold, sigma_threshold);
        real log_F_right = normal_lcdf(min_to_right | threshold, sigma_threshold);
        return log_diff_exp(log_F_left, fmin(log_F_right, log_F_left - 1e-12));
      } else {
        return normal_lcdf(min_to_right | threshold, sigma_threshold);
      }
    }
  }
}

data {
  int<lower=1> N_obs;
  int<lower=1> N_pat;
  int<lower=1> N_grid;

  array[N_obs] int<lower=1, upper=N_pat> id_obs;
  vector<lower=0>[N_obs] t_obs;
  vector[N_obs] sample_bm;

  // obs_type: 1 exact continuous, 2 left-censored, 3 interval-censored.
  array[N_obs] int<lower=1, upper=3> obs_type;
  vector[N_obs] y;
  vector[N_obs] y_lower;
  vector[N_obs] y_upper;

  // A patient-specific grid should include baseline, observed visit times,
  // event-interval endpoints, censoring times, and optional intermediate points.
  array[N_grid] int<lower=1, upper=N_pat> id_grid;
  vector<lower=0>[N_grid] t_grid;
  array[N_pat] int<lower=1> grid_start;
  array[N_pat] int<lower=1> idx_left_end;
  array[N_pat] int<lower=1> idx_right_end;

  array[N_pat] int<lower=0, upper=1> dmr_observed;

  real dmr_threshold;
  real floor_value;
  real<lower=0> softmin_temperature;
  real<lower=0> deterministic_eps;
  int<lower=1, upper=2> threshold_model;
}

parameters {
  real beta0;
  real beta_time;
  real beta_time2;
  real beta_source_bm;
  real<lower=0> sigma_y;

  vector<lower=0>[2] tau_b;
  matrix[2, N_pat] z_b;

  real<lower=0> sigma_threshold;
}

transformed parameters {
  matrix[N_pat, 2] b;
  vector[N_grid] eta_grid;

  for (i in 1:N_pat) {
    b[i, 1] = tau_b[1] * z_b[1, i];
    b[i, 2] = tau_b[2] * z_b[2, i];
  }

  for (m in 1:N_grid) {
    int i = id_grid[m];
    real lt = log1p(t_grid[m]);
    eta_grid[m] = beta0 + beta_time * lt + beta_time2 * square(lt) +
      b[i, 1] + b[i, 2] * lt;
  }
}

model {
  beta0 ~ normal(-2.5, 2);
  beta_time ~ normal(-1, 1);
  beta_time2 ~ normal(0, 0.5);
  beta_source_bm ~ normal(0, 0.75);
  sigma_y ~ exponential(1);

  tau_b ~ exponential(1);
  to_vector(z_b) ~ normal(0, 1);

  // Weakly informative but regularizing: threshold uncertainty is on the
  // log10 MRD scale and should not absorb ordinary measurement error.
  sigma_threshold ~ normal(0, 0.20);

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real lt = log1p(t_obs[n]);
    real eta_bio = beta0 + beta_time * lt + beta_time2 * square(lt) +
      b[i, 1] + b[i, 2] * lt;
    real mu_obs = eta_bio + beta_source_bm * sample_bm[n];

    if (obs_type[n] == 1) {
      y[n] ~ normal(mu_obs, sigma_y);
    } else if (obs_type[n] == 2) {
      target += normal_lcdf(floor_value | mu_obs, sigma_y);
    } else {
      target += log_diff_exp(
        normal_lcdf(y_upper[n] | mu_obs, sigma_y),
        normal_lcdf(y_lower[n] | mu_obs, sigma_y)
      );
    }
  }

  for (i in 1:N_pat) {
    int n_left = idx_left_end[i] - grid_start[i] + 1;
    int n_right = idx_right_end[i] - grid_start[i] + 1;
    real min_left = soft_min_segment(
      segment(eta_grid, grid_start[i], n_left),
      softmin_temperature
    );
    real min_right = soft_min_segment(
      segment(eta_grid, grid_start[i], n_right),
      softmin_temperature
    );

    target += threshold_event_lpmf(
      dmr_observed[i] |
      min_left,
      min_right,
      dmr_threshold,
      sigma_threshold,
      deterministic_eps,
      threshold_model
    );
  }
}

generated quantities {
  vector[N_pat] posterior_crossed_by_right;

  for (i in 1:N_pat) {
    int n_right = idx_right_end[i] - grid_start[i] + 1;
    real min_right = soft_min_segment(
      segment(eta_grid, grid_start[i], n_right),
      softmin_temperature
    );
    posterior_crossed_by_right[i] = inv_logit((dmr_threshold - min_right) / deterministic_eps);
  }
}
