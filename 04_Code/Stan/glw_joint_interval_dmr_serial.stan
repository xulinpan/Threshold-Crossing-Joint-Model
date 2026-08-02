// =====================================================================
// glw_joint_interval_dmr_serial.stan
//
// Audit item A2/B1: the primary model absorbs within-patient biological
// fluctuation into sigma_y (~1.8 on log10, far above assay CV). This model
// adds a within-patient continuous-time AR(1) / Ornstein-Uhlenbeck latent
// residual u(t) to the OBSERVATION mean only, so sigma_y becomes pure
// measurement error while u(t) carries the serial fluctuation. The latent
// biological trajectory that drives the interval hazard is left smooth.
//
// Observations MUST be sorted by (patient, time) so that prev_obs[n] < n.
// The companion script 14_fit_serial_correlation.R builds first_in_patient,
// prev_obs and dt_prev and reorders all obs-level vectors accordingly.
//
// STATUS: NEW MODEL — compile and check on simulated data before trusting
// the real-cohort results. Recommended adapt_delta >= 0.99.
// =====================================================================
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

  // serial-correlation bookkeeping (built in R; obs sorted by patient,time)
  array[N_obs] int<lower=0, upper=1> first_in_patient;
  array[N_obs] int<lower=0, upper=N_obs> prev_obs;   // 0 if first in patient
  vector<lower=0>[N_obs] dt_prev;                     // years since prev visit
}

parameters {
  real beta0;
  real beta_time;
  real beta_time2;
  real beta_bm;
  real<lower=0> sigma_y;         // now pure measurement error

  real gamma0;
  real gamma_time;
  real gamma_gap;
  real alpha_mrd;

  vector<lower=0>[2] tau_b;
  matrix[2, N_pat] z_b;

  real<lower=0> sigma_u;         // OU marginal SD (serial fluctuation)
  real<lower=0> ell;             // OU correlation length (years)
  vector[N_obs] z_u;             // non-centered OU innovations
}

transformed parameters {
  matrix[N_pat, 2] b;
  vector[N_obs] u;               // within-patient OU residual
  for (i in 1:N_pat) {
    b[i, 1] = tau_b[1] * z_b[1, i];
    b[i, 2] = tau_b[2] * z_b[2, i];
  }
  for (n in 1:N_obs) {
    if (first_in_patient[n] == 1) {
      u[n] = sigma_u * z_u[n];
    } else {
      real rho = exp(-dt_prev[n] / ell);
      u[n] = rho * u[prev_obs[n]] + sigma_u * sqrt(1 - square(rho)) * z_u[n];
    }
  }
}

model {
  beta0 ~ normal(-2.5, 2);
  beta_time ~ normal(-1, 1);
  beta_time2 ~ normal(0, 0.5);
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(2);      // encourage small pure-assay error
  gamma0 ~ normal(-2, 2);
  gamma_time ~ normal(0, 1);
  gamma_gap ~ normal(0, 1);
  alpha_mrd ~ normal(-0.5, 0.75);
  tau_b ~ exponential(1);
  to_vector(z_b) ~ normal(0, 1);

  sigma_u ~ exponential(1);
  ell ~ lognormal(0, 1);         // ~ years; prior mass on months-to-years
  z_u ~ normal(0, 1);

  // longitudinal: smooth mean + serial residual, then measurement error
  for (n in 1:N_obs) {
    int i = id_obs[n];
    real lt = log1p(t_obs[n]);
    real mu = beta0 + beta_time * lt + beta_time2 * square(lt) +
      beta_bm * sample_bm[n] + b[i, 1] + b[i, 2] * lt + u[n];
    if (is_floor[n] == 1)
      target += normal_lcdf(floor_value | mu, sigma_y);
    else
      y[n] ~ normal(mu, sigma_y);
  }

  // interval hazard uses the SMOOTH latent trajectory (no serial residual)
  for (m in 1:N_int) {
    int i = id_int[m];
    real t_mid = 0.5 * (t_start[m] + t_end[m]);
    real lt_mid = log1p(t_mid);
    real mu_mid = beta0 + beta_time * lt_mid + beta_time2 * square(lt_mid) +
      b[i, 1] + b[i, 2] * lt_mid;
    real eta = gamma0 + gamma_time * lt_mid + gamma_gap * log1p(gap[m]) +
      alpha_mrd * mu_mid;
    real cum_hazard = exp(eta) * fmax(gap[m], 1e-6);
    if (event_interval[m] == 1)
      target += log1m_exp(-cum_hazard);
    else
      target += -cum_hazard;
  }
}

generated quantities {
  real sigma_total = sqrt(square(sigma_y) + square(sigma_u)); // marginal obs SD
}
