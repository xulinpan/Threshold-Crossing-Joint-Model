// =====================================================================
// glw_unified_joint.stan  —  ROUTE B: one model, flag-switched likelihood
//
// A single latent trajectory and one set of random effects serve every
// endpoint in the paper. What used to be three separate models are now
// data-flag options of this one program:
//
//   link_type = 0 : interval cloglog hazard with association alpha
//                   (= the "joint longitudinal-interval" model)
//   link_type = 1 : threshold crossing with width sigma_thr
//                   (= the "latent first-hitting-time" model)
//   + N_dur > 0 and/or N_rel > 0 add durability / relapse transitions
//                   (= the "multi-state" extension)
//   use_spline    : patient-specific slope change after the knot
//
// Trajectory (l = log(1+t), t in YEARS):
//   m_i(l) = b0 + b1*l + b2*l^2 + u0_i + u1_i*l + zeta_i*(l-kappa_l)_+
// with b2 > 0 so the running minimum is the (clamped) vertex value and
// the running maximum is attained at a segment endpoint -- both closed
// form, so no grid approximation is needed.
//
// Sign conventions for the crossing device (C_i ~ N(c_D, sigma_thr^2)):
//   P(DMR attained by t)  = Phi((c_D - Mmin(t)) / sigma_thr)
//   P(not yet attained)   = Phi((Mmin(t) - c_D) / sigma_thr)
//   P(relapsed by t)      = Phi((Mmax(t) - c_D) / sigma_thr)
// Differences of these are formed with log_diff_exp for stability.
// =====================================================================
functions {
  // log(exp(la) - exp(lb)) with a finite floor.
  // The crossing likelihood is a difference of two CDFs evaluated at the
  // running minimum at L and at R. Because the running minimum is
  // non-increasing, Mmin(L) >= Mmin(R); but if the trajectory attains its
  // minimum BEFORE L, the two are numerically equal and the event
  // probability is exactly zero -> log_diff_exp returns -inf and the chain
  // cannot start or move. Flooring keeps the target finite (an event in
  // that window is merely astronomically unlikely, not impossible).
  real safe_log_diff(real la, real lb) {
    if (la <= lb) return -700;
    return log_diff_exp(la, lb);
  }

  // trajectory value at l
  real m_at(real c0, real a1, real b2, real zeta, real kap,
            real l, int spline) {
    real v = c0 + a1 * l + b2 * square(l);
    if (spline == 1 && l > kap) v += zeta * (l - kap);
    return v;
  }

  // running MINIMUM over [0, T]: candidates are segment endpoints, the
  // knot, and each segment's vertex clamped into that segment.
  real run_min(real c0, real a1, real b2, real zeta, real kap,
               real T, int spline) {
    real hi1 = (spline == 1) ? fmin(T, kap) : T;
    real m = m_at(c0, a1, b2, zeta, kap, 0.0, spline);
    real v1 = fmin(fmax(-a1 / (2 * b2), 0.0), hi1);          // vertex, seg 1
    m = fmin(m, m_at(c0, a1, b2, zeta, kap, hi1, spline));
    m = fmin(m, m_at(c0, a1, b2, zeta, kap, v1, spline));
    if (spline == 1 && T > kap) {
      real a1b = a1 + zeta;                                   // slope after knot
      real v2 = fmin(fmax(-a1b / (2 * b2), kap), T);          // vertex, seg 2
      m = fmin(m, m_at(c0, a1, b2, zeta, kap, kap, spline));
      m = fmin(m, m_at(c0, a1, b2, zeta, kap, T,   spline));
      m = fmin(m, m_at(c0, a1, b2, zeta, kap, v2,  spline));
    }
    return m;
  }

  // running MAXIMUM over [lo, hi]: each piece is convex (b2 > 0), so the
  // maximum sits at a segment endpoint (or the knot).
  real run_max(real c0, real a1, real b2, real zeta, real kap,
               real lo, real hi, int spline) {
    real m = m_at(c0, a1, b2, zeta, kap, lo, spline);
    m = fmax(m, m_at(c0, a1, b2, zeta, kap, hi, spline));
    if (spline == 1 && lo < kap && kap < hi)
      m = fmax(m, m_at(c0, a1, b2, zeta, kap, kap, spline));
    return m;
  }
}

data {
  // ---- switches ------------------------------------------------------
  int<lower=0, upper=1> link_type;    // 0 = interval hazard, 1 = crossing
  int<lower=0, upper=1> use_spline;
  real knot_years;                    // clinical knot (e.g. 2 years)
  real c_D;                           // DMR threshold, e.g. -4.5
  real floor_value;                   // left-censoring point, e.g. -5.0 or -4.5

  // ---- longitudinal --------------------------------------------------
  int<lower=1> N_obs;
  int<lower=1> N_pat;
  array[N_obs] int<lower=1, upper=N_pat> id_obs;
  vector[N_obs] y;
  array[N_obs] int<lower=0, upper=1> is_floor;
  vector<lower=0>[N_obs] t_obs;       // YEARS
  vector[N_obs] sample_bm;

  // ---- event: interval-hazard form (used when link_type == 0) --------
  int<lower=0> N_int;
  array[N_int] int<lower=1, upper=N_pat> id_int;
  vector<lower=0>[N_int] t_start;     // YEARS
  vector<lower=0>[N_int] t_end;       // YEARS
  vector<lower=0>[N_int] gap;         // YEARS
  array[N_int] int<lower=0, upper=1> event_interval;

  // ---- event: crossing form (used when link_type == 1) ---------------
  int<lower=0> N_ons;
  array[N_ons] int<lower=1, upper=N_pat> id_ons;
  vector<lower=0>[N_ons] ons_L;       // YEARS, 0 if left end is baseline
  vector<lower=0>[N_ons] ons_R;       // YEARS (censoring time if no event)
  array[N_ons] int<lower=0, upper=1> ons_event;

  // ---- multi-state: durability (0 rows disables) ---------------------
  int<lower=0> N_dur;
  array[N_dur] int<lower=1, upper=N_pat> id_dur;
  vector<lower=0>[N_dur] dur_lo;      // YEARS (onset time)
  vector<lower=0>[N_dur] dur_hi;      // YEARS (onset + window w)
  array[N_dur] int<lower=0, upper=1> dur_ok;   // 1 = stayed below c_D

  // ---- multi-state: relapse (0 rows disables) ------------------------
  int<lower=0> N_rel;
  array[N_rel] int<lower=1, upper=N_pat> id_rel;
  vector<lower=0>[N_rel] rel_L;       // YEARS (post-onset)
  vector<lower=0>[N_rel] rel_R;       // YEARS
  array[N_rel] int<lower=0, upper=1> rel_event;
}

transformed data {
  int n_haz = (link_type == 0) ? 1 : 0;
  int n_thr = (link_type == 1) ? 1 : 0;
  int n_spl = use_spline;
  int n_spl_pat = use_spline ? N_pat : 0;
  real kap_l = log1p(knot_years);     // knot on the log(1+t) scale
}

// NOTE ON UPPER BOUNDS: the positive parameters below carry generous
// upper bounds. Stan samples them on the log scale, so an extreme proposal
// during early warmup can exponentiate to Inf, making mu = c0 + a1*l + b2*l^2
// infinite and triggering "normal_lcdf: Location parameter is inf".
// All posterior mass sits far below these bounds (beta_time2 ~ 0.5,
// sigma_y ~ 1.8, tau_b ~ 1.0/2.2, sigma_thr ~ 0.2-0.6), so they prevent
// overflow without truncating the posterior. Check the summaries: if any
// parameter piles up against its bound, widen it rather than ignore it.
parameters {
  real beta0;
  real beta_time;
  real<lower=0, upper=20> beta_time2;   // convexity: closed-form extrema
  real beta_bm;
  real<lower=0, upper=20> sigma_y;

  vector<lower=0, upper=20>[2] tau_b;
  matrix[2, N_pat] z_b;

  // link-specific (size 0 when unused)
  array[n_haz] real gamma0;
  array[n_haz] real gamma_time;
  array[n_haz] real gamma_gap;
  array[n_haz] real alpha_mrd;
  array[n_thr] real<lower=0, upper=20> sigma_thr;

  // spline (size 0 when unused)
  array[n_spl] real<lower=0, upper=20> sigma_zeta;
  vector[n_spl_pat] z_zeta;
}

transformed parameters {
  matrix[N_pat, 2] b;
  vector[N_pat] zeta = rep_vector(0.0, N_pat);
  for (i in 1:N_pat) {
    b[i, 1] = tau_b[1] * z_b[1, i];
    b[i, 2] = tau_b[2] * z_b[2, i];
  }
  if (use_spline == 1)
    for (i in 1:N_pat) zeta[i] = sigma_zeta[1] * z_zeta[i];
}

model {
  // ---- priors --------------------------------------------------------
  beta0 ~ normal(-2.5, 2);
  beta_time ~ normal(-1, 1);
  beta_time2 ~ normal(0, 0.5);
  beta_bm ~ normal(0, 1);
  sigma_y ~ exponential(1);
  tau_b ~ exponential(1);
  to_vector(z_b) ~ normal(0, 1);

  if (link_type == 0) {
    gamma0[1] ~ normal(-2, 2);
    gamma_time[1] ~ normal(0, 1);
    gamma_gap[1] ~ normal(0, 1);
    alpha_mrd[1] ~ normal(-0.5, 0.75);
  } else {
    sigma_thr[1] ~ exponential(2);
  }
  if (use_spline == 1) {
    // regularised: the unconstrained fit gave an implausibly large
    // sigma_zeta (~6.4) on only a handful of relapse events.
    sigma_zeta[1] ~ normal(0, 1);
    z_zeta ~ normal(0, 1);
  }

  // ---- longitudinal sub-model (left-censored at the assay floor) -----
  for (n in 1:N_obs) {
    int i = id_obs[n];
    real l = log1p(t_obs[n]);
    real mu = m_at(beta0 + b[i, 1], beta_time + b[i, 2], beta_time2,
                   zeta[i], kap_l, l, use_spline)
              + beta_bm * sample_bm[n];
    if (is_floor[n] == 1) target += normal_lcdf(floor_value | mu, sigma_y);
    else                  y[n] ~ normal(mu, sigma_y);
  }

  // ---- event sub-model: interval cloglog hazard ----------------------
  if (link_type == 0) {
    for (m in 1:N_int) {
      int i = id_int[m];
      real lmid = log1p(0.5 * (t_start[m] + t_end[m]));
      real mu_mid = m_at(beta0 + b[i, 1], beta_time + b[i, 2], beta_time2,
                         zeta[i], kap_l, lmid, use_spline);
      real eta = gamma0[1] + gamma_time[1] * lmid
                 + gamma_gap[1] * log1p(gap[m]) + alpha_mrd[1] * mu_mid;
      real ch = exp(eta) * fmax(gap[m], 1e-6);
      if (event_interval[m] == 1) target += log1m_exp(-ch);
      else                        target += -ch;
    }
  }

  // ---- event sub-model: threshold crossing (onset) -------------------
  if (link_type == 1) {
    for (m in 1:N_ons) {
      int i = id_ons[m];
      real c0 = beta0 + b[i, 1];
      real a1 = beta_time + b[i, 2];
      real MR = run_min(c0, a1, beta_time2, zeta[i], kap_l,
                        log1p(ons_R[m]), use_spline);
      if (ons_event[m] == 1) {
        real ML = run_min(c0, a1, beta_time2, zeta[i], kap_l,
                          log1p(ons_L[m]), use_spline);
        // P(not by L) - P(not by R), both = Phi((Mmin - c_D)/sigma)
        target += safe_log_diff(normal_lcdf(ML | c_D, sigma_thr[1]),
                                normal_lcdf(MR | c_D, sigma_thr[1]));
      } else {
        target += normal_lcdf(MR | c_D, sigma_thr[1]);   // not yet attained
      }
    }
  }

  // ---- multi-state: durability (running max stays below c_D) ---------
  for (m in 1:N_dur) {
    int i = id_dur[m];
    real MX = run_max(beta0 + b[i, 1], beta_time + b[i, 2], beta_time2,
                      zeta[i], kap_l, log1p(dur_lo[m]), log1p(dur_hi[m]),
                      use_spline);
    if (dur_ok[m] == 1) target += normal_lccdf(MX | c_D, sigma_thr[1]);
    else                target += normal_lcdf(MX | c_D, sigma_thr[1]);
  }

  // ---- multi-state: relapse (first up-crossing after onset) ----------
  for (m in 1:N_rel) {
    int i = id_rel[m];
    real c0 = beta0 + b[i, 1];
    real a1 = beta_time + b[i, 2];
    real XR = run_max(c0, a1, beta_time2, zeta[i], kap_l,
                      log1p(rel_L[m]), log1p(rel_R[m]), use_spline);
    if (rel_event[m] == 1) {
      real XL = m_at(c0, a1, beta_time2, zeta[i], kap_l,
                     log1p(rel_L[m]), use_spline);
      target += safe_log_diff(normal_lcdf(XR | c_D, sigma_thr[1]),
                              normal_lcdf(XL | c_D, sigma_thr[1]));
    } else {
      target += normal_lccdf(XR | c_D, sigma_thr[1]);
    }
  }
}

generated quantities {
  // pointwise log-lik for the LONGITUDINAL block (comparable across all
  // variants; the event block differs in dimension between link types, so
  // LOO across link_type is not like-for-like -- see the fitting script).
  vector[N_obs] log_lik;
  // vertex reparameterisation: clinically interpretable nadir
  real nadir_l = -beta_time / (2 * beta_time2);
  real nadir_years = exp(nadir_l) - 1;
  real nadir_depth = beta0 + beta_time * nadir_l + beta_time2 * square(nadir_l);

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real l = log1p(t_obs[n]);
    real mu = m_at(beta0 + b[i, 1], beta_time + b[i, 2], beta_time2,
                   zeta[i], kap_l, l, use_spline) + beta_bm * sample_bm[n];
    log_lik[n] = (is_floor[n] == 1)
      ? normal_lcdf(floor_value | mu, sigma_y)
      : normal_lpdf(y[n] | mu, sigma_y);
  }
}
