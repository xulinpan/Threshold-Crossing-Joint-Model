# Longitudinal Time-Function Audit and Recommendation

## A. Recommended Primary Time Function

For a renewed Bayesian analysis, the preferred primary time function is a low-rank penalized spline for the population longitudinal trend, optionally combined with a soft monotone-decline regularization term on the population curve.

This recommendation is conditional on a full Bayesian refit with assay-floor censoring. The screening benchmark fitted in `04_Code/R/10_time_function_comparison.R` found that a 4-df spline proxy had the lowest AIC (1723.4), whereas the current log-quadratic model had the lowest non-floor RMSE (1.095) and lowest patient-level Brier score (0.179). Therefore, the spline should not be claimed as superior until the full joint model is refitted and posterior predictive checks and DMR calibration are recomputed.

If the manuscript is submitted without refitting the Bayesian joint model, the current log-quadratic time function remains the most defensible primary analysis because it is parsimonious, interpretable, and already has acceptable HMC diagnostics. In that case, the spline and piecewise-linear specifications should be reported as sensitivity analyses or a planned model extension.

## B. Recommended Sensitivity Time Functions

1. Current log-quadratic function:
   Use as a parsimony sensitivity analysis and as the benchmark against which any flexible curve must improve.

2. Clinically knotted piecewise-linear function:
   Use knots at 3, 6, 12, 24, and 60 months. This model is easy for clinicians to interpret because the knots match early molecular-response and longer-term monitoring windows. It carries more fixed-effect parameters than the log-quadratic model and should not include patient-level random slopes for every segment in this cohort.

3. Soft monotone-decline spline:
   Use as a biological sensitivity analysis. The constraint should be applied to the population mean trend, not to every individual trajectory, because individual patients can experience assay variation, nonadherence, resistance, or late molecular increase. A soft constraint is preferable to a hard constraint unless there is a strong clinical reason to exclude rebound.

## C. Mathematical Specification

Let \(t_{ij}^{(m)}\) be treatment time in months and \(x_{ij}=t_{ij}^{(m)}/12\) be treatment time in years. The likelihood uses years for numerical scaling, while tables and figures may report months for clinical readability. Let
\[
\ell(x)=\log(1+x).
\]

The general longitudinal model is
\[
\eta_i(x_{ij}) =
\beta_0 + f(x_{ij}) + \beta_{\mathrm{BM}} I(\mathrm{bone\ marrow}_{ij})
+ b_{0i} + b_{1i}\ell(x_{ij}),
\]
where \(f(x)\) is the population time trend. Independent random effects are recommended for the primary model:
\[
b_{0i}\sim N(0,\tau_0^2), \qquad b_{1i}\sim N(0,\tau_1^2).
\]

For non-floor observations,
\[
y_{ij}\mid \eta_i(x_{ij}),\sigma_y \sim N\{\eta_i(x_{ij}),\sigma_y^2\}.
\]
For assay-floor observations at \(c=-5.0\),
\[
\Pr(y_{ij}^{\mathrm{obs}}=c)
=
\Pr(y_{ij}^{\ast}\le c)
=
\Phi\left\{\frac{c-\eta_i(x_{ij})}{\sigma_y}\right\}.
\]

### Current Log-Quadratic Time Function

\[
f_{\log}(x)=\beta_1\ell(x)+\beta_2\ell(x)^2.
\]

This model is parsimonious and stable, but its shape is globally determined by two parameters. It can represent rapid early decline and flattening, but it may impose an unrealistic late curvature if the follow-up tail is sparse.

### Piecewise-Linear Time Function

Let \(\kappa=(0.25,0.5,1,2,5)\) years, corresponding to 3, 6, 12, 24, and 60 months. Then
\[
f_{\mathrm{PL}}(x)=\theta_1 x+\sum_{r=1}^{5}\theta_{r+1}(x-\kappa_r)_+,
\qquad
(z)_+=\max(z,0).
\]

The slope between adjacent knots is interpretable as the average molecular-response rate during the corresponding clinical period. In this cohort, the model should use a random intercept plus at most one patient-level random slope, not random slopes for all segments.

### Penalized Spline Time Function

Let \(B_1(x),\ldots,B_K(x)\) be centered spline basis functions with \(K\) kept small, for example \(K=4\) to \(6\). The time trend is
\[
f_{\mathrm{S}}(x)=\sum_{k=1}^{K}\theta_k B_k(x).
\]
Stability is obtained by a second-order random-walk penalty:
\[
\Delta^2\theta_k=\theta_k-2\theta_{k-1}+\theta_{k-2}
\sim N(0,\sigma_f^2), \qquad k=3,\ldots,K,
\]
with \(\sigma_f\) given a shrinkage prior.

### Soft Monotone-Decline Regularization

Define grid points \(0=x_1<\cdots<x_G\) covering the observed follow-up range. A soft population-level monotonicity penalty can be written as
\[
\log p_{\mathrm{mono}} \propto
-\frac{1}{2\sigma_{\mathrm{mono}}^2}
\sum_{g=2}^{G}\left[\max\{0, f(x_g)-f(x_{g-1})\}\right]^2.
\]
This discourages upward drift in the population curve while allowing late increases when supported by the data. A hard monotone alternative can be built with decreasing I-spline bases,
\[
f_{\mathrm{M}}(x)=\theta_0-\sum_{k=1}^{K} a_k I_k(x), \qquad a_k>0,
\]
but this is more restrictive and should be treated as sensitivity analysis in CML because true loss of molecular control is clinically possible.

The interval component can retain the current joint longitudinal--interval structure:
\[
p_{ik}=1-\exp\{-\Delta_{ik}\exp(\gamma_0+\gamma_1\ell(x_{ik}^{mid})
+\gamma_2\log(1+\Delta_{ik})+\alpha\eta_i(x_{ik}^{mid}))\}.
\]

## D. Stan Implementation Guidance

Build the spline basis in R and pass it to Stan rather than constructing splines inside Stan. Use the same centered basis for observed longitudinal visits and interval midpoints:

```stan
data {
  int<lower=1> N_obs;
  int<lower=1> N_int;
  int<lower=1> N_pat;
  int<lower=2> K_time;
  matrix[N_obs, K_time] B_obs;
  matrix[N_int, K_time] B_mid;
  vector[N_obs] log_time_obs;
  vector[N_int] log_time_mid;
  ...
}
parameters {
  real beta0;
  vector[K_time] theta_time;
  real beta_bm;
  real<lower=0> sigma_y;
  real<lower=0> sigma_f;
  vector<lower=0>[2] tau_b;
  matrix[2, N_pat] z_b;
  ...
}
model {
  beta0 ~ normal(-2.5, 2);
  theta_time[1] ~ normal(0, 1);
  theta_time[2] ~ normal(0, 1);
  for (k in 3:K_time) {
    theta_time[k] - 2 * theta_time[k - 1] + theta_time[k - 2] ~ normal(0, sigma_f);
  }
  sigma_f ~ exponential(2);
  tau_b ~ exponential(1);
  to_vector(z_b) ~ normal(0, 1);

  for (n in 1:N_obs) {
    int i = id_obs[n];
    real mu = beta0 + dot_product(B_obs[n], theta_time) +
      beta_bm * sample_bm[n] + tau_b[1] * z_b[1, i] +
      tau_b[2] * z_b[2, i] * log_time_obs[n];
    if (is_floor[n] == 1) {
      target += normal_lcdf(floor_value | mu, sigma_y);
    } else {
      y[n] ~ normal(mu, sigma_y);
    }
  }
}
```

For the soft monotone sensitivity, pass a grid basis matrix `B_grid` and add a weak penalty on positive increments of the population curve. A smooth positive-part function is preferable to a hard `fmax` for HMC:

```stan
functions {
  real softplus_scaled(real x, real eps) {
    return eps * log1p_exp(x / eps);
  }
}
...
for (g in 2:G_grid) {
  real delta = dot_product(B_grid[g], theta_time) -
    dot_product(B_grid[g - 1], theta_time);
  target += normal_lpdf(softplus_scaled(delta, 0.05) | 0, sigma_mono);
}
```

Recommended priors for the renewed spline model:

- \(\beta_0\sim N(-2.5,2^2)\)
- \(\beta_{\mathrm{BM}}\sim N(0,1^2)\)
- \(\sigma_y\sim \mathrm{Exponential}(1)\)
- \(\tau_0,\tau_1\sim \mathrm{Exponential}(1)\)
- \(\sigma_f\sim \mathrm{Exponential}(2)\) or half-normal \((0,0.5)\)
- \(\sigma_{\mathrm{mono}}\sim \mathrm{Exponential}(5)\) for a moderate soft monotone penalty
- \(\gamma_0\sim N(-2,2^2)\), \(\gamma_1,\gamma_2\sim N(0,1^2)\), and \(\alpha\sim N(-0.5,0.75^2)\)

Validation after refitting should include R-hat, bulk and tail ESS, divergences, maximum treedepth, E-BFMI, posterior predictive checks for non-floor observations, posterior predicted floor probability, interval-level calibration, patient-level calibration, and comparison with the current log-quadratic model.

## E. Interpretation for the Current Cohort

The screening comparison does not prove that the spline joint model is clinically superior. It shows that added flexibility may improve longitudinal fit by AIC, but calibration proxies are not uniformly better. In a cohort of 87 patients with 48.3% assay-floor observations, a flexible time trend should be used only with shrinkage and with the same floor-censoring likelihood as the primary model. The most statistically defensible statement is that the low-rank penalized spline is the preferred renewed model candidate, while the current log-quadratic model remains the reference primary model until the full Bayesian refit is completed.
