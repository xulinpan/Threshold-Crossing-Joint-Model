# Latent Threshold-Crossing Model for DMR

Prepared: 2026-07-09

Manuscript: *A Bayesian Joint Longitudinal--Interval Model for Deep Molecular Response Under Irregular Molecular Monitoring in Chronic Myeloid Leukemia*

## Executive Recommendation

The most coherent redesign is a latent threshold-crossing model. In this model, DMR is not a separate event process predicted by latent MRD. Instead, DMR onset is the first time the latent biological MRD trajectory crosses the DMR threshold:

```text
T_i^DMR = inf{t >= 0 : eta_i(t) <= c_DMR},  c_DMR = -4.5.
```

This directly addresses the endpoint-circularity concern in the current interval-hazard model, where latent log-MRD is used as a predictor of an endpoint that is itself defined by log-MRD <= -4.5.

Important implementation point: the assay-source effect should be treated as an observation-level measurement shift, not as part of the biological trajectory that defines DMR onset. Therefore I recommend the notation below:

- `m_i(t)`: latent true biological log-MRD trajectory.
- `mu_ij = m_i(t_ij) + beta_source I(BM_ij)`: expected observed log-MRD on the assay/source scale.
- DMR crossing uses `m_i(t)`, not `mu_ij`.

If the manuscript keeps the symbol `eta_i(t)`, define it as the biological trajectory and put source adjustment in the observation model.

## A. Full Mathematical Model

### Time Scale

Let time be measured in years since treatment initiation for all model components. Descriptive tables may continue to report months.

### Latent Biological Trajectory

For patient `i` at time `t`,

```latex
m_i(t) =
\beta_0 + f(t) + b_{0i} + b_{1i} g(t),
```

where a practical choice matching the current code is

```latex
f(t) = \beta_1 \log(1+t) + \beta_2 \{\log(1+t)\}^2,
\qquad
g(t) = \log(1+t).
```

Independent patient effects are recommended for the small cohort:

```latex
b_{0i} \sim N(0,\tau_{b0}^2), \qquad
b_{1i} \sim N(0,\tau_{b1}^2).
```

The observation-scale mean is

```latex
\mu_{ij} = m_i(t_{ij}) + \beta_{\mathrm{BM}} I(\mathrm{bone\ marrow}_{ij}).
```

### Longitudinal Observation Model

For an exact non-floor log-MRD observation,

```latex
y_{ij} \mid m_i(t_{ij}) \sim
N(\mu_{ij}, \sigma_y^2).
```

For a floor-limited observation reported as log-MRD <= -5.0,

```latex
P(y_{ij} \leq -5.0 \mid m_i(t_{ij}))
=
\Phi\left\{ \frac{-5.0-\mu_{ij}}{\sigma_y} \right\}.
```

For optional interval-reported molecular observations, for example `a < y <= b`,

```latex
P(a < y_{ij} \leq b \mid m_i(t_{ij}))
=
\Phi\left\{ \frac{b-\mu_{ij}}{\sigma_y} \right\}
-
\Phi\left\{ \frac{a-\mu_{ij}}{\sigma_y} \right\}.
```

### DMR Onset as Latent Threshold Crossing

Let

```latex
c_{\mathrm{DMR}} = -4.5.
```

The latent DMR onset time is

```latex
T_i^{\mathrm{DMR}}
=
\inf\{t \geq 0: m_i(t) \leq c_{\mathrm{DMR}}\}.
```

For irregular visits, define the observed DMR interval as `(L_i, R_i]`, where `L_i` is the last visit before documented DMR and `R_i` is the first visit with documented DMR. For a censored patient, let `C_i` be the last observed follow-up time without DMR.

For a general nonmonotone trajectory, define the running minimum

```latex
M_i(t) = \min_{0 \leq u \leq t} m_i(u).
```

Then crossing by time `t` occurs when `M_i(t) <= c_DMR`.

## Likelihood Contributions

### 1. Patients With Observed First DMR

For deterministic threshold crossing,

```latex
\mathcal{L}_{i,\mathrm{DMR}}
=
I\{M_i(L_i) > c_{\mathrm{DMR}},\; M_i(R_i) \leq c_{\mathrm{DMR}}\}.
```

If `m_i(t)` is constrained or assumed monotone decreasing, this reduces to

```latex
\mathcal{L}_{i,\mathrm{DMR}}
=
I\{m_i(L_i) > c_{\mathrm{DMR}},\; m_i(R_i) \leq c_{\mathrm{DMR}}\}.
```

This hard deterministic likelihood is mathematically coherent but not ideal for Hamiltonian Monte Carlo because it creates discontinuous target-density boundaries.

### 2. Patients Censored Without Observed DMR

For deterministic threshold crossing,

```latex
\mathcal{L}_{i,\mathrm{cens}}
=
I\{M_i(C_i) > c_{\mathrm{DMR}}\}.
```

Under monotonicity, this becomes `I{m_i(C_i) > c_DMR}`.

### 3. Patients With Floor-Limited Molecular Observations

Floor-limited observations contribute to the longitudinal likelihood, not to a separate event hazard:

```latex
\mathcal{L}_{ij,\mathrm{floor}}
=
\Phi\left\{ \frac{-5.0-\mu_{ij}}{\sigma_y} \right\}.
```

Because the floor threshold -5.0 is below the DMR threshold -4.5, a floor observation strongly informs the crossing model by making `m_i(t_ij) <= -4.5` plausible, but it does not need to be treated as an exact value.

## B. Two Threshold-Crossing Versions

### A. Deterministic Threshold Crossing

Definition:

```latex
T_i^{\mathrm{DMR}} = \inf\{t: M_i(t) \leq c_{\mathrm{DMR}}\}.
```

Event interval:

```latex
I\{M_i(L_i)>c_{\mathrm{DMR}}\}
I\{M_i(R_i)\leq c_{\mathrm{DMR}}\}.
```

Censoring:

```latex
I\{M_i(C_i)>c_{\mathrm{DMR}}\}.
```

Stan recommendation: avoid literal hard indicators in the primary Stan implementation. Use a near-deterministic smooth approximation:

```latex
\log \mathcal{L}_{i,\mathrm{DMR}}
\approx
\logit^{-1}_{\log}\left\{\frac{M_i(L_i)-c_{\mathrm{DMR}}}{\epsilon}\right\}
+
\logit^{-1}_{\log}\left\{\frac{c_{\mathrm{DMR}}-M_i(R_i)}{\epsilon}\right\},
```

where `logit^{-1}_{log}(x) = log(inv_logit(x))`, and `epsilon` is a fixed small scale, such as 0.03 to 0.05 log10 units.

Censoring:

```latex
\log \mathcal{L}_{i,\mathrm{cens}}
\approx
\logit^{-1}_{\log}\left\{\frac{M_i(C_i)-c_{\mathrm{DMR}}}{\epsilon}\right\}.
```

This is a computational approximation to deterministic crossing, not a new clinical assumption.

### B. Probabilistic Threshold Crossing With Threshold Uncertainty

Let the individual or classification threshold be

```latex
C_i \sim N(c_{\mathrm{DMR}}, \sigma_{\mathrm{threshold}}^2).
```

Given the trajectory, crossing by time `t` occurs if `C_i >= M_i(t)`. Integrating over `C_i` gives a smooth likelihood.

For an observed event interval `(L_i, R_i]`:

```latex
\mathcal{L}_{i,\mathrm{DMR}}
=
P\{M_i(R_i) \leq C_i < M_i(L_i)\}
=
\Phi\left\{\frac{M_i(L_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{threshold}}}\right\}
-
\Phi\left\{\frac{M_i(R_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{threshold}}}\right\}.
```

For a censored patient:

```latex
\mathcal{L}_{i,\mathrm{cens}}
=
P\{C_i < M_i(C_i)\}
=
\Phi\left\{\frac{M_i(C_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{threshold}}}\right\}.
```

This version is smoother and more Stan-friendly. However, `sigma_threshold` can be weakly identified in a small cohort, so it should be regularized and reported as a sensitivity model unless the manuscript has strong assay or clinical justification for threshold uncertainty.

## B. Stan-Ready Likelihood Structure

The accompanying Stan skeleton is:

```text
04_Code/Stan/glw_latent_threshold_crossing_skeleton.stan
```

Key Stan implementation choices:

1. Use `m_i(t)` for biological crossing.
2. Use `m_i(t) + beta_source_bm * sample_bm` for observed log-MRD.
3. For exact observations:

```stan
y[n] ~ normal(mu_obs, sigma_y);
```

4. For floor observations:

```stan
target += normal_lcdf(floor_value | mu_obs, sigma_y);
```

5. For interval-censored molecular values:

```stan
target += log_diff_exp(
  normal_lcdf(y_upper[n] | mu_obs, sigma_y),
  normal_lcdf(y_lower[n] | mu_obs, sigma_y)
);
```

6. Approximate `M_i(t)` on a patient-specific time grid containing baseline, observed visits, event endpoints, censoring times, and optional intermediate points. For smoother HMC behavior use a soft minimum:

```stan
real soft_min_segment(vector x, real temperature) {
  return -temperature * log_sum_exp(-x / temperature);
}
```

7. Deterministic version uses smooth indicators:

```stan
target += log_inv_logit((min_to_left - threshold) / deterministic_eps) +
          log_inv_logit((threshold - min_to_right) / deterministic_eps);
```

8. Probabilistic threshold version uses:

```stan
target += log_diff_exp(
  normal_lcdf(min_to_left | threshold, sigma_threshold),
  normal_lcdf(min_to_right | threshold, sigma_threshold)
);
```

9. Censored patients use:

```stan
target += normal_lcdf(min_to_right | threshold, sigma_threshold);
```

where `min_to_right` is the running minimum through the last follow-up time.

## C. Recommended Priors

For 87 patients, use weakly informative but regularizing priors:

```latex
\beta_0 \sim N(-2.5, 2^2),
\qquad
\beta_1 \sim N(-1, 1^2),
\qquad
\beta_2 \sim N(0, 0.5^2),
```

```latex
\beta_{\mathrm{BM}} \sim N(0, 0.75^2),
\qquad
\sigma_y \sim \mathrm{Exponential}(1).
```

For random effects:

```latex
\tau_{b0}, \tau_{b1} \sim \mathrm{Exponential}(1),
\qquad
z_{0i}, z_{1i} \sim N(0,1),
```

with

```latex
b_{0i}=\tau_{b0} z_{0i}, \qquad b_{1i}=\tau_{b1} z_{1i}.
```

For deterministic crossing:

```latex
\epsilon = 0.03 \text{ to } 0.05
```

fixed as a computational smoothing parameter, not estimated.

For probabilistic threshold crossing:

```latex
\sigma_{\mathrm{threshold}} \sim \mathrm{half}\text{-}N(0,0.20^2)
```

or run fixed sensitivity values:

```text
sigma_threshold = 0.05, 0.10, 0.20 log10 units.
```

Recommendation: treat the deterministic or near-deterministic model as the primary conceptual model, and the probabilistic threshold model as a sensitivity model unless threshold uncertainty is scientifically justified by assay data.

## D. Why This Avoids Circularity Better Than the Interval Hazard Model

The current interval-hazard model contains a conceptual loop:

```text
latent log-MRD -> predicts first DMR,
but DMR is defined by log-MRD <= -4.5.
```

This is not automatically invalid, but it can look circular to reviewers because the event process is modeled as if it were distinct from the biomarker threshold.

The threshold-crossing model avoids this issue by making the endpoint a deterministic or probabilistic functional of the latent MRD trajectory:

```text
latent trajectory -> threshold crossing time -> observed interval.
```

There is no separate regression coefficient claiming that latent MRD "causes" DMR. Instead, DMR is what happens when the latent trajectory crosses the clinical threshold. The model therefore aligns the statistical endpoint with the clinical endpoint definition.

Important caveat: if the DMR interval is derived from the same exact observed log-MRD values already included in the longitudinal likelihood, the interval should not be treated as an independent additional data source without careful justification. The cleanest approach is to view the interval as coarsened endpoint information, or to use the model to derive posterior crossing-time summaries from the longitudinal likelihood. If exact values and derived threshold labels are both used, the manuscript should explain that the threshold likelihood is a coarsening constraint and not an independent event process.

## E. Simulation Plan

### Objectives

Evaluate whether the latent threshold-crossing model recovers:

1. patient-specific latent DMR onset time;
2. interval-censored DMR status under irregular monitoring;
3. population trajectory and random-effect parameters;
4. calibration of predicted crossing probabilities;
5. robustness to assay-floor censoring.

### Data-Generating Mechanism

1. Set `N=87` patients.
2. Generate visit schedules by resampling empirical visit times and gaps from the real cohort, preserving irregularity.
3. Generate latent trajectories:

```latex
m_i(t)=\beta_0+\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2+b_{0i}+b_{1i}\log(1+t).
```

4. Generate true crossing time:

```latex
T_i^{\mathrm{DMR}}=\inf\{t:m_i(t)\leq -4.5\}.
```

Compute this on a fine grid, such as every 0.01 years.

5. Generate observed log-MRD:

```latex
y_{ij}=m_i(t_{ij})+\beta_{\mathrm{BM}} I(\mathrm{BM}_{ij})+\epsilon_{ij},
\qquad
\epsilon_{ij}\sim N(0,\sigma_y^2).
```

6. Apply assay floor:

```text
if y_ij <= -5.0, record as floor/left-censored at -5.0.
```

7. Derive observed DMR interval from the visit sequence:

```text
L_i = last visit before first documented DMR,
R_i = first visit with documented DMR.
```

Use both latent-crossing-derived and observed-measurement-derived endpoint variants to quantify misclassification from measurement error.

### Simulation Scenarios

Use at least these scenarios:

1. regular 3-month monitoring versus empirical irregular monitoring;
2. low, moderate, and high floor frequency, approximately 20%, 50%, 70%;
3. measurement error `sigma_y = 0.5, 1.0, 1.5`;
4. deterministic threshold and probabilistic threshold with `sigma_threshold = 0.05, 0.10, 0.20`;
5. monotone trajectories and mildly nonmonotone trajectories;
6. complete follow-up and administrative censoring.

### Fitted Models to Compare

1. deterministic threshold-crossing model;
2. probabilistic threshold-crossing model;
3. current interval-hazard model;
4. exact-floor sensitivity model.

### Performance Metrics

For each scenario and replicate:

1. bias and RMSE of posterior median `T_i^DMR`;
2. coverage of 90% and 95% posterior intervals for `T_i^DMR`;
3. accuracy of patient-level ever-DMR probability;
4. interval-level Brier score;
5. patient-level Brier score;
6. calibration intercept and slope;
7. bias in fixed effects and random-effect scales;
8. convergence diagnostics: R-hat, ESS, divergences, E-BFMI.

### Recommended Number of Replicates

For a manuscript supplement:

```text
100 simulation replicates per main scenario
```

is a reasonable starting target. For a lighter revision:

```text
25 to 50 replicates
```

may be acceptable if clearly labeled as an illustrative operating-characteristic check.

## F. Manuscript-Ready Methods Text

```latex
\subsection{Latent Threshold-Crossing Model for DMR}

As a sensitivity and conceptual refinement, we formulated DMR onset as a latent threshold-crossing time rather than as a separate interval hazard predicted by latent MRD. Let \(m_i(t)\) denote the latent true biological log-MRD trajectory for patient \(i\) at treatment time \(t\), measured in years. The trajectory was modeled as
\[
m_i(t)=\beta_0+\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2
+b_{0i}+b_{1i}\log(1+t),
\]
with independent patient-specific random intercept and slope terms. The expected observed log-MRD value additionally included a sample-source adjustment,
\[
\mu_{ij}=m_i(t_{ij})+\beta_{\mathrm{BM}}I(\mathrm{bone\ marrow}_{ij}).
\]
Thus, sample source was treated as an observation-level measurement effect and did not alter the biological trajectory used to define DMR onset.

For exact non-floor observations, \(y_{ij}\sim N(\mu_{ij},\sigma_y^2)\). Observations reported at the assay floor contributed the left-censored likelihood
\[
P(y_{ij}\leq -5.0)=
\Phi\{(-5.0-\mu_{ij})/\sigma_y\}.
\]
If molecular values were available only as intervals, their likelihood contribution was the corresponding normal probability over the reporting interval.

DMR was defined by the clinical threshold \(c_{\mathrm{DMR}}=-4.5\). The latent DMR onset time was
\[
T_i^{\mathrm{DMR}}=\inf\{t\geq 0:m_i(t)\leq c_{\mathrm{DMR}}\}.
\]
For irregularly monitored patients with observed first DMR, the event was represented as \(L_i<T_i^{\mathrm{DMR}}\leq R_i\), where \(L_i\) was the last visit before documented DMR and \(R_i\) was the first visit with documented DMR. Patients without documented DMR were treated as censored after their last observed follow-up time. For a general trajectory, the event likelihood was expressed through the running minimum \(M_i(t)=\min_{0\leq u\leq t}m_i(u)\). An observed event interval contributed the condition \(M_i(L_i)>c_{\mathrm{DMR}}\) and \(M_i(R_i)\leq c_{\mathrm{DMR}}\), whereas a censored patient contributed \(M_i(C_i)>c_{\mathrm{DMR}}\).

We considered two threshold-crossing versions. The deterministic version used a near-deterministic smooth approximation to the threshold indicators for stable Hamiltonian Monte Carlo computation. The probabilistic version introduced threshold uncertainty \(C_i\sim N(c_{\mathrm{DMR}},\sigma_{\mathrm{threshold}}^2)\), giving
\[
P(L_i<T_i^{\mathrm{DMR}}\leq R_i)
=
\Phi\{(M_i(L_i)-c_{\mathrm{DMR}})/\sigma_{\mathrm{threshold}}\}
-
\Phi\{(M_i(R_i)-c_{\mathrm{DMR}})/\sigma_{\mathrm{threshold}}\}.
\]
For censored patients, the corresponding probability was
\[
\Phi\{(M_i(C_i)-c_{\mathrm{DMR}})/\sigma_{\mathrm{threshold}}\}.
\]
Weakly informative priors were used for all fixed effects, variance components, and threshold-uncertainty parameters. Because the model was developed in a modest single-cohort dataset, all threshold-crossing results were interpreted as monitoring-oriented model development rather than as externally validated clinical decision rules.
```

## G. Manuscript-Ready Discussion Paragraph

```latex
The latent threshold-crossing formulation addresses a conceptual limitation of modeling DMR as a separate event process predicted by latent MRD. Because DMR is clinically defined by log-MRD \(\leq -4.5\), a model in which latent MRD predicts first DMR can appear circular: the biomarker is used to predict an endpoint defined by the same biomarker. The threshold-crossing model resolves this by defining DMR onset as a functional of the latent biological MRD trajectory itself. In this framework, irregular monitoring produces interval-censored observation of the latent crossing time, and assay-floor values contribute through the longitudinal measurement model rather than as exact deep-response values. This provides a closer alignment between the statistical model and the clinical endpoint definition. The approach should still be interpreted as a development-stage monitoring model, because threshold-crossing probabilities and patient-specific crossing times require external validation and calibration before they can be used for treatment decisions or treatment-free-remission counseling.
```

## Practical Recommendation for This Manuscript

Use this model in one of two ways:

1. As a conceptual replacement for the interval-hazard joint model in a future full reanalysis.
2. As a strong sensitivity/model-redesign section in the current manuscript, explaining that it resolves the endpoint-circularity concern and should be the next primary model if reviewers object to the interval-hazard formulation.

Do not claim that this model is ready for treatment decision-making without external validation.
