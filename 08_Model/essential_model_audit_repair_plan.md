# Essential Model Audit and Repair Plan

Prepared: 2026-07-09

Manuscript: A Bayesian Joint Longitudinal--Interval Model for Deep Molecular Response Under Irregular Molecular Monitoring in Chronic Myeloid Leukemia

## Executive Summary

The current analysis is broadly coherent as a monitoring-oriented Bayesian development model, but it needs clearer endpoint auditing, a standardized time scale, stronger floor-value language, and a more conservative primary-model claim. The most important finding from this audit is that the identical DMR and CMR counts are not caused by a mismatch between stored endpoint flags and the stated thresholds. In the processed analysis data, DMR and CMR are identical because there are no retained observations with log-MRD in the DMR-only range `(-5.0, -4.5]`. All 239 retained DMR-positive observations are exact `-5.0` floor values.

One raw record should be adjudicated before final submission: raw zero-based row 500, spreadsheet row 502 including the header, has `IS = 0.002`, missing ABL, missing sample type, missing laboratory date, and missing `New LOG-MRD`. If log-MRD were derived from `log10(IS / 100)`, this record would be approximately -4.699, which is DMR-positive but CMR-negative. It was excluded from the processed analysis because essential modeling fields were incomplete. This is not a current coding error, but it is a data-adjudication item.

The fitted Stan model uses years, not months, for `t_obs`, `t_start`, `t_end`, and `gap`; these values are created by dividing processed month variables by 12. Descriptive cohort summaries may remain in months, but the model notation and coefficient interpretation must be stated in years.

The left-censoring likelihood for assay-floor observations is correctly implemented as `normal_lcdf(floor_value | mu, sigma_y)` for observations with `log_mrd <= -5.0`. This is appropriate if `-5.0` represents an undetectable/floor report. It should not be described as an exact measurement. Interval-censoring would require a finite lower and upper measurement interval; here the current information is `y <= -5.0`, so the model is left-censored.

The main overparameterization concern is not the random slope itself. The data include repeated monitoring for all patients, 80 of 87 patients have at least 3 visits, and 52 patients have at least 24 months of observed time range. The weak component is the random-effect correlation: `Omega_b[1,2]` is near zero, has R-hat 1.019, and has bulk ESS 191. The primary model should either use independent random intercept and random slope or retain the correlated structure but explicitly treat the correlation as a nuisance parameter. If the model is simplified, all posterior estimates and figures must be regenerated.

## A. Detailed Critique of the Current Model

### Endpoint Definitions

The endpoint definitions are implemented correctly in the processed longitudinal file:

- DMR: `log_mrd <= -4.5`.
- CMR: `log_mrd <= -5.0`.
- Stored DMR mismatches versus recomputation: 0.
- Stored CMR mismatches versus recomputation: 0.

The identical DMR and CMR counts are a data feature of the retained analysis set:

| Quantity | Value |
|---|---:|
| Retained longitudinal observations | 495 |
| Stored DMR-positive observations | 239 |
| Stored CMR-positive observations | 239 |
| Recomputed `log_mrd <= -4.5` | 239 |
| Recomputed `log_mrd <= -5.0` | 239 |
| Observations in `(-5.0, -4.5]` | 0 |
| Patients ever DMR | 68 |
| Patients ever CMR | 68 |
| Patients ever DMR but never CMR | 0 |

Across follow-up windows:

| Window, months | Observations | DMR | CMR | DMR-only |
|---|---:|---:|---:|---:|
| 0-3 | 67 | 2 | 2 | 0 |
| >3-6 | 36 | 8 | 8 | 0 |
| >6-12 | 79 | 26 | 26 | 0 |
| >12-24 | 95 | 51 | 51 | 0 |
| >24-60 | 135 | 90 | 90 | 0 |
| >60 | 83 | 62 | 62 | 0 |

Raw-data check:

- Raw `New LOG-MRD` nonmissing records: 495.
- Raw `New LOG-MRD <= -4.5`: 239.
- Raw `New LOG-MRD <= -5.0`: 239.
- Raw `New LOG-MRD` in `(-5.0, -4.5]`: 0.
- Raw exact `New LOG-MRD == -5.0`: 239.
- Among raw `New LOG-MRD == -5.0`, BCR-ABL copy number was zero in 239 rows and IS was zero in 239 rows.

Interpretation: the identical DMR/CMR count is mainly a floor/reporting feature, not an endpoint coding error. However, the excluded raw record with `IS = 0.002` should be reviewed because it could become a DMR-only record if clinically valid and if enough essential information can be recovered.

### Time Scale

The manuscript currently mixes months and years. The processed data store time in months:

- `t_months`.
- `gap_months`.
- interval `t_start`, `t_end`, and `gap_months`.

The Stan data-preparation script converts these to years:

- `t_obs = t_months / 12`.
- `t_start = interval$t_start / 12`.
- `t_end = interval$t_end / 12`.
- `gap = gap_months / 12`.

Therefore, the model uses years. The notation should define `u_{ij}` as years from imatinib initiation and should state that descriptive summaries are reported in months for clinical readability.

The coefficient `beta_time` is not a per-month effect and not a simple per-year linear effect. It is the effect of `log(1 + years)` in a nonlinear trajectory with a quadratic `log(1 + years)` term and patient-specific random slope on the same transformed scale.

### Assay-Floor Handling

The current Stan likelihood for floor observations is:

```stan
if (is_floor[n] == 1) {
  target += normal_lcdf(floor_value | mu, sigma_y);
} else {
  y[n] ~ normal(mu, sigma_y);
}
```

This is mathematically correct for left-censored observations with observed information `y <= -5.0`.

Recommended wording:

- Use "left-censored at the assay/reporting floor of -5.0".
- Avoid saying floor observations were "equal to -5.0" in the likelihood.
- Describe exact `-5.0` values as "reported floor values" or "floor-coded values".

When to use other treatments:

- Exact treatment is not preferred for the primary model because 239 observations are floor-coded and treating them as exact understates uncertainty below the floor.
- Left-censoring is appropriate when the only measurement information is `log-MRD <= -5.0`.
- Interval-censoring would be appropriate only if the laboratory report provided a finite interval, such as `a < log-MRD <= -5.0`. The current dataset does not provide such a lower finite bound.

### Model Complexity

The sample size is modest but not too small for a random-intercept/random-slope longitudinal component:

- 87 patients.
- 495 longitudinal observations.
- Median 5 visits per patient.
- 80 patients have at least 3 visits.
- 52 patients have at least 24 months of observed time range.
- 275 at-risk intervals.
- 68 first-DMR events.

Random intercept and random slope are defensible because the data show repeated longitudinal trajectories and posterior standard deviations are nontrivial:

- `tau_b[1] = 0.964`, 95% posterior interval 0.332 to 1.521.
- `tau_b[2] = 2.212`, 95% posterior interval 1.550 to 2.963.

The random-effect correlation is not well supported:

- `Omega_b[1,2] = 0.026`, 95% posterior interval -0.417 to 0.621.
- R-hat 1.019.
- Bulk ESS 191.

Recommendation:

1. Primary conservative option: keep random intercept and random slope but remove the random-effect correlation by fitting independent random effects. This requires reanalysis.
2. Acceptable current-reporting option: retain the current correlated model but explicitly state that the correlation is weakly estimated and not interpreted clinically.
3. Sensitivity option: compare correlated random intercept/slope, independent random intercept/slope, and random-intercept-only models using convergence, PPC, calibration, and predictive performance. This requires reanalysis.

### Interval DMR Component

The interval component is coherent for first observed DMR:

\[
P(d_{ik}=1 \mid \mathrm{at\ risk}) = 1 - \exp\{-h_{ik}\Delta_{ik}\}.
\]

However, interpretation of `gamma_gap` must be careful because `gap` appears both in the linear predictor through `log(1 + gap)` and in the cumulative hazard multiplier. A negative `gamma_gap` does not necessarily mean longer gaps reduce cumulative DMR probability in an intuitive clinical sense; it means the hazard-rate component decreases with gap after the interval-length multiplier is also included. This should be described as a nuisance adjustment for monitoring interval length, not as a causal monitoring-frequency effect.

## B. Corrected Model Specification

Let \(u_{ij}\) be years from imatinib initiation for patient \(i\) at visit \(j\), and let \(y_{ij}\) be the observed log-MRD value. The processed data are stored in months, but all fitted model time variables are divided by 12 and analyzed in years.

Longitudinal latent trajectory:

\[
\mu_{ij} =
\beta_0 + \beta_1 \log(1+u_{ij})
+ \beta_2\{\log(1+u_{ij})\}^2
+ \beta_3 I(\mathrm{bone\ marrow}_{ij})
+ b_{0i} + b_{1i}\log(1+u_{ij}).
\]

For non-floor observations:

\[
y_{ij} \sim N(\mu_{ij}, \sigma_y^2).
\]

For floor-coded observations at \(c=-5.0\):

\[
P(y_{ij} \le c \mid \mu_{ij}, \sigma_y) =
\Phi\left(\frac{c-\mu_{ij}}{\sigma_y}\right).
\]

For interval \(k\) of patient \(i\), define interval start \(L_{ik}\), end \(R_{ik}\), and gap \(\Delta_{ik}=R_{ik}-L_{ik}\), all in years. Let \(u_{ik}^{mid}=(L_{ik}+R_{ik})/2\). The interval hazard-rate component is

\[
h_{ik} =
\exp\{\gamma_0 + \gamma_1\log(1+u_{ik}^{mid})
+\gamma_2\log(1+\Delta_{ik})
+\alpha\mu_i(u_{ik}^{mid})\}.
\]

The first-DMR interval probability is

\[
P(d_{ik}=1 \mid \mathrm{at\ risk}) =
1-\exp(-h_{ik}\Delta_{ik}).
\]

Current priors used in Stan:

| Parameter | Prior |
|---|---|
| \(\beta_0\) | Normal(-2.5, 2) |
| \(\beta_1\) | Normal(-1, 1) |
| \(\beta_2\) | Normal(0, 0.5) |
| \(\beta_3\) | Normal(0, 1) |
| \(\sigma_y\) | Exponential(1) |
| \(\gamma_0\) | Normal(-2, 2) |
| \(\gamma_1\) | Normal(0, 1) |
| \(\gamma_2\) | Normal(0, 1) |
| \(\alpha\) | Normal(-0.5, 0.75) |
| \(\tau_{b0}, \tau_{b1}\) | Exponential(1) |
| \(L_b\) | LKJ Cholesky correlation prior with shape 2 |
| \(z_b\) | Standard normal |

Recommended primary-model repair:

- Keep the random intercept and random slope.
- Prefer independent random effects for primary reporting unless the correlated model is retained as the already-fitted model with the correlation declared nuisance.
- If independent random effects are adopted, re-fit the model and regenerate posterior summaries, PPC, calibration, and figures. Current estimates cannot be reused after this structural change.

## C. Required Coding Checks

Minimum required checks before submission:

1. Recompute `dmr = as.integer(log_mrd <= -4.5)` and `cmr = as.integer(log_mrd <= -5.0)` from the final analysis CSV and assert zero mismatches.
2. Count observations and patients in `(-5.0, -4.5]`; report if zero rather than leaving readers to infer why DMR and CMR counts are identical.
3. Audit raw records with nonmissing IS but missing `New LOG-MRD`, especially raw row 500 / spreadsheet row 502 with `IS = 0.002`.
4. Decide whether any excluded row can be clinically recovered. If recovered, rebuild longitudinal, interval, and patient-level datasets, then re-fit all models. This requires reanalysis.
5. Confirm that floor-coded retained observations are exactly the rows with `log_mrd == -5.0`, BCR-ABL copy number 0, and IS 0.
6. Confirm no retained positive-IS row has derived `log10(IS / 100) <= -4.5` while `log_mrd > -4.5` or missing.
7. Assert that Stan model time variables are in years and stored descriptive variables are in months.
8. Print ranges and summaries of `t_obs`, `gap`, `t_start`, and `t_end` after conversion to years.
9. Add a priors table to the manuscript or supplement.
10. Add a diagnostic table with post-warmup draws, divergences, maximum treedepth, minimum E-BFMI, maximum R-hat, and low-ESS parameters.
11. Report random-effect correlation diagnostics separately and avoid clinical interpretation.
12. Run posterior predictive checks for non-floor log-MRD and floor probability.
13. Run interval-level and patient-level calibration summaries.
14. If simplifying the random-effects correlation, re-fit and regenerate every posterior and predictive output.
15. Add a reduced-model sensitivity analysis: independent random intercept/slope and random-intercept-only. This requires reanalysis.

## D. Revised Manuscript-Ready Methods Text

### Data Structure and Endpoints

The raw molecular monitoring file was processed into privacy-preserving analysis datasets using coded patient identifiers. Treatment follow-up was stored descriptively in months from imatinib initiation. DMR was defined as log-MRD \(\leq -4.5\), and CMR was defined as log-MRD \(\leq -5.0\). These endpoint indicators were recomputed directly from the final log-MRD values as part of the analysis audit. First observed DMR was represented as an interval-observed endpoint because the biological onset of response was not continuously observed but was detected at laboratory monitoring visits. For each patient, intervals before the first DMR-positive visit were coded as at risk without event, the first interval ending in DMR was coded as the event interval, and subsequent intervals were excluded from first-event modeling. Patients without observed DMR by last follow-up were censored.

### Time Scale

For clinical description, follow-up time and visit gaps are reported in months. For the Bayesian model, all time variables were converted to years before fitting. Thus, \(u_{ij}\) denotes years from imatinib initiation at visit \(j\) for patient \(i\), and interval start, interval end, and interval gap were also analyzed in years. Model coefficients involving time therefore correspond to functions of \(\log(1+\mathrm{years})\), not to linear per-month effects.

### Joint Longitudinal--Interval Model

Let \(y_{ij}\) denote the observed log-MRD value and \(u_{ij}\) denote years from imatinib initiation. The latent longitudinal trajectory was modeled as
\[
\mu_{ij} =
\beta_0 + \beta_1 \log(1+u_{ij})
+ \beta_2\{\log(1+u_{ij})\}^2
+ \beta_3 I(\mathrm{bone\ marrow}_{ij})
+ b_{0i} + b_{1i}\log(1+u_{ij}).
\]
For non-floor observations, \(y_{ij}\sim N(\mu_{ij},\sigma_y^2)\). Observations reported at the assay floor \(c=-5.0\) were treated as left-censored, contributing
\[
P(y_{ij}\le c)=\Phi\{(c-\mu_{ij})/\sigma_y\}
\]
to the likelihood. This censoring formulation was used because floor-coded deep-response values indicate that the true molecular burden is at or below the reporting floor rather than exactly equal to the floor.

For the DMR process, let \(L_{ik}\), \(R_{ik}\), and \(\Delta_{ik}=R_{ik}-L_{ik}\) be the start, end, and length in years of at-risk interval \(k\) for patient \(i\), and let \(u_{ik}^{mid}=(L_{ik}+R_{ik})/2\). The interval hazard-rate component was specified as
\[
h_{ik} =
\exp\{\gamma_0 + \gamma_1\log(1+u_{ik}^{mid})
+\gamma_2\log(1+\Delta_{ik})
+\alpha\mu_i(u_{ik}^{mid})\}.
\]
The probability that first DMR occurred during the interval was
\[
P(d_{ik}=1\mid \mathrm{at\ risk}) =
1-\exp(-h_{ik}\Delta_{ik}).
\]
The association parameter \(\alpha\) links the latent molecular trajectory to interval DMR probability. Because lower log-MRD indicates deeper molecular response, a negative \(\alpha\) indicates that lower latent molecular burden is associated with a higher interval probability of DMR.

### Priors, Computation, and Model Checking

The priors were weakly informative and chosen to regularize the nonlinear longitudinal trajectory, interval event model, and patient-specific effects. The current fitted model used \(N(-2.5,2^2)\) for \(\beta_0\), \(N(-1,1^2)\) for the linear log-time coefficient, \(N(0,0.5^2)\) for the quadratic log-time coefficient, \(N(0,1^2)\) for the bone marrow sample-source coefficient, and an exponential(1) prior for \(\sigma_y\). The interval model used \(N(-2,2^2)\) for \(\gamma_0\), \(N(0,1^2)\) for time and gap coefficients, and \(N(-0.5,0.75^2)\) for the latent MRD association parameter. Random-effect standard deviations had exponential(1) priors, standardized random effects had standard normal priors, and the random-effect correlation matrix used an LKJ Cholesky prior with shape 2.

The model was fitted using four Markov chains. Convergence and computation were assessed using divergent transitions, maximum treedepth, approximate E-BFMI, R-hat, effective sample size, and inspection of weakly estimated parameters. Model adequacy was assessed using posterior predictive checks for non-floor longitudinal log-MRD values, posterior floor probabilities for assay-floor observations, and calibration summaries comparing estimated DMR probabilities with observed DMR at the interval and patient levels.

## E. Revised Manuscript-Ready Results Text

The final processed analysis dataset contained 87 patients, 495 longitudinal molecular observations, and 275 at-risk intervals for first observed DMR. Endpoint definitions were audited by recomputing DMR and CMR directly from the final log-MRD values. DMR was defined as log-MRD \(\leq -4.5\), and CMR was defined as log-MRD \(\leq -5.0\). The stored endpoint indicators matched the recomputed values exactly. In the retained analysis data, DMR and CMR counts were identical because no observations fell in the DMR-only range \((-5.0,-4.5]\). All 239 DMR-positive longitudinal observations were also CMR-positive floor-coded observations at log-MRD = -5.0. Thus, the identical DMR and CMR summaries reflect the retained data distribution and assay-floor reporting pattern rather than a mismatch in endpoint coding. One excluded raw record had IS = 0.002 with missing ABL, sample type, laboratory date, and New LOG-MRD; adjudication of this record would be required before any attempt to recover it for analysis.

The fitted Bayesian joint model used time in years, after conversion from the processed month-scale variables. Descriptive summaries are reported in months for clinical readability. The Bayesian sampler produced 4000 post-warmup draws from four chains. There were no divergent transitions, the maximum treedepth was 7, the minimum approximate E-BFMI was 0.500, and the maximum R-hat was 1.019. The main fixed-effect and association parameters showed acceptable convergence. The random-effect correlation was weakly estimated and was not interpreted clinically.

The longitudinal submodel showed a strong nonlinear decline in latent log-MRD over follow-up. The posterior mean for the log-time coefficient was -3.539, with 95% posterior interval -4.497 to -2.613. The quadratic log-time coefficient was positive on average, with posterior mean 0.499 and 95% posterior interval -0.005 to 1.012, supporting nonlinear response dynamics. The bone marrow sample-source coefficient was positive on average but uncertain, with posterior mean 0.613 and 95% posterior interval -0.826 to 2.065; it was therefore treated as a measurement-source adjustment rather than as a confirmed clinical effect.

The interval DMR component showed a strong association between latent molecular burden and first observed DMR. The association parameter was negative, with posterior mean -1.256 and 95% posterior interval -1.851 to -0.809. Under the fitted parameterization, lower latent log-MRD was associated with higher interval probability of DMR. Because the model uses a nonlinear time transformation and an interval hazard formulation, this association should be interpreted as a model-based monitoring association rather than as a validated clinical decision threshold.

Posterior predictive checks supported cautious interpretation of the longitudinal component. Among non-floor observations, posterior predictive coverage was 0.961 for the 90% predictive interval and 0.988 for the 95% predictive interval. The observed floor rate was 0.483, while the posterior mean floor probability was 0.414. Calibration was descriptive because it was assessed in the development cohort. At the interval level, the observed event rate was 0.247 and the mean predicted probability was 0.247, with Brier score 0.082. At the patient level, the observed DMR rate was 0.782 and the mean predicted cumulative DMR probability was 0.581, with Brier score 0.124. This patient-level underprediction indicates that model-estimated probabilities should not be used as treatment-decision thresholds without recalibration and external validation.

## F. Short Limitations Paragraph

This study should be interpreted as development and internal checking of a monitoring-oriented model rather than validation of a treatment-decision rule. The cohort was modest, with 87 patients, incomplete core covariates, and no external validation cohort. In the retained analysis data, DMR and CMR indicators were identical because all observed DMR-positive measurements were floor-coded at log-MRD = -5.0; this limits the ability to distinguish DMR from deeper response in the current dataset. One excluded raw record with low IS but missing essential modeling fields requires adjudication before final analysis. The current model includes patient-specific random intercepts and slopes, but the random-effect correlation was weakly estimated and should not be interpreted clinically. Patient-level calibration showed underprediction of cumulative DMR probability, so model-based probabilities require recalibration and external validation before clinical deployment.

## Repair Priority

1. Add endpoint audit language to the manuscript immediately.
2. Standardize Methods notation to years for the fitted model and months for descriptive reporting.
3. Add the priors table and diagnostics table.
4. Add a data-adjudication note for the excluded IS = 0.002 raw record.
5. Decide whether to re-fit with independent random intercept and slope. If yes, regenerate all posterior estimates, PPCs, calibration summaries, and figures.
6. Keep all clinical claims framed as monitoring support and model development, not as treatment-free remission decision support.
