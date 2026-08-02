# Clinical-Statistical Interpretation of the Joint Model Results

Prepared: 2026-07-09

## Executive Conclusion

The model results support a meaningful clinical application contribution: the study converts serial molecular monitoring data into a dynamic estimate of molecular response and links that latent trajectory to the probability of deep molecular response (DMR). The strongest evidence is not that the model is already a clinical decision rule, but that it provides a statistically coherent framework for individualized monitoring in real-world CML care.

The clinically defensible conclusion is:

> The joint model shows a marked decline in latent MRD over follow-up and a strong association between lower latent MRD burden and higher interval probability of DMR. This supports the clinical value of using the full longitudinal MRD history, rather than isolated measurements or exact-time approximations, to characterize evolving molecular response.

## Clinical Estimand

The clinically relevant question is not simply whether a baseline marker predicts DMR. The better question is:

> Given a patient's serial MRD pattern, irregular monitoring schedule, and current latent molecular burden, what is the evolving probability that the patient reaches DMR during follow-up?

This is why the joint longitudinal-interval model is clinically useful. It targets the monitoring process that clinicians actually observe:

- repeated BCR-ABL1 or MRD measurements;
- irregular visit timing;
- assay-floor values at deep response levels;
- DMR onset known between visits rather than at an exact date;
- patient-to-patient heterogeneity in baseline burden and response speed.

## Data Structure Supporting the Model

The real-data application includes:

| Quantity | Value |
|---|---:|
| Patients | 87 |
| Longitudinal MRD observations | 495 |
| At-risk DMR intervals | 275 |
| DMR event patients | 68 |
| Censored patients | 19 |
| Median follow-up | 39.0 months |
| Median visits per patient | 5 |
| Median visit gap | 6.0 months |
| Visit gaps over 6 months | 50.3% |
| Visit gaps over 12 months | 20.0% |
| Floor observations at log-MRD <= -5 | 239 |
| Complete core covariates | 71.3% |

This structure justifies the model choice. Nearly half of the longitudinal observations are at or below the assay floor, and DMR onset is naturally interval-observed. A standard exact-time survival model would ignore both of these clinically important measurement features.

## Methodological Interpretation

The longitudinal submodel estimates a latent log-MRD trajectory:

- fixed nonlinear time trend;
- bone marrow versus other sample source effect;
- patient-specific random intercept and slope;
- left-censoring for assay-floor log-MRD observations.

The DMR submodel estimates interval DMR probability:

- interval midpoint time;
- visit gap length;
- the current latent molecular burden;
- patient-specific latent trajectory carried from the longitudinal submodel.

The key statistical parameter is `alpha_mrd`, the association between latent log-MRD and DMR interval hazard. Since lower log-MRD means deeper molecular response, a negative `alpha_mrd` means that lower molecular burden is associated with higher DMR probability.

## Computation and Diagnostic Quality

The current Stan fit is adequate for substantive interpretation:

- 4 chains;
- 1000 warmup and 1000 sampling iterations per chain;
- 4000 post-warmup draws;
- 0 divergent transitions;
- maximum treedepth 7;
- minimum approximate E-BFMI 0.500;
- maximum R-hat 1.019.

The main fixed-effect and event association parameters have acceptable R-hat and effective sample sizes. The only notable weakness is the random-effect correlation estimate:

| Parameter | Mean | 95% posterior interval | R-hat | Bulk ESS |
|---|---:|---:|---:|---:|
| `Omega_b[1,2]` | 0.026 | -0.417 to 0.621 | 1.019 | 191 |

This should be treated as a weakly identified nuisance parameter, not a clinical finding.

## Main Posterior Results

The table below reports posterior means and 95% posterior intervals from the saved chain CSVs.

| Parameter | Mean | 95% posterior interval | Clinical interpretation |
|---|---:|---:|---|
| `beta_time` | -3.539 | -4.497 to -2.613 | Strong decline in latent log-MRD over follow-up. |
| `beta_time2` | 0.499 | -0.005 to 1.012 | Curvature is positive on average; 95% interval nearly excludes zero. |
| `beta_bm` | 0.613 | -0.826 to 2.065 | Bone marrow measurements are higher on average, but uncertainty is substantial. |
| `sigma_y` | 1.815 | 1.633 to 2.013 | Residual measurement variability remains meaningful. |
| `gamma_time` | -0.215 | -1.062 to 0.601 | Calendar-time effect in the interval model is uncertain after accounting for latent MRD. |
| `gamma_gap` | -1.837 | -2.890 to -0.790 | Visit-gap effect is negative in the interval hazard component, but should be interpreted cautiously because interval length also enters the cumulative hazard. |
| `alpha_mrd` | -1.256 | -1.851 to -0.809 | Lower latent MRD is strongly associated with higher interval probability of DMR. |
| `tau_b[1]` | 0.964 | 0.332 to 1.521 | Meaningful heterogeneity in patient-specific baseline MRD burden. |
| `tau_b[2]` | 2.212 | 1.550 to 2.963 | Meaningful heterogeneity in patient-specific response slope. |
| `Omega_b[1,2]` | 0.026 | -0.417 to 0.621 | Random-effect correlation is weakly estimated and should not be emphasized. |

## Longitudinal MRD Finding

The longitudinal submodel indicates a clear decline in latent log-MRD over follow-up. The time coefficient is strongly negative (`beta_time = -3.539`, 95% posterior interval -4.497 to -2.613), while the quadratic log-time term is positive on average (`beta_time2 = 0.499`, 95% posterior interval -0.005 to 1.012). Clinically, this suggests a rapid early decline followed by a changing rate of response over time rather than a simple linear trajectory.

The 95% interval for `beta_time2` nearly touches zero, so the safest claim is that the dominant feature is strong decline over time, with evidence of nonlinear curvature. If using the 5th to 95th posterior interval from the Stan summary, the curvature is positive (0.084 to 0.917).

The sample-source effect is uncertain. The posterior mean for bone marrow sampling is positive (`beta_bm = 0.613`), but the 95% interval crosses zero (-0.826 to 2.065). This should be described as adjustment for measurement source rather than a confirmed clinical effect.

## DMR Association Finding

The strongest clinical result is the association parameter:

| Quantity | Estimate |
|---|---:|
| `alpha_mrd` mean | -1.256 |
| 95% posterior interval | -1.851 to -0.809 |
| Multiplicative effect for 1-unit higher log-MRD | 0.285 |
| 95% interval for 1-unit higher log-MRD effect | 0.157 to 0.445 |
| Multiplicative effect for 1-unit lower log-MRD | 3.51 |
| 95% interval for 1-unit lower log-MRD effect | 2.25 to 6.36 |

This means that, under the model parameterization, a one-unit lower latent log-MRD value is associated with an approximately 3.5-fold higher interval DMR hazard component. This is clinically coherent: deeper molecular response is associated with higher DMR probability.

This result is the central bridge between the statistical method and clinical application. It shows that the latent biomarker trajectory is not merely descriptive; it is prognostically linked to clinically meaningful response.

## Predicted Interval DMR Probability

The generated interval probabilities show substantial clinical heterogeneity:

| Summary of interval-level posterior mean DMR probabilities | Value |
|---|---:|
| Number of intervals | 275 |
| Mean | 0.247 |
| Median | 0.136 |
| Interquartile range | 0.039 to 0.383 |
| Minimum | near 0 |
| Maximum | 0.964 |
| Intervals with mean probability < 0.01 | 51 |
| Intervals with mean probability > 0.80 | 20 |

This is important clinically because it suggests the model can distinguish low-probability and high-probability monitoring intervals. That supports future use as a monitoring-support or risk-stratification framework. However, calibration and external validation are still needed before using these probabilities for treatment decisions.

## Patient Heterogeneity

The random-effect standard deviations are nontrivial:

- `tau_b[1] = 0.964`, indicating patient-level heterogeneity in baseline latent log-MRD;
- `tau_b[2] = 2.212`, indicating patient-level heterogeneity in response trajectory slope.

Clinically, this supports individualized modeling. Patients do not share one common response curve. Some have faster molecular decline and others have slower or more variable decline. This is exactly the setting where dynamic longitudinal modeling adds value over a single baseline predictor.

The random-effect correlation is close to zero and weakly estimated, so it should not be interpreted as evidence that higher baseline burden is or is not linked to faster decline.

## Benchmark Model Context

The benchmark models are consistent with the joint-model story.

The interval-censored Weibull benchmark model using baseline log-MRD showed that higher baseline log-MRD was associated with longer time to DMR:

- molecular-only model: baseline log-MRD coefficient 0.247, p = 0.025;
- complete-case clinical model: baseline log-MRD coefficient 0.181, p = 0.067;
- complete-case clinical model used 62 patients.

This supports the clinical direction of the joint model: greater molecular burden is associated with delayed or less likely DMR.

The longitudinal mixed model also supports substantial time-varying MRD change:

- 495 observations;
- 87 patients;
- random intercept and random slope;
- nonlinear time terms with strong evidence for MRD decline;
- sample source effect positive but not conventionally significant.

The benchmark models are not the primary contribution. They serve as quality-control and triangulation evidence showing that the data support the broader joint-model conclusions.

## Clinical Application Value

The model contributes to clinical application in five ways.

1. It uses the full MRD history.

   The model avoids reducing a patient's course to baseline MRD or a single landmark value. It uses repeated molecular measurements to estimate a latent disease-control trajectory.

2. It respects real-world monitoring.

   Real CML follow-up is irregular. DMR is detected between clinic visits, not at a precisely observed biological onset time. The interval model directly represents that observation process.

3. It handles assay-floor measurements.

   Deep responses often produce floor-limited values. Treating these as exact measurements can distort the longitudinal trajectory. The left-censoring component is clinically and statistically appropriate.

4. It supports individualized monitoring intensity.

   The model can separate intervals with very low predicted DMR probability from intervals with high predicted DMR probability. This is relevant for designing follow-up intensity, patient counseling, and future decision-support tools.

5. It provides a foundation for dynamic prediction.

   The current model can be developed into a tool that updates DMR probability as new MRD measurements arrive. That is a clear clinical application pathway.

## Claim Boundaries

Strong claims supported by the current results:

- Serial MRD trajectories contain clinically meaningful information about DMR.
- A joint longitudinal-interval model is appropriate for irregular CML monitoring data.
- Lower latent MRD is strongly associated with higher interval DMR probability.
- Patient-specific heterogeneity supports individualized molecular-response modeling.
- The model is a clinically oriented statistical framework for monitoring and risk stratification.

Claims to avoid:

- The model is ready to guide treatment changes.
- The model is validated for TKI discontinuation decisions.
- The estimated random-effect correlation has clinical meaning.
- The model proves causal effects of monitoring frequency or visit gaps.
- The single-cohort result is externally generalizable without validation.

## Expert Clinical-Statistical Interpretation

From a clinical-statistical perspective, the value of this work is the alignment between the estimand, the data-generating process, and the clinical monitoring problem. The model does not force the data into an exact-event survival framework. Instead, it recognizes that DMR is observed through scheduled molecular monitoring and that the underlying molecular state is only partially observed through noisy and floor-limited measurements.

The dominant clinical signal is the latent MRD trajectory. Patients with lower estimated molecular burden have substantially higher interval probability of DMR. This finding is statistically stable, clinically coherent, and consistent with established CML response biology. The model therefore provides a stronger clinical interpretation than a baseline-only analysis: it shows how evolving molecular response, not just initial status, relates to deep response.

The most important limitation is validation. The current analysis demonstrates clinical plausibility and internal coherence, but it does not yet establish transportability to other cohorts or readiness for prospective clinical decision-making. The appropriate current framing is a clinically motivated joint modeling framework with potential for dynamic monitoring support.

## Manuscript-Ready Results Paragraph

Among 87 patients, the analysis included 495 longitudinal MRD observations and 275 at-risk DMR intervals, with 68 patients reaching DMR and 19 censored. The cohort showed substantial real-world monitoring complexity, including a median follow-up of 39.0 months, median visit gap of 6.0 months, 50.3% of visit gaps longer than 6 months, and 239 observations at or below the log-MRD assay floor. These features supported use of a joint longitudinal-interval model with assay-floor censoring.

The Bayesian joint model showed adequate computation for substantive interpretation, with 4000 post-warmup draws, no divergent transitions, maximum treedepth of 7, and minimum approximate E-BFMI of 0.500. The longitudinal submodel showed a strong decline in latent log-MRD over follow-up (`beta_time = -3.539`, 95% posterior interval -4.497 to -2.613), with evidence of nonlinear curvature (`beta_time2 = 0.499`, 95% posterior interval -0.005 to 1.012). The sample-source effect was positive on average but uncertain (`beta_bm = 0.613`, 95% posterior interval -0.826 to 2.065).

The interval DMR submodel showed a strong association between latent molecular burden and DMR probability. The association parameter was negative (`alpha_mrd = -1.256`, 95% posterior interval -1.851 to -0.809), indicating that lower latent log-MRD was associated with higher interval probability of DMR. A one-unit lower latent log-MRD value corresponded to an approximately 3.5-fold higher interval DMR hazard component. Interval-level posterior mean DMR probabilities varied widely, with median 0.136 and interquartile range 0.039 to 0.383, supporting clinically meaningful heterogeneity in dynamic response probability.

## Manuscript-Ready Discussion Paragraph

These findings support the clinical value of a joint modeling approach for CML molecular monitoring. The model integrates serial MRD measurements, irregular visit timing, assay-floor behavior, and interval-observed DMR onset into a single framework. This structure is better aligned with real-world clinical follow-up than analyses based only on baseline MRD, landmark response, or exact-time approximations of DMR onset. The strong association between lower latent MRD and higher DMR probability supports the biological coherence of the model and suggests potential use for individualized monitoring and future dynamic prediction. However, the current results should be interpreted as evidence for a clinically oriented modeling framework rather than as a validated treatment-decision tool. External validation, calibration assessment, and prospective evaluation are needed before using model-based probabilities to guide therapeutic decisions.

