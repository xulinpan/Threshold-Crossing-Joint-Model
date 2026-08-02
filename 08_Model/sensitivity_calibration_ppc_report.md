# Sensitivity, Calibration, and Posterior Predictive Checks

Generated: 2026-07-09

## Bottom Line

The added checks support the manuscript's main conclusion: the fitted joint model is clinically coherent and adequate for cautious applied interpretation. The posterior predictive checks show reasonable longitudinal behavior, and calibration summaries show useful separation of low- and high-risk DMR intervals/patients. These results strengthen the case for submission to BMC Medical Research Methodology, provided the claims remain focused on monitoring support rather than validated treatment decisions.

## Posterior Predictive Check for Longitudinal MRD

- Observations: 495
- Floor observations: 239
- Observed floor rate: 0.483
- Posterior mean floor probability: 0.414
- Non-floor RMSE to posterior mean trajectory: 1.492
- Non-floor 90% predictive coverage: 0.961
- Non-floor 95% predictive coverage: 0.988

Interpretation: the model captures the dominant longitudinal pattern and explicitly recognizes that many deep-response observations are floor-limited. The floor probability check is especially important because treating floor values as exact values would understate measurement uncertainty at deep response.

## DMR Probability Calibration

- Interval-level observed event rate: 0.247
- Interval-level mean predicted probability: 0.247
- Interval-level Brier score: 0.082
- Patient-level observed DMR rate: 0.782
- Patient-level mean predicted DMR probability: 0.581
- Patient-level Brier score: 0.124

Interpretation: calibration is descriptive because the same data were used for fitting and checking. Interval-level aggregate calibration is close: expected and observed interval events are nearly identical. Patient-level cumulative DMR probability is lower than the observed patient-level DMR rate, which should be reported as a calibration limitation and a reason not to use these probabilities as clinical decision thresholds without recalibration or external validation. The grouped calibration tables are still useful for reviewers because they show whether predicted probabilities separate lower- and higher-risk monitoring intervals and patients.

## Assay-Floor Sensitivity

Three longitudinal benchmark variants were fitted: floor values treated as exact -5, floor values shifted to -5.5, and non-floor observations only. The goal is not to replace the Bayesian censored model, but to show whether the qualitative longitudinal signal is robust to simple floor-handling choices.

- Floor sensitivity outputs: `sensitivity_assay_floor_summary.csv`, `sensitivity_assay_floor_longitudinal_coefficients.csv`, and `sensitivity_assay_floor_fixed_predictions.csv`.

Interpretation: the two full-data floor variants support the same qualitative time-response pattern. The non-floor-only variant is a stress test rather than a direct replacement model: after removing floor observations, the richer mixed model is singular and the script falls back to a fixed-effect model, with sample source not estimable because of rank deficiency. This should be described transparently if included in the supplement.

## Visit-Gap and Interval-Model Sensitivity

Discrete-time logistic sensitivity models were fitted to the interval records with and without visit-gap terms and baseline/clinical adjustment. These are benchmark sensitivity models and should be presented as supporting analyses, not as replacements for the primary joint model.

- Best interval-logistic sensitivity model by AIC: `complete_case_clinical` (AIC 189.716).
- Interval sensitivity outputs: `sensitivity_interval_logistic_summary.csv` and `sensitivity_interval_logistic_coefficients.csv`.

## Generated Figures

- `figure_08_posterior_predictive_longitudinal.pdf` and `.png`
- `figure_09_dmr_calibration.pdf` and `.png`

## Manuscript-Ready Text

Posterior predictive checks supported the adequacy of the longitudinal component of the joint model. Among non-floor observations, predictive coverage was acceptable, and the model explicitly represented assay-floor behavior through posterior floor probabilities rather than treating deep-response values as exact. Calibration summaries were then used to compare model-estimated DMR probabilities with observed DMR frequencies at both interval and patient levels. These calibration analyses are descriptive because they are based on the development cohort. Interval-level aggregate calibration was close, whereas patient-level cumulative DMR probability was lower than the observed DMR rate, reinforcing that the current model should be used for monitoring-oriented interpretation rather than uncalibrated decision thresholds.

Sensitivity analyses evaluated whether the main interpretation depended on simple modeling choices. Longitudinal benchmark models were refitted under alternative assay-floor treatments, and interval-level logistic benchmark models were fitted with and without visit-gap terms and baseline clinical adjustment. These analyses supported the robustness of the central conclusion that serial MRD trajectories contain clinically meaningful information about DMR, while also reinforcing that the fitted model should be interpreted as a monitoring-support framework requiring external validation before clinical decision use.
