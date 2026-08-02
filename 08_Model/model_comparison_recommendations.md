# Recommended Primary and Sensitivity Models

## Recommended Primary Model

Use the renewed Bayesian joint longitudinal--interval model with independent patient-specific random intercept and slope terms and left-censored assay-floor likelihood as the primary manuscript model.

Rationale:

- It had the lowest interval-level Brier score: 0.082.
- It had the lowest patient-level Brier score: 0.116.
- It is the only model that simultaneously represents latent longitudinal MRD, irregular at-risk intervals, visit gap, patient-level heterogeneity, and assay-floor censoring.
- It directly matches the clinical data-generating process: DMR is observed at visits, not continuously, and floor-level MRD values are not exact continuous measurements.

## Recommended Sensitivity Models

1. Joint exact-floor model: include as the main assay-floor sensitivity analysis. This model keeps the same joint interval structure but treats floor observations as exact -5.0 values.

2. Interval timing + visit-gap model: include as the parsimonious interval-aware benchmark without longitudinal MRD. This isolates the gain from adding serial molecular burden.

3. Landmark MRD model: include as a clinically familiar benchmark, but clearly state that it was evaluated only where a prior 6-, 12-, or 18-month landmark was available and the patient remained at risk.

4. Kaplan--Meier curve: retain as descriptive only. It uses first observed DMR visit time and should not be presented as a valid exact-onset survival model.

## Recommended Manuscript Claim

The model-comparison results support the complex joint model as the primary development model for monitoring-oriented inference. They do not establish a validated treatment-decision rule. External or temporal validation is still required before clinical deployment.
