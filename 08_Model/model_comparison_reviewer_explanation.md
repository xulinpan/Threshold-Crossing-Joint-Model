# Reviewer-Facing Model Comparison and Validation Explanation

## Bottom Line

The primary Bayesian joint longitudinal--interval model remains justified as the main analysis because it is the only fitted model that simultaneously handles serial latent MRD, interval-observed DMR onset, irregular visit gaps, patient-level heterogeneity, and left-censored assay-floor observations. Its interval Brier score was 0.082 and its patient-level Brier score was 0.116.

## What Each Comparator Tests

- Kaplan--Meier: useful descriptive benchmark only; it treats first observed DMR visit time as the event time and therefore does not respect interval-observed onset.
- Interval timing + visit gap: tests whether interval-aware modeling alone is sufficient without longitudinal MRD.
- Landmark MRD: tests whether a simpler clinically familiar landmark approach captures enough information; evaluated only when prior landmark MRD was actually available and the patient remained at risk.
- Joint exact-floor model: tests whether the assay-floor censoring likelihood matters by fitting the same joint structure while treating floor observations as exact -5.0 values.
- Primary joint model: tests the full proposed monitoring framework with latent MRD, random patient effects, interval timing, visit gap, and left-censored floor observations.

## Internal Validation

The refittable simpler models used patient-cluster bootstrap optimism correction with 200 bootstrap resamples. For the two Stan joint models, the script reports patient-cluster bootstrap stability of the fixed posterior predictions; it does not claim full optimism correction because that would require repeated Stan refitting.

## Recommendation

Use the left-censored Bayesian joint longitudinal--interval model as the primary model. Present Kaplan--Meier as descriptive, the timing-gap interval model and landmark MRD model as simpler clinical benchmarks, and the exact-floor joint model as a sensitivity analysis showing the effect of assay-floor handling. Avoid claiming that the model is a validated treatment-decision rule; describe it as a monitoring-oriented development model requiring external validation.

## Files Generated

- `model_comparison_performance.csv`
- `table_07_model_comparison_validation.tex`
- `figure_10_model_comparison_interval_calibration.pdf/.png`
- `figure_11_model_comparison_patient_calibration.pdf/.png`
- `model_comparison_internal_validation.tex`
- `model_comparison_reviewer_explanation.md`
