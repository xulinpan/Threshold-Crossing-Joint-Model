# Dynamic DMR Prediction and Recalibration

Generated: 2026-07-13 22:10:16 EDT

## Scope

This analysis uses posterior draws from the renewed Bayesian joint longitudinal--interval model to generate landmark-level probabilities of documented DMR by 6, 12, and 24 months after each landmark. Predictions are evaluated among patients who were DMR-free at the landmark, had at least one MRD measurement at or before the landmark, and had observable outcome status at the prediction horizon.

Important limitation: these are development-cohort dynamic predictions using the already fitted posterior. The patient-specific random effects in the saved posterior were estimated using the full longitudinal dataset, so the numeric results should be interpreted as an apparent/prototype analysis. A strict prospective dynamic prediction analysis requires landmark-specific posterior updating using only MRD history available at time s, or leave-one-patient-out refitting.

## Landmark and Horizon Definitions

- Landmarks: 6, 12, 18, 24 months.
- Horizons: 6, 12, 24 months.
- Eligible records: DMR-free at landmark, at least one prior MRD measurement, and documented DMR by horizon or follow-up beyond the horizon.

## Pooled Validation Summary

- 6-month horizon: n=108, events=39, original mean predicted=0.506 vs observed=0.361, original Brier=0.179; recalibrated mean predicted=0.361, recalibrated Brier=0.157, optimism-corrected recalibrated Brier=0.167.
- 12-month horizon: n=99, events=57, original mean predicted=0.626 vs observed=0.576, original Brier=0.162; recalibrated mean predicted=0.576, recalibrated Brier=0.161, optimism-corrected recalibrated Brier=0.172.
- 24-month horizon: n=93, events=77, original mean predicted=0.721 vs observed=0.828, original Brier=0.081; recalibrated mean predicted=0.828, recalibrated Brier=0.055, optimism-corrected recalibrated Brier=0.063.

## Output Files

- `dynamic_prediction_landmark_predictions.csv`
- `dynamic_prediction_recalibration_coefficients.csv`
- `dynamic_prediction_validation.csv`
- `dynamic_prediction_bootstrap_optimism.csv`
- `dynamic_prediction_calibration_bins.csv`
- `dynamic_prediction_decision_curve.csv`
- `table_09_dynamic_prediction_validation.tex`
- `figure_13_dynamic_prediction_calibration.pdf/.png`
- `figure_14_dynamic_prediction_decision_curve.pdf/.png`

