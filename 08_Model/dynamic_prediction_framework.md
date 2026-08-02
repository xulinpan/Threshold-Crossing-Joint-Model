# Dynamic DMR Prediction and Recalibration Framework

## A. Dynamic Prediction Algorithm

The prediction target is:
\[
\Pr(T_i^{\mathrm{DMR}}\le s+H \mid T_i^{\mathrm{DMR}}>s,\mathcal H_i(s)),
\]
where \(s\) is a landmark time, \(H\) is the prediction horizon, and \(\mathcal H_i(s)\) denotes all MRD measurements available up to \(s\).

### Landmark Times

Recommended landmark times:

- Primary: 6 and 12 months.
- Exploratory: 18 and 24 months.

The later landmarks are clinically meaningful but sparse in this cohort after requiring that patients remain DMR-free and have observable horizon status.

### Prediction Horizons

Recommended horizons:

- 6 months;
- 12 months;
- 24 months.

All three horizons are supported when pooling across eligible landmarks. Landmark-specific estimates after 18 and 24 months should be interpreted cautiously because some cells contain fewer than 20 patients.

### Eligibility at Landmark

A patient-landmark-horizon record is included only if:

1. the patient had not yet documented DMR by landmark \(s\);
2. the patient had at least one MRD measurement at or before \(s\);
3. the binary outcome was observable, meaning either DMR was documented by \(s+H\) or follow-up extended beyond \(s+H\).

The observed outcome is documented DMR by \(s+H\), not true biological DMR onset, because true onset remains interval-observed.

### Prediction Steps

1. Select landmark \(s\) and horizon \(H\).
2. Identify eligible patients.
3. For each posterior draw \(m\), compute the patient-specific latent MRD trajectory using the joint model.
4. Convert the predicted latent trajectory into horizon-specific DMR probability.
5. Average over posterior draws to obtain the individualized original dynamic prediction:
   \[
   \hat p_i(s,H)=M^{-1}\sum_{m=1}^{M}p_i^{(m)}(s,H).
   \]
6. Fit horizon-specific logistic recalibration:
   \[
   \mathrm{logit}\{p_{i,\mathrm{cal}}(s,H)\}
   =
   a_H+b_H\mathrm{logit}\{\hat p_i(s,H)\}.
   \]
7. Evaluate original and recalibrated predictions using Brier score, calibration intercept, calibration slope, grouped calibration plots, and exploratory decision-curve analysis.
8. Use clustered bootstrap resampling by patient to estimate optimism in the recalibration layer.

## B. Mathematical Prediction Formula

Let \(x=t/12\) denote treatment time in years when time is recorded in months. The latent biological trajectory excluding sample-source adjustment is:
\[
\eta_i(t)=
\beta_0+\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2
+b_{0i}+b_{1i}\log(1+t).
\]

For a continuous-time prediction target, the ideal dynamic probability is:
\[
p_i(s,H\mid \theta,b_i)
=
1-\exp\left\{-\int_s^{s+H}h_i(u\mid\theta,b_i)\,du\right\}.
\]

Using the current interval-hazard parameterization, the implementation approximates the horizon as one prediction interval:
\[
t^{mid}=s+\frac{H}{2},
\]
\[
h_i(s,H\mid\theta,b_i)
=
\exp\{\gamma_0+\gamma_1\log(1+t^{mid})
+\gamma_2\log(1+H)+\alpha\eta_i(t^{mid})\},
\]
\[
p_i(s,H\mid\theta,b_i)
=
1-\exp\{-Hh_i(s,H\mid\theta,b_i)\}.
\]

For planned clinical monitoring schedules, the same formula can be generalized to multiple future intervals:
\[
p_i(s,H)=
1-\prod_{k: (L_k,R_k]\subset(s,s+H]}
\{1-p_{ik}\},
\]
where each \(p_{ik}\) uses the interval midpoint and planned visit gap.

## C. R Implementation Plan

The implementation is in:

`04_Code/R/12_dynamic_prediction_recalibration.R`

The script:

1. reads the processed patient-level and longitudinal datasets;
2. reads only the needed posterior columns from the renewed CmdStan draw CSVs;
3. constructs eligible landmark-horizon records for \(s=6,12,18,24\) months and \(H=6,12,24\) months;
4. computes posterior mean individualized DMR probabilities;
5. fits horizon-specific recalibration models;
6. computes dynamic Brier score, calibration intercept, and calibration slope;
7. generates grouped calibration plots;
8. computes an exploratory decision-curve analysis;
9. estimates optimism for the recalibration layer using clustered bootstrap resampling by patient.

Important technical note: the numeric analysis uses the existing fitted posterior. The patient-specific random effects in these draws were estimated from the full longitudinal record, not refitted using only MRD history available at landmark \(s\). Therefore, the generated metrics are development/prototype metrics. A strict prospective dynamic prediction analysis requires landmark-specific posterior updating of \(b_i\) based only on \(\mathcal H_i(s)\), or leave-one-patient-out refitting.

## D. Validation Table

The generated validation table is:

`08_Model/table_09_dynamic_prediction_validation.tex`

The corresponding CSV is:

`08_Model/dynamic_prediction_validation.csv`

Pooled horizon results from the development analysis were:

| Horizon | Model | n | Events | Observed | Mean predicted | Brier | Optimism-corrected Brier | Calibration intercept | Calibration slope |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 6 months | Original posterior | 108 | 39 | 0.361 | 0.506 | 0.179 | 0.179 | -1.00 | 0.85 |
| 6 months | Recalibrated | 108 | 39 | 0.361 | 0.361 | 0.157 | 0.167 | 0.00 | 1.00 |
| 12 months | Original posterior | 99 | 57 | 0.576 | 0.626 | 0.162 | 0.162 | -0.36 | 0.81 |
| 12 months | Recalibrated | 99 | 57 | 0.576 | 0.576 | 0.161 | 0.172 | 0.00 | 1.00 |
| 24 months | Original posterior | 93 | 77 | 0.828 | 0.721 | 0.081 | 0.081 | 1.02 | 1.81 |
| 24 months | Recalibrated | 93 | 77 | 0.828 | 0.828 | 0.055 | 0.063 | 0.00 | 1.00 |

Optimism correction applies only to the recalibration layer, because the full Bayesian joint model was not refitted in each bootstrap sample.

## E. Calibration Plot Design

The calibration plot is:

`08_Model/figure_13_dynamic_prediction_calibration.pdf`

and the 600 dpi PNG is:

`08_Model/figure_13_dynamic_prediction_calibration.png`

Design:

- Panel 1: original posterior dynamic predictions.
- Panel 2: recalibrated dynamic predictions.
- Separate curves for 6-, 12-, and 24-month horizons.
- Points are quantile-based prediction bins.
- Point size reflects bin size.
- The 45-degree line indicates perfect calibration.

The calibration-bin data are in:

`08_Model/dynamic_prediction_calibration_bins.csv`

Decision-curve outputs are exploratory:

- `08_Model/dynamic_prediction_decision_curve.csv`
- `08_Model/figure_14_dynamic_prediction_decision_curve.pdf`
- `08_Model/figure_14_dynamic_prediction_decision_curve.png`

Decision-curve analysis should not be emphasized in the manuscript unless a clinically defensible intervention threshold is specified.

## F. Manuscript-Ready Results Section

Use the companion LaTeX insert:

`08_Model/dynamic_prediction_methods_results.tex`

## G. Limitation Statement

The dynamic prediction analysis remains an internal development analysis. First, the saved posterior draws use patient-specific random effects estimated from the full longitudinal record, so the numeric analysis may be optimistic relative to a truly prospective landmark prediction setting. Second, recalibration was performed in the same cohort and only bootstrap optimism correction for the recalibration layer was available; the full Bayesian model was not refitted in bootstrap or leave-one-patient-out samples. Third, DMR onset is documented only at monitoring visits, so the horizon outcome is documented DMR by \(s+H\), not exact biological DMR onset. Fourth, the cohort is small, especially at later landmarks, and no external validation cohort is available. Therefore, recalibrated probabilities should be interpreted as development-stage monitoring-support estimates and should not be used as clinical decision thresholds without external validation and prospective calibration.
