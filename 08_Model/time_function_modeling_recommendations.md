# Time-Function Comparison and Recommendation

Generated: 2026-07-09 16:03:52 EDT

## Scope

This file compares candidate longitudinal time functions for the CML molecular-monitoring model. The empirical comparison is a screening benchmark using Gaussian mixed models; it is not a replacement for the Bayesian joint model, because assay-floor observations are treated as observed at -5.0 in this screening step. Any final manuscript claim about posterior predictive performance or DMR calibration requires refitting the selected Bayesian joint model with the left-censoring likelihood.

## Empirical Screening Results

- Best AIC in this screening benchmark: Low-rank spline proxy (4 df) (AIC 1723.4).
- Lowest non-floor RMSE: Log-quadratic in log(1+t) (RMSE 1.095).
- Lowest interval-level Brier score in the two-stage calibration proxy: Piecewise linear knots 3, 6, 12, 24, 60 months (Brier 0.123).
- Lowest patient-level Brier score in the two-stage calibration proxy: Log-quadratic in log(1+t) (Brier 0.179).

## Recommendation

Recommended fitted primary time function for the current manuscript version: retain the current log-quadratic function unless the Bayesian joint model is fully refitted with the spline time trend. It is parsimonious, interpretable, has already passed HMC diagnostics in the renewed Stan fit, and performed best for non-floor RMSE and patient-level Brier score in this screening benchmark.

Recommended renewal candidate if the model is refitted: a low-rank penalized spline for the population time trend, with independent patient-level random intercept and random log-time slope. The spline is the best flexibility candidate because it had the lowest AIC in the screening benchmark and can represent rapid early decline and later plateau without forcing one global quadratic curvature. It should replace the current primary model only if the full Bayesian refit preserves convergence, posterior predictive performance, and interval- and patient-level DMR calibration.

The clinically knotted piecewise-linear model should be reported as an interpretability sensitivity analysis. A monotone-decline-constrained spline should be framed as a biological sensitivity analysis, preferably with a soft rather than hard population-level constraint so that genuine late increases in MRD are not masked.

## Output Files

- `time_function_comparison_summary.csv`
- `time_function_interval_predictions.csv`
- `time_function_population_trends.csv`
- `table_08_time_function_comparison.tex`
- `figure_12_time_function_population_trends.pdf`
- `figure_12_time_function_population_trends.png`

## Required Reanalysis Before Submission

1. Refit the full Bayesian joint model using the selected penalized spline time function.
2. Keep the assay-floor likelihood as left-censored at log-MRD <= -5.0.
3. Compare posterior predictive checks against the current log-quadratic model, especially early follow-up, floor probability, and late follow-up behavior.
4. Recompute interval-level and patient-level DMR calibration using posterior predicted event probabilities.
5. Prefer the spline as primary only if convergence diagnostics are acceptable and the calibration/PPC checks are at least as good as the current primary model.

