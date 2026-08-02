# Renewed Stan Result Check

Checked: 2026-07-09

## Bottom Line

The renewed primary model replaces the correlated random-intercept/random-slope structure with independent patient-specific random intercept and slope terms. This directly removes the weakly estimated correlation parameter from the original primary model while preserving the longitudinal floor-censoring likelihood and interval DMR event likelihood.

The final renewed fit used longer chains and is suitable for manuscript reporting: no divergent transitions were observed, no parameters exceeded R-hat 1.01, and no bulk or tail ESS values were below 400. The model should still be described as a monitoring-oriented development model, not a validated treatment-decision rule.

## Current Artifacts

- Stan fit RDS: `08_Model/stan_joint_interval_dmr_independent_renewed_fit.rds`
- Stan summary CSV: `08_Model/stan_joint_interval_dmr_independent_renewed_summary.csv`
- CmdStan chain CSVs: `08_Model/stan_joint_interval_dmr_independent_renewed_draws/`
- Renewed Stan code: `04_Code/Stan/glw_joint_interval_dmr_independent.stan`
- Previous correlated-model outputs archived in: `08_Model/archive_correlated_model_20260709_renewal/`

## Sampling Setup

- Chains: 4
- Warmup per chain: 2000
- Sampling draws per chain: 2000
- Total post-warmup draws: 8000
- Save warmup: false
- HMC metric: diagonal Euclidean

## Sampler Diagnostics

- Divergences: 0
- Maximum treedepth reached: 7
- Mean acceptance by chain: 0.939, 0.949, 0.946, 0.935
- Minimum approximate E-BFMI across chains: 0.610
- Max R-hat in summary: 1.003
- Minimum bulk ESS: 1359.258
- Minimum tail ESS: 1395.463
- Parameters with R-hat > 1.01: 0
- Parameters with bulk ESS < 400: 0
- Parameters with tail ESS < 400: 0

No parameters had R-hat > 1.01, bulk ESS < 400, or tail ESS < 400.

## Core Posterior Results

Intervals below are posterior 5th and 95th percentiles from the renewed summary CSV.

| Parameter | Mean | 5th pct | 95th pct | R-hat | Bulk ESS |
|---|---:|---:|---:|---:|---:|
| `beta_time` | -3.561 | -4.351 | -2.778 | 1.001 | 4623.043 |
| `beta_time2` | 0.502 | 0.078 | 0.914 | 1.000 | 5019.407 |
| `beta_bm` | 0.607 | -0.565 | 1.815 | 1.000 | 4293.185 |
| `sigma_y` | 1.813 | 1.660 | 1.976 | 1.000 | 6909.979 |
| `gamma_time` | -0.228 | -0.973 | 0.507 | 1.001 | 10428.085 |
| `gamma_gap` | -1.831 | -2.699 | -0.973 | 1.000 | 9709.047 |
| `alpha_mrd` | -1.275 | -1.724 | -0.899 | 1.001 | 3060.939 |
| `tau_b[1]` | 0.968 | 0.589 | 1.339 | 1.001 | 1359.258 |
| `tau_b[2]` | 2.203 | 1.657 | 2.815 | 1.002 | 1777.439 |

## Interpretation Check

- The longitudinal trajectory terms remain strongly supported: `beta_time` is negative and `beta_time2` is positive, indicating nonlinear decline in log-MRD over follow-up.
- The bone-marrow sample-source adjustment remains uncertain because its posterior interval crosses zero.
- The joint association remains clinically coherent: lower latent MRD is associated with higher interval probability of DMR through the negative `alpha_mrd` estimate.
- The visit-gap coefficient remains negative in this parameterization, so interpretation should be tied to the discrete interval hazard and cumulative interval length rather than read as a simple marginal effect.
- The random-effect correlation is no longer part of the primary model; this is the intended repair for the original overparameterized correlation component.

## Posterior Predictive and Calibration Checks

- Non-floor RMSE to posterior mean: 1.494
- Non-floor 90% predictive coverage: 0.957
- Non-floor 95% predictive coverage: 0.984
- Observed floor rate: 0.483
- Posterior mean floor probability: 0.414
- Interval observed event rate versus mean predicted probability: 0.247 versus 0.247
- Patient observed DMR rate versus mean predicted DMR probability: 0.782 versus 0.581

## Recommendation

Use the renewed independent-random-effects model as the primary manuscript model. Present the original correlated structure, if mentioned at all, as an audited predecessor that motivated simplification. Keep clinical claims focused on monitoring support, descriptive calibration, and methodological development under irregular molecular monitoring.
