# Stan Result Check

Checked: 2026-07-09

## Bottom Line

The full Stan model has been sampled and the current outputs are usable for interpretation, with one caution.

The main fixed-effect and event-model parameters have good convergence diagnostics. The only material warning is the random-effect correlation estimate, which has R-hat 1.019 and bulk ESS 191. This is not a blocker for the main substantive conclusions, but it should be reported as a weakly estimated nuisance parameter or improved with a longer run if the correlation itself matters.

## Current Artifacts

- Stan fit RDS: `08_Model/stan_joint_interval_dmr_fit.rds`
- Stan summary CSV: `08_Model/stan_joint_interval_dmr_summary.csv`
- CmdStan chain CSVs: `08_Model/stan_joint_interval_dmr_draws/`
- Stan input data: `03_Data/Processed/stan_data_real_joint_interval_dmr.rds`
- Model code: `04_Code/Stan/glw_joint_interval_dmr.stan`

The chain CSVs were recovered from the live temporary paths referenced by the fit object and copied into the durable project folder because the expected draws directory was missing.

## Sampling Setup

- Chains: 4
- Warmup per chain: 1000
- Sampling draws per chain: 1000
- Total post-warmup draws: 4000
- Save warmup: false
- HMC metric: diagonal Euclidean

## Sampler Diagnostics

- Divergences: 0
- Maximum treedepth reached: 7
- Mean acceptance by chain: 0.915, 0.944, 0.930, 0.933
- Minimum approximate E-BFMI across chains: 0.500
- Max R-hat in summary: 1.019
- Parameters with R-hat > 1.01: 4 rows
- Parameters with R-hat > 1.05: 0 rows
- Parameters with bulk ESS < 400: 3 rows
- Parameters with tail ESS < 400: 0 rows

Rows triggering convergence or ESS warnings:

| Variable | R-hat | Bulk ESS | Tail ESS |
|---|---:|---:|---:|
| `L_b[2,1]` | 1.019 | 191 | 858 |
| `Omega_b[1,2]` | 1.019 | 191 | 858 |
| `Omega_b[2,1]` | 1.019 | 191 | 858 |
| `lp__` | 1.012 | 430 | 1235 |

## Core Posterior Results

Intervals below are posterior 5th and 95th percentiles from the current summary CSV.

| Parameter | Mean | 5th pct | 95th pct | R-hat | Bulk ESS |
|---|---:|---:|---:|---:|---:|
| `beta0` | -2.085 | -3.366 | -0.853 | 1.000 | 2137 |
| `beta_time` | -3.539 | -4.330 | -2.768 | 1.002 | 2559 |
| `beta_time2` | 0.499 | 0.084 | 0.917 | 1.001 | 2583 |
| `beta_bm` | 0.613 | -0.591 | 1.858 | 1.000 | 2281 |
| `sigma_y` | 1.815 | 1.659 | 1.982 | 1.001 | 3345 |
| `gamma0` | -4.036 | -6.082 | -2.254 | 1.001 | 1778 |
| `gamma_time` | -0.215 | -0.934 | 0.467 | 1.000 | 4982 |
| `gamma_gap` | -1.837 | -2.694 | -0.979 | 1.002 | 4359 |
| `alpha_mrd` | -1.256 | -1.733 | -0.864 | 1.002 | 1068 |
| `tau_b[1]` | 0.964 | 0.472 | 1.416 | 1.006 | 587 |
| `tau_b[2]` | 2.212 | 1.646 | 2.844 | 1.008 | 637 |
| `Omega_b[1,2]` | 0.026 | -0.368 | 0.523 | 1.019 | 191 |

## Interpretation Check

- The longitudinal trajectory terms are strongly supported: `beta_time` is negative and `beta_time2` is positive, suggesting a nonlinear decline pattern in log MRD over time.
- The bone-marrow sample effect `beta_bm` is positive on average but its 90 percent posterior interval crosses zero.
- The interval event model supports a negative association between latent MRD level and DMR interval hazard through `alpha_mrd`, with the current parameterization: lower predicted MRD corresponds to higher DMR probability.
- Longer visit gaps have a negative interval coefficient in this fitted parameterization (`gamma_gap`), so interpretation should be tied carefully to how `event_interval` and interval risk are encoded.
- The random-effect correlation is close to zero but weakly estimated; do not emphasize it as a substantive finding.

## Benchmark Models

Benchmark outputs also exist and are internally consistent:

- Interval model patients: 87
- Complete-case adjusted interval model patients: 62
- Longitudinal observations: 495
- Longitudinal patients: 87
- Molecular interval Weibull AIC: 290.85
- Complete-case clinical interval Weibull AIC: 192.80
- Longitudinal mixed model AIC: 1741.06

## Data Verification

The existing verification report passes all listed checks:

- Longitudinal rows: 495
- Patients: 87
- Interval rows: 275
- Patient-level rows: 87
- DMR event patients: 68
- Required fields, time ordering, coded IDs, privacy scan, figures, tables, PDF report, and Stan RDS checks: PASS

## Recommendation

Use the current Stan results for the manuscript's main model interpretation, but keep the claims focused on the well-estimated fixed effects and event association. If final submission depends on the random-effect correlation, rerun with more sampling iterations or a stronger correlation prior and regenerate the summary.
