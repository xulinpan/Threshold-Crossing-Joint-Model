# A Bayesian Threshold-Crossing Joint Model for Irregular, Assay-Floor-Censored Molecular Monitoring in Chronic Myeloid Leukemia

Reproducibility materials for the manuscript of the same name.

The model treats deep molecular response (DMR) in chronic myeloid leukemia as a
**threshold crossing of a latent log-MRD trajectory** rather than as a separately
coded event. Because DMR is *defined* as log-MRD ≤ −4.5, predicting a coded DMR
event from the same biomarker is circular; reading the endpoint as a
first-hitting time of the trajectory removes that circularity. The longitudinal
sub-model treats assay-floor values as **left-censored**, and the endpoint is
**interval-observed** because DMR is ascertained only at monitoring visits.

A conventional interval-censored hazard (with association parameter α) and a
multi-state extension (onset, durability, relapse) are options of the *same*
Stan program, sharing one trajectory and one set of random effects.

---

## ⚠️ Patient data are **not** included

The clinical cohort cannot be shared. It comprises identifiable patient health
information (the source records contain patient names) and is subject to
institutional and national privacy regulation.

**Deliberately excluded from this repository:**

| Excluded | Reason |
|---|---|
| `03_Data/Raw/*` (`glw.csv`, `glw_data.xlsx`, `PH+染色体-Table 1.csv`) | contain a patient-name column |
| `03_Data/Processed/private/patient_key.csv` | maps coded IDs → patient names |
| `real_longitudinal_analysis.csv`, `real_interval_survival_analysis.csv`, `real_patient_level_analysis.csv`, `real_patient_level_complete_covariates.csv` | de-identified but still the real cohort |
| `stan_data_real_*.rds` | serialised real-cohort model input |
| patient-level result files (per-record predictions, calibration detail) | derived from real patients |

**Included instead:** a fully **synthetic** dataset of identical structure
(`data/simulated_*.csv`, patient IDs `S0001…`), generated from the fitted
model's posterior predictive distribution, plus the generating parameters
(`data/simulation_true_parameters.csv`). The entire pipeline runs end-to-end on
the synthetic data and reproduces the qualitative features of the analysis.

Aggregate results (posterior summaries, simulation operating characteristics,
calibration bins, model comparison) are included in `results/` because they
contain no patient-level rows.

---

## Contents

```
R/                          analysis pipeline (00–17)
  00–03  inspect, clean, build Stan data
  04–07  prepare joint-model data, benchmarks, fit, sensitivity/PPC
  08–12  reporting, model comparison, time-function, visit process, dynamic prediction
  13     c_F = −4.5 censoring-point sensitivity
  14     serial-correlation (OU) residual model
  15     multi-state HMC credible intervals
  16     unified single-program model, fitted as a ladder (M1→M2→M3→M3s)
  17     LOO diagnostics (r_eff, Pareto-k)
stan/                       Stan models
  glw_joint_interval_dmr*.stan          interval-hazard joint models
  glw_joint_interval_dmr_serial.stan    + OU serial residual
  glw_unified_joint.stan                unified program (flag-switched likelihood)
  glw_latent_threshold_crossing_skeleton.stan
threshold_crossing_model/   standalone threshold-crossing study (R, Stan, Python)
data/                       synthetic data, data dictionary, builders
results/                    aggregate fitted results and generated LaTeX tables
figures/                    manuscript figures (vector PDF)
```

## Requirements

- **R** ≥ 4.2 with `cmdstanr`, `posterior`, `loo`, `dplyr`, `tidyr`, `ggplot2`, `ggpubr`, `survival`, `survminer`
- **CmdStan** ≥ 2.32 (a working C++ toolchain is required)
- **Python** ≥ 3.9 with `numpy`, `pandas` (for the synthetic-data generator only)

Compiled Stan binaries (`*.exe`) are intentionally **not** distributed; models
compile on first use.

## Reproducing the analysis (synthetic data)

```r
setwd("<repo root>")
source("R/03_build_stan_data.R")     # build Stan inputs from data/simulated_*.csv
source("R/06_fit_stan_joint_model_independent.R")
source("R/16_fit_unified_model.R")   # unified ladder: M1 → M2 → M3 → M3s
source("R/17_loo_diagnostics.R")     # LOO with r_eff and Pareto-k audit
```

Scripts resolve paths relative to the project root. Sampling the full ladder
takes on the order of hours on four cores.

**Note.** Because the shipped data are synthetic, numerical output will *not*
equal the manuscript values; the manuscript's fitted quantities are provided as
aggregate summaries in `results/`.

## Model ladder (`R/16_fit_unified_model.R`)

| Fit | Link | Spline | Transitions |
|---|---|---|---|
| M1 | interval hazard (α) | – | onset |
| M2 | threshold crossing (σ_thr) | – | onset |
| M3 | threshold crossing | – | onset + durability + relapse |
| M3s | threshold crossing | ✓ | onset + durability + relapse |

M1 reproduces the conventional joint model, which is what establishes that the
specifications differ in what they assume rather than in how they were fitted.

## Citation

See `CITATION.cff`. Please cite both the manuscript and the archived software
release (Zenodo DOI).

## License

Code is released under the MIT License (`LICENSE`). Documentation, figures and
aggregate result tables are released under CC-BY-4.0.
