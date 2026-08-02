# A Bayesian Threshold-Crossing Joint Model for Irregular, Assay-Floor-Censored Molecular Monitoring in Chronic Myeloid Leukemia

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21753041.svg)](https://doi.org/10.5281/zenodo.21753041)

Code, models and reproducible workflow for the manuscript of the same name.

Deep molecular response (DMR) in CML is *defined* as the BCR::ABL1 transcript
crossing −4.5 on the log₁₀ scale. Modelling it as a separate, conditionally
independent survival outcome is therefore circular: the same measurements enter
the likelihood twice, once as biomarker and once as event. This repository
implements a joint model that reads DMR instead as a **first-hitting time of the
latent trajectory**, together with the simulation studies that establish its
operating characteristics.

---

## ⚠ Data availability

**The clinical data are not in this repository and cannot be shared.** The
cohort is identifiable patient health information from a single centre and is
subject to institutional and national privacy regulation. `.gitignore` excludes
the raw extracts, the analysed patient records, and — critically — the
re-identification keys.

To make the pipeline runnable, `04_Code/R/19_make_synthetic_cohort.R` generates a
**fully synthetic cohort with the identical column schema**, from the posterior
means of the fitted model. No patient data is read or used. Every script accepts
`prefix = "synthetic"` wherever it accepts `prefix = "real"`.

The synthetic data verify that the code runs and produces output of the right
shape. **They will not reproduce the reported estimates.** The reported numbers
are recorded in `06_Tables/` and `10_Final/reproducibility_package/results/`.

---

## Quick start

```bash
# 1. dependencies (R >= 4.1, CmdStan >= 2.30)
Rscript 04_Code/R/00_install_dependencies.R      # if present; else install cmdstanr, posterior, loo, dplyr

# 2. generate the synthetic cohort
Rscript 04_Code/R/19_make_synthetic_cohort.R

# 3. fit the primary joint model
Rscript 04_Code/R/06_fit_stan_joint_model_independent.R
```

Run everything from the **project root**, not from `04_Code/R/` — several scripts
resolve paths relative to the root.

---

## What is where

| Path | Contents |
|---|---|
| `04_Code/Stan/` | Stan models for the application |
| `04_Code/R/` | Data preparation, fitting, sensitivity and reporting scripts |
| `04_Code/threshold_crossing_model/` | Self-contained simulation study (see below) |
| `06_Tables/`, `05_Figures/` | Generated tables and figures |
| `02_LaTeX/` | Manuscript source |
| `08_Model/` | Fit objects and posterior draws (excluded from the archive; regenerable) |

### Stan models

| File | Role |
|---|---|
| `glw_joint_interval_dmr_independent.stan` | **Primary model.** Independent random effects, Gaussian errors. Produces Table 5. |
| `glw_joint_interval_dmr_independent_studentt.stan` | Student-*t* observation model, ν estimated. Table 16. |
| `glw_joint_interval_dmr.stan` | Correlated random effects (LKJ). Not the primary analysis. |
| `glw_joint_interval_dmr_exact_floor.stan` | Floor values as exact, for the like-for-like censoring contrast. |
| `glw_joint_interval_dmr_serial.stan` | Adds an OU/AR(1) residual. |
| `glw_unified_joint.stan` | Unified single-process formulation. |

### Simulation study

`04_Code/threshold_crossing_model/` is a **simulation-only sandbox**. It never
reads patient data; its "primary cohort" is synthetic by construction. All
data-generating code lives in one file, `R/dgm_common.R`, whose behaviour is set
by four environment variables (`DGM_VERSION`, `TRUTH_SOURCE`, `EVENT_MECHANISM`,
`PRIOR_SET`). The chosen combination is written into the checkpoint directory
name so that runs under different mechanisms can never be pooled.

| Script | Produces |
|---|---|
| `06_simulation_redesign.R` | Central + misspecification cells (Tables 11, 13) |
| `16_floor_vs_exact.R` | Paired assay-floor contrast (Table 10, Figure 15) |
| `17_calibration_study.R` | Patient-level calibration (Table 12, Figure 16) |
| `15_prior_sensitivity.R` | Prior-widening analysis (Table 14) |
| `10_sbc.R` | Simulation-based calibration (Table 15) |
| `13_dgm_sanity_check.R` | Compares the mechanism against cohort marginals |
| `07_aggregate_sim.R`, `14_make_tables.R` | Aggregation and LaTeX table generation |

Long runs are sharded and checkpointed per replicate; `./run_shards*.ps1`
launches them and resumes cleanly after interruption.

`04_simulation_study.R` and `12_misspecification_study.R` are retained but
**deprecated** — their headers explain why, and no result in the manuscript
comes from them.

---

## Reproducing the reported analyses

| Result | Command (from project root) |
|---|---|
| Primary fit, Table 5 | `Rscript 04_Code/R/06_fit_stan_joint_model_independent.R` |
| Student-*t* sensitivity, Table 16 | edit `stan_file` in the above to `..._studentt.stan`, re-run |
| Model comparison by LOO | `Rscript 04_Code/R/18_loo_gaussian_vs_studentt.R` |
| Simulation study | `cd 04_Code/threshold_crossing_model/R && ./run_shards.ps1` |
| Floor contrast | `./run_shards_floor.ps1` |
| Calibration study | `./run_shards_calib.ps1` |
| All simulation tables | `Rscript 14_make_tables.R` |

Expect roughly 10 hours for the main simulation study on 16 cores, 3–4 for the
floor contrast, 4 for the calibration study. The primary model fit takes minutes.

---

## Environment

Developed with R 4.5, CmdStan ≥ 2.30, `cmdstanr`, `posterior`, `loo`, `dplyr`.
Simulation scripts additionally use `plyr` (deprecated scripts only). PowerShell
launchers assume Windows; the R scripts are platform-independent.

Record your own session with `sessionInfo()` when reporting a discrepancy.

---

## Citation

See `CITATION.cff`. If you use the framework, please cite the manuscript; if you
use this code specifically, please also cite the archived release DOI.

## Licence

Code and Stan models: MIT (`LICENSE`). Manuscript text and figures: © the
authors, all rights reserved pending publication. No patient data is licensed
here because none is included.
