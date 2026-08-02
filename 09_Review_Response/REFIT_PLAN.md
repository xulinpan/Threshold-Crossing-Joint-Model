# Refit plan and status — working in `paper01_glw` (real source)

Date: 26 July 2026. Target manuscript: `02_LaTeX/glw_full_manuscript_smmr_rev7b.tex`.

## Scopes 1 & 2 — already satisfied in the real source (confirmed)

- **Real vector figures**: rev7b uses `\graphicspath{{../05_Figures/}{../08_Model/}}` and
  `\includegraphics{figure_01…17.pdf}` — the authoritative vector figures, not raster crops.
- **Exact numbers**: rev7b wires tables via `\input{table_0x…}`, resolved from `06_Tables/`
  and `08_Model/` (the generated `.tex` tables built from the fitted results).
- **Build verified**: compiled cleanly to **26 pages, 0 undefined references, no missing
  figures** with `TEXINPUTS="./:../06_Tables//:../08_Model//:"`.

No action needed here — the raster/approximate versions were only ever in the `Paper_02`
copy. Ethics/TRIPOD live in the submission package (`submission_SMMR/01_title_page.tex`,
`TRIPOD_checklist.md`), which is why the anonymized main text has no inline ethics
paragraph — that is correct for double-blind.

## Scope 3 — the three refits (scripts written; must be run in your R/Stan)

R/Stan cannot execute in this environment, so these are ready-to-run scripts wired to the
real data and your Stan conventions. Your repo already has cmdstan configured (compiled
`.exe` models are present).

### 1. Censoring point c_F = −4.5 — `04_Code/R/13_sensitivity_censoring_point.R`
- Loads `stan_data_real_joint_interval_dmr.rds`, sets `floor_value = -4.5`, refits the
  primary independent model. `is_floor` is unchanged (no observations in (−5.0, −4.5]), so
  only the left-censored term moves — the change is exactly one data value.
- Output: `08_Model/sensitivity_cF_minus4p5_summary.csv` (+ diagnostics, fit rds).
- Wire into: the Prior-sensitivity / censoring discussion in §4; report β_time, τ_b1,
  σ_y, α_MRD, and the floor-probability PPC at −4.5 vs −5.0.

### 2. Serial-correlation residual — `04_Code/Stan/glw_joint_interval_dmr_serial.stan` + `04_Code/R/14_fit_serial_correlation.R`
- Adds a within-patient continuous-time AR(1)/Ornstein–Uhlenbeck residual u(t) to the
  **observation** mean only; the biological trajectory driving the interval hazard stays
  smooth. This lets σ_y become pure measurement error.
- The R script reorders obs by (patient, time) and builds `first_in_patient`, `prev_obs`,
  `dt_prev`, then fits at `adapt_delta = 0.995`.
- Output: `08_Model/serial_correlation_summary.csv` (params include σ_y, σ_u, ℓ, σ_total).
- **Success check**: σ_y should drop toward the assay CV (~0.2–0.3 on log10) while τ_b1 and
  α_MRD stay stable. If so, this is the fix for the audit's "σ_y is not measurement error".
- **STATUS: new model — validate on simulated data before trusting real-cohort results.**

### 3. Multi-state HMC intervals — `04_Code/R/15_multistate_hmc_intervals.R`
- A full HMC threshold fit already exists (`04_Code/threshold_crossing_model/outputs/
  fit_threshold.rds`). This script extracts proper **95%** credible intervals to replace the
  asymptotic intervals in the multi-state table.
- Output: `08_Model/multistate_hmc_intervals_95.csv`.
- **Two caveats flagged by the script**:
  1. `hmc_diagnostics.csv` shows **219 divergent transitions** on the threshold fit —
     re-run at `adapt_delta ≥ 0.995` (and reparameterise) before quoting the intervals.
  2. `posterior_summary_threshold.csv` carries a `truth` column, i.e. it is a
     **simulation-recovery** run. Confirm `fit_threshold.rds` is the real-cohort fit
     (`primary_cohort.rds`) before wiring into the manuscript; refit on the real cohort if not.

## Suggested run order

```r
setwd("<path to>/paper01_glw")
source("04_Code/R/13_sensitivity_censoring_point.R")   # c_F = -4.5
source("04_Code/R/14_fit_serial_correlation.R")        # serial residual (validate first)
source("04_Code/R/15_multistate_hmc_intervals.R")      # multi-state 95% CrIs + audit
```

Then update the corresponding `06_Tables/`/`08_Model/` `.tex` tables (or the inline
numbers) and recompile rev7b with the TEXINPUTS above.

## Remaining audit items not addressed by these scripts

- Colour-vision safety of Figures (re-export `figure_*.pdf` from `04_Code/R/02…05` with a
  CVD-safe palette + redundant line/shape).
- Broken KM/interval-only comparators (negative scaled Brier) — re-examine
  `04_Code/R/05_fit_benchmark_models.R`.
- LOO-ELPD for quadratic-vs-spline and correlated-vs-independent (add to `09_model_comparison_validation.R`).
