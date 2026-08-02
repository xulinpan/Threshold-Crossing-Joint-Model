# Bayesian Hierarchical Threshold-Crossing Model — R/Stan pipeline

Reproducible R + Stan implementation of the model in
`bayesian_hierarchical_threshold_crossing_model.tex`: a latent patient-specific
log-MRD trajectory with assay-floor left-censoring, a probabilistic latent
threshold-crossing likelihood for interval-observed deep molecular response
(DMR), plus dynamic prediction and a simulation study.

This folder is **self-contained** and does not overwrite the existing
`04_Code/R` pipeline or the manuscript figures/tables.

## Which model is primary? (read this first)

The **primary model** is the **multi-state latent threshold-crossing joint
model with a linear-spline trajectory**:

- Stan: [`multistate/stan/multistate_spline.stan`](multistate/stan/multistate_spline.stan)
- Fit to the real CML cohort: [`multistate/R/03_fit_spline.R`](multistate/R/03_fit_spline.R)
- Trace/diagnostics: [`multistate/R/04_trace_plot.R`](multistate/R/04_trace_plot.R)

It defines DMR **onset** (first downward crossing), **durable response**
(running-maximum staying below threshold over a window) and **molecular relapse**
(subsequent upward crossing) through crossings of one latent trajectory, and — because
the trajectory is piecewise-quadratic — computes the running minimum and maximum
in **closed form** (no soft-min grid, divergence-free sampling).

Everything else is secondary:

| Role | Model / file |
|------|--------------|
| **Primary** | `multistate/stan/multistate_spline.stan` (multi-state, spline trajectory) |
| Base multi-state (quadratic) | `multistate/stan/multistate_threshold.stan` |
| **Comparator** | `stan/interval_hazard_joint.stan` (interval-hazard) |
| Development, **superseded** | `stan/threshold_crossing_joint.stan` (single-state, soft-min grid; ~219 divergences) |
| Exploratory, **not fitted** | `Stan/glw_joint_interval_dmr_visit_process_skeleton.stan` (informative visits) |

The single-state soft-min model (`stan/threshold_crossing_joint.stan`) and the
interval-hazard comparator are fit to a **simulated** cohort in
`R/02_fit_models.R`/`R/03_numeric_results.R` as a recovery study; they are
**not** the proposed method. Use the multi-state spline model for applied work.

## Layout

```
threshold_crossing_model/
├── stan/
│   ├── threshold_crossing_joint.stan   # primary model (soft-min running minimum)
│   └── interval_hazard_joint.stan      # documented-DMR comparator
├── R/
│   ├── 00_setup.R            # packages, palettes (RColorBrewer), paths, TRUTH
│   ├── 01_simulate_data.R    # generative model + Stan-data builders
│   ├── 02_fit_models.R       # (1) fit both Stan models with cmdstanr/rstan
│   ├── 03_numeric_results.R  # (2) posterior summaries, diagnostics, calibration, dynamic prediction
│   ├── 04_simulation_study.R # (3) repeated simulation: bias / Brier / calibration
│   ├── 05_figures.R          # (4) 600-dpi figures with RColorBrewer
│   └── run_all.R             # orchestrator
├── python/                   # earlier fast NumPy/SciPy simulation (reference)
└── outputs/                  # created on run: csv results + figures/ (600 dpi)
```

## Requirements

- R (>= 4.1)
- **cmdstanr** (preferred) with a working CmdStan (>= 2.30), or **rstan**
- `posterior`, `RColorBrewer` (and optionally `ragg` for PNG output)

Install CmdStan once:

```r
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
cmdstanr::install_cmdstan()
install.packages(c("posterior", "RColorBrewer", "ragg"))
```

## Run

```bash
cd 04_Code/threshold_crossing_model/R
Rscript run_all.R
```

Or step by step (`00_setup.R` is sourced automatically by each script):

```bash
Rscript 02_fit_models.R        # deliverable 1: Bayesian fit in Stan
Rscript 03_numeric_results.R   # deliverable 2: numeric results
Rscript 04_simulation_study.R  # deliverable 3: simulation study
Rscript 05_figures.R           # deliverable 4: high-quality figures
```

## Deliverables produced

1. **Stan Bayesian model** — `stan/threshold_crossing_joint.stan` (and the
   `interval_hazard_joint.stan` comparator), fitted by HMC in `02_fit_models.R`.
2. **Numeric results** (`outputs/`): `posterior_summary_threshold.csv`,
   `posterior_summary_interval.csv`, `hmc_diagnostics.csv`, `ppc_summary.csv`,
   `calibration_summary.csv`, `dynamic_prediction.csv`.
3. **Simulation study** (`outputs/`): `simulation_raw.csv`,
   `simulation_summary.csv` — bias, RMSE, Brier and calibration comparing
   assay-floor left-censoring vs exact-floor coding across n and monitoring.
4. **Figures** (`outputs/figures/`, 600 dpi, RColorBrewer):
   `figR_01_trajectories.png`, `figR_02_ppc_floor.png`,
   `figR_03_dmr_calibration.png`, `figR_04_dynamic_prediction.png`,
   `figR_05_sim_floor_bias.png`, `figR_06_sim_calibration.png`.

## Troubleshooting

**"Rejecting initial value: Log probability evaluates to log(0), i.e. negative
infinity."** Stan is starting from *random* parameter values at which the
threshold-crossing / interval-hazard likelihood underflows to zero (log = −∞),
so it cannot begin. This is fixed by supplying sensible initial values, which
`02_fit_models.R` and `04_simulation_study.R` now do via `init_threshold()` /
`init_interval()` (z0/z1 = 0, `sigma_thr` started at 0.4 to avoid Φ-difference
underflow). If you call the models yourself, always pass `init = ...`.

If rejections persist, in order of effectiveness: (i) confirm you passed the
`init` function; (ii) start `sigma_thr` a little larger (e.g. 0.5); (iii) lower
the soft-min sharpness `KAPPA` (60 → 30) in `00_setup.R`; (iv) as a last resort
add a small lower bound in the Stan file, `real<lower=0.02> sigma_thr;`.
Rejections that appear only briefly during warmup (not at initialization) are
harmless and can be ignored.

**"External command failed... No such file or directory ...Rtmp.../*.csv"**
(from `read_cmdstan_csv` in `03`/`05`). A cmdstanr fit object only stores
*pointers* to temporary CSV files, which are deleted between R sessions — so a
`saveRDS(fit)` cannot be reloaded later. `02_fit_models.R` now avoids this by
**materializing** the posterior into a `posterior::draws_df` and a diagnostics
data frame before saving (`res$draws`, `res$diag`); `03`/`05` read those
directly. If you save your own cmdstanr fits, use `fit$save_object(file=...)`
or `posterior::as_draws_df(fit$draws())` before `saveRDS()`.

## Notes

- Data are simulated from the fitted-model parameters (`TRUTH` in `00_setup.R`)
  so the whole pipeline is reproducible without patient data. To fit the real
  cohort, replace the `simulate_cohort()` call in `02_fit_models.R` with a
  loader that returns the same `list(long, pat, n)` structure.
- The simulation study uses `cmdstanr$optimize()` (posterior mode) per replicate
  for speed; set `USE_MCMC <- TRUE` in `04_simulation_study.R` for full sampling
  and interval coverage. Reduce `N_REP` for a quick check.
- Colour figures are designed to remain legible in grayscale (Dark2 palette +
  distinct point/line types), per SMMR artwork guidance.
```
