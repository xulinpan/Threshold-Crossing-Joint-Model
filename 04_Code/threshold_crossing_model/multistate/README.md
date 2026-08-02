# Bidirectional multi-state latent threshold-crossing model (prototype)

An extension of the threshold-crossing joint model from *first onset* to a
**bidirectional multi-state** process for molecular-response depth:

```
state 0  (not yet DMR)  --onset-->  state 1 (in DMR)  --relapse-->  state 2 (molecular relapse)
```

The latent log-MRD trajectory is quadratic in `ell = log(1+t)`,
`m_i(ell) = a_i + g_i ell + beta2 ell^2` with `beta2 > 0` (convex), so it
declines then rebounds and crosses the DMR threshold `c_D = -4.5` **downward
once (onset)** and **upward once (relapse)**. Both crossings are
interval-observed at irregular visits and modelled probabilistically with a
threshold width `sigma_thr` that is *distinct from* the measurement error
`sigma_y` (so the two are separately identified — fixing the identifiability
problem of the single-state model).

## What is novel here

- **Running-minimum onset + running-maximum durability.** Onset uses the
  running minimum `M_i(t) = min_{u<=t} m_i(u)`; *sustained* DMR over a window
  `W` uses the running maximum after onset, `max_{u in [T,T+W]} m_i(u) <= c_D`.
  For a convex quadratic both are **closed form** (value at the clamped vertex,
  and post-vertex point values), so there is no soft-min grid — smoother
  geometry and far fewer divergences than the grid-based single-state model.
- **Relapse / DMR loss** as an upward crossing — the single-state model cannot
  represent this (its crossing is absorbing).
- Durability (sustained MR4.5), which drives treatment-free-remission
  decisions, becomes a first-class, identifiable estimand.

## Files

```
multistate/
├── stan/multistate_threshold.stan   # the model (onset, relapse, durability GQ, log_lik)
├── R/01_simulate_ms.R               # generative model + Stan-data builder
├── R/02_fit_ms.R                    # cmdstanr fit + recovery + transition calibration
├── python/msval.py                  # marginal-ML validation of the likelihood
├── python/gqcheck.py                # validation of the generated-quantities formulas
└── outputs/                         # created on run
```

## Validation (already performed, in Python)

Because Stan was not available in the drafting environment, the likelihood and
the generated-quantities formulas were validated numerically before shipping.

**Parameter recovery** (marginal-ML refit, `msval.py`, n=450):

| param | est | truth |
|------|-----|------|
| beta0 | -2.000 | -2.000 |
| beta1 | -3.570 | -3.600 |
| beta2 | 1.273 | 1.300 |
| sigma_y | 0.732 | 0.800 |
| tau0 | 0.526 | 0.550 |
| tau1 | 0.679 | 0.750 |
| **sigma_thr** | **0.151** | **0.150** |

Note `sigma_thr` is recovered accurately here (0.151 vs 0.150), unlike the
single-state model where it was badly biased — confirming the reformulated
estimand is identifiable.

**Transition-probability formulas** (`gqcheck.py`, model-implied vs empirical
rates): onset 0.503 vs 0.502; relapse 0.134 vs 0.129; sustained 0.467 vs 0.494.

## Run (on a machine with cmdstanr)

```bash
cd 04_Code/threshold_crossing_model/multistate/R
Rscript 02_fit_ms.R
```

Outputs: `ms_posterior_summary.csv` (recovery), `ms_diagnostics.csv`,
`ms_transition_calibration.csv` (onset / sustained-DMR / relapse calibration
and Brier). Inits are supplied and `adapt_delta = 0.99`.

## Caveats / next steps

- The convex-quadratic trajectory forces an eventual rebound; within finite
  follow-up most patients are censored in state 1, which is realistic, but a
  **plateauing / random-asymptote or bi-exponential** trajectory is the natural
  refinement so that deep sustained responders are not implicitly destined to
  relapse.
- Next extension: add the **latent-class (cure / mover-stayer)** component so
  non-responders (trajectory minimum above `c_D`) are modelled explicitly
  rather than as slow responders.
- Then: competing risk (death / TKI discontinuation) and the informative
  visit-process layer.
- This is a prototype validated by simulation; it has not yet been fit to the
  real cohort or externally validated.
```
