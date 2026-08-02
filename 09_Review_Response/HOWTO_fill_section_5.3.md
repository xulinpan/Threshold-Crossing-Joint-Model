# How to run the misspecification scenarios and fill §5.3

The only genuinely incomplete section of the manuscript. Script:
`04_Code/threshold_crossing_model/R/12_misspecification_study.R`

## The design (why it answers the question)

Data are generated under a **departure** from the working model; the **working
model is then fitted unchanged** (piecewise-quadratic trajectory, Gaussian
errors, single assay floor, schedule-conditional). Everything reported is
therefore the cost of fitting the working model when it is wrong — which is
exactly what "robustness to misspecification" means.

| Scenario | What is broken |
|---|---|
| `correct` | nothing — reference row |
| `traj_monotone` | quadratic term removed; monotone decline, no upturn |
| `traj_exp` | exponential approach to an asymptote (sharper early decline) |
| `threshold_m5` | DMR threshold at −5.0 instead of −4.5 |
| `hetero_floor` | patient-specific assay floors; lab still *reports* −5.0 |
| `heavy_t3` | *t*₃ measurement errors **and** *t*₃ random effects |
| `informative` | response-adaptive visit intensity (visits depend on latent MRD) |

## Step 1 — smoke test (~10 minutes)

Edit the header of the script:

```r
N_REP    <- 5        # tiny
USE_MCMC <- FALSE    # fast MAP path, no coverage
```

then

```r
setwd("<repo>/04_Code/threshold_crossing_model/R")
source("12_misspecification_study.R")
```

You are checking that all seven scenarios generate and fit, and that
`outputs/misspec_summary.csv` appears with sensible bias columns. Coverage will
be meaningless at this stage — ignore it.

**If `fit_one()` is not found**, `source("04_simulation_study.R")` first (that is
where it is defined) or move the definition into `02_fit_models.R`.

## Step 2 — production run

```r
N_REP    <- 200      # matches the 200 replicates used elsewhere
USE_MCMC <- TRUE     # required for honest credible-interval coverage
ADAPT    <- 0.99
```

**Cost.** 7 scenarios × 200 replicates × 4 chains ≈ 5,600 HMC fits. On 4 cores
this is realistically **overnight to a couple of days**. Two ways to cut it:

- **Shard it.** Your repo already has `run_shards.ps1`; split the replicate range
  across processes and `rbind` the `misspec_raw.csv` files afterwards.
- **Reduce to 100 replicates.** Monte Carlo SE on a coverage near 0.95 is then
  ≈ 0.022 instead of 0.015 — still adequate, and the script reports the MCse so
  the precision is explicit.

Do **not** drop below ~100: with 30–50 replicates the coverage MCse (~0.03–0.04)
is too wide to distinguish 0.95 from 0.88, which is the comparison that matters.

## Step 3 — insert the table

The script writes `outputs/table_misspec.tex`. Copy it to `08_Model/` (where the
other generated tables live), then in §5.3 replace the "runs are ongoing"
paragraph with:

```latex
Table~\ref{tab:sim-misspec} reports the results. [One or two sentences of
interpretation — see the template below.]

\input{table_misspec.tex}
```

## Step 4 — restore the stronger claims

Once the table exists, the hedges I inserted can be reverted **if the results
support them**. Currently three places say the association is stable "across the
censoring point and across prior families". If α holds up under all seven
scenarios, restore:

> stable across the censoring point, prior families and every misspecification examined

in the **Abstract** (~line 83), **Introduction** (~line 215), and **Conclusions**
(~line 1366). Also delete the sentence in §5.3 beginning *"we therefore make no
claim in this paper about robustness to structural misspecification"*.

**Only do this if the numbers actually support it.** If α coverage drops under
some departure, say so — a robustness section that reports a failure is more
credible than one that reports uniform success.

## Interpretation template

Fill from the table rather than from expectation:

> The biomarker–event association $\alpha_{\mathrm{MRD}}$ was [robust / degraded]
> across the departures studied, with $95\%$ credible-interval coverage between
> [x] and [y]. The longitudinal variance components were, as expected for a
> Gaussian piecewise-quadratic working model, sensitive to the departures that
> directly violate their assumptions: under [scenario] the random-slope SD
> $\tau_{b1}$ was under-covered ([z]), and under heavy-tailed $t_3$ errors the
> residual SD $\sigma_y$ [result]. In short, [the association is well calibrated
> under all departures studied / the following departures compromise it: …],
> whereas the variance components should be interpreted cautiously when the
> trajectory shape or the error distribution is grossly wrong.

## What to watch for

1. **`heavy_t3` will likely hurt σ_y most.** That is the policy-relevant result,
   given the manuscript already reports σ_y ≈ 1.81 as implausibly large for pure
   assay error. If coverage collapses here, it belongs in the Abstract.
2. **`informative` is the one that could undermine a headline claim.** §6.4
   asserts that conditioning on the observed visit schedule is reasonable. If α
   coverage degrades under informative monitoring, §6.4 must be rewritten.
3. **Report the convergence rate honestly.** An earlier draft of this table
   carried a mislabelled failure-rate column that implied 100% failure. The
   script now emits `conv_rate` = proportion meeting all diagnostics; do not
   invert it.
4. **Do not fabricate.** If the runs do not finish before submission, leave the
   current honest hedging in place — it is defensible as written.
