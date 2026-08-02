# Feasibility for Statistics in Medicine

Prepared: 2026-07-09

Question: does the GLW dataset support model innovation and submission to *Statistics in Medicine*?

## Short Answer

**Current state: not yet. Conditional future state: yes, if the paper becomes a true statistical-methods paper rather than an applied data-analysis paper.**

The GLW dataset is valuable and supports the motivation for a new or adapted model. It contains irregular longitudinal molecular monitoring, interval-observed first DMR, assay-floor log-MRD values, and incomplete covariates. These are exactly the kinds of complications that can motivate methodological work. However, *Statistics in Medicine* would likely expect more than a single-cohort application of a joint model. The manuscript would need a clear methodological innovation, simulation evidence, comparison with existing methods, and a fully fitted real-data example.

## Dataset Strength for a Methods Paper

The dataset supports a methodological paper because it has:

- 87 patients with model-ready longitudinal data.
- 495 longitudinal log-MRD observations.
- 275 at-risk survival intervals.
- 68 DMR event patients and 19 censored patients.
- Median follow-up of 39.0 months.
- Median visit gap of 6.0 months.
- 50.3% of visit gaps longer than 6 months.
- 20.0% of visit gaps longer than 12 months.
- 239 log-MRD observations at or below the reporting floor of -5.
- 62 patients with complete core covariates.

These features justify the model problem: **joint modeling of a longitudinal biomarker with interval-observed first response under irregular visit timing and assay-floor censoring.**

## Why the Current Dataset Alone Is Not Enough

The current dataset is a single real-world cohort. For *Statistics in Medicine*, this is usually not enough to claim a model innovation because:

1. The sample is modest: 87 patients.
2. Complete covariates are available for only 62 patients.
3. There is no external validation cohort.
4. The full Bayesian Stan model is specified but not yet fitted.
5. Benchmark models are fitted, but no full comparison against the proposed joint model is available yet.
6. There is no simulation study showing bias, coverage, calibration, or efficiency.
7. The proposed model is currently an adaptation/integration of known ideas: joint modeling, interval censoring, random effects, and left-censoring.

Therefore, the current paper is strong for **BMC Medical Research Methodology**, **BMC Bioinformatics**, or **PLOS ONE**, but it is still a stretch for *Statistics in Medicine*.

## What Would Make It Suitable for Statistics in Medicine

To make *Statistics in Medicine* realistic, the paper should be reframed as a statistical-methods contribution with the GLW cohort as the motivating and applied example.

### Required Upgrade 1: Define the Statistical Innovation

The innovation should be stated as something like:

> We propose an assay-floor-aware joint longitudinal-interval response model for irregularly monitored biomarker data, where the event time is interval-observed and the longitudinal biomarker contains lower-limit/floor observations.

The novelty should not be only "we applied a joint model to CML." It should be:

- a specific likelihood for interval-observed response events;
- left-censoring for assay-floor longitudinal measurements;
- explicit visit-gap dependence;
- shared random effects connecting latent biomarker trajectory to event probability;
- an implementation strategy for small, irregularly monitored clinical cohorts.

### Required Upgrade 2: Complete the Full Bayesian Model

The Stan model must be fitted and reported with:

- posterior means/medians;
- 95% credible intervals;
- convergence diagnostics, including R-hat and effective sample size;
- posterior predictive checks for log-MRD trajectories;
- posterior interval-DMR probability calibration;
- sensitivity analysis for priors and assay-floor handling.

### Required Upgrade 3: Add Simulation Studies

The manuscript needs simulation studies that mimic GLW-like data:

- small sample sizes, such as \(N=50, 100, 200\);
- irregular visit schedules;
- different censoring rates;
- different percentages of floor log-MRD observations;
- different strengths of association between latent log-MRD and DMR;
- missing covariates.

Key performance metrics:

- bias of association parameter \(\alpha\);
- coverage of credible intervals;
- root mean squared error;
- calibration of predicted DMR probability;
- sensitivity to interval width;
- comparison of exact-time mis-specified models vs interval-aware models.

### Required Upgrade 4: Compare Against Existing Approaches

At minimum, compare:

1. naive exact-time Cox/Kaplan-Meier using observed DMR visit date;
2. interval-censored Weibull model;
3. two-stage longitudinal plus survival model;
4. standard joint model ignoring assay-floor censoring;
5. proposed joint longitudinal-interval model with floor censoring.

The key result should show when and why the proposed model improves estimation or prediction.

### Required Upgrade 5: Clarify Generalizability

The method should be presented for a general class of problems:

> irregular biomarker monitoring with interval-observed threshold-crossing events and lower-limit measurement behavior.

CML should be the application, not the only reason the method exists.

## Possible Statistics in Medicine Manuscript Title

**A Joint Longitudinal-Interval Model for Irregular Biomarker Monitoring with Assay-Floor Measurements**

Subtitle/application:

**Application to Deep Molecular Response in Chronic Myeloid Leukemia**

## Recommended Paper Structure for Statistics in Medicine

1. Introduction: methodological problem, not only CML background.
2. Motivating GLW dataset.
3. Proposed model and likelihood.
4. Estimation and computation.
5. Simulation study design.
6. Simulation results.
7. GLW real-data application.
8. Sensitivity analyses.
9. Discussion: method strengths, limits, and generalizability.

## Go/No-Go Decision

| Submission state | Recommendation |
|---|---|
| Current EDA + proposed model + benchmark models | **No for Statistics in Medicine** |
| Full Bayesian fit but no simulation | **Still risky** |
| Full Bayesian fit + simulation + method comparison | **Possible** |
| Full Bayesian fit + simulation + comparison + reusable code | **Reasonable target** |

## Final Recommendation

Do not submit the current version directly to *Statistics in Medicine*. The dataset supports model innovation as a motivating application, but the paper must be upgraded into a methodological statistics paper. If the immediate goal is publication, target **BMC Medical Research Methodology** first. If the goal is specifically *Statistics in Medicine*, the next step should be to add simulation studies and complete the Bayesian joint model fitting.

## Sources Checked

- *Statistics in Medicine* journal page: https://onlinelibrary.wiley.com/journal/10970258
- Prior journal recommendation file: `08_Model/journal_recommendations.md`
- GLW model summary: `08_Model/joint_interval_model_data_summary_real.csv`
- GLW academic contribution evaluation: `08_Model/academic_highlights_novelty_contribution.md`

