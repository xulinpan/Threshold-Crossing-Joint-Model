# Current Journal Submission Readiness and Target Recommendation

Prepared: 2026-07-09

## Short Decision

Yes, the current results are enough for a peer-reviewed journal if the paper is framed as an applied clinical-statistical modeling study. They are not enough for a top statistical-methods journal or for a strong clinical decision-tool claim.

The best current submission target is:

> **BMC Medical Research Methodology**

The best practical fallback is:

> **PLOS ONE**

The best alternative if the manuscript emphasizes reproducible computational modeling is:

> **BMC Bioinformatics**

## Why the Results Are Publishable

The current analysis has a coherent publishable unit:

- a real CML molecular monitoring cohort;
- 87 patients;
- 495 longitudinal MRD observations;
- 275 at-risk interval-DMR records;
- 68 observed DMR events and 19 censored patients;
- irregular follow-up and interval-observed response;
- 239 assay-floor log-MRD observations;
- fitted Bayesian joint longitudinal-interval model;
- clean MCMC diagnostics for the main fixed effects and event association;
- clinically coherent association between lower latent MRD and higher DMR probability;
- reproducible data-processing, benchmark, Stan, and reporting workflow.

The strongest result is the association between latent MRD and DMR:

| Result | Estimate |
|---|---:|
| `alpha_mrd` | -1.256 |
| 95% posterior interval | -1.851 to -0.809 |
| Effect of 1-unit lower latent log-MRD | about 3.5-fold higher interval DMR hazard component |

This result is clinically meaningful and statistically stable. It supports the paper's central claim that dynamic MRD history has value for monitoring DMR probability.

## Main Limits Reviewers Will Notice

The results are not yet enough for a high-claim clinical prediction paper because:

- the cohort is modest;
- there is no external validation cohort;
- complete clinical covariates are available for only 62 patients;
- model-based DMR probabilities are not externally calibrated;
- no prospective decision-impact analysis is available;
- the random-effect correlation is weakly estimated and should not be interpreted.

These limitations are acceptable for an applied methodology or reproducible modeling paper if clearly stated.

## Recommended Journal Ranking

| Rank | Journal | Current fit | Submit readiness | Best framing |
|---:|---|---|---|---|
| 1 | BMC Medical Research Methodology | Strong | Good after adding sensitivity and calibration checks | Healthcare research methodology for irregular biomarker monitoring |
| 2 | PLOS ONE | Strong fallback | Good if conclusions are conservative | Scientifically valid reproducible applied modeling study |
| 3 | BMC Bioinformatics | Strong if computational workflow is emphasized | Good after code/workflow emphasis | Reproducible statistical model for biological monitoring data |
| 4 | BMC Medical Informatics and Decision Making | Moderate to strong | Needs decision-support framing | Dynamic monitoring and decision-ready MRD modeling |
| 5 | Scientific Reports | Moderate | Needs stronger validation/sensitivity | Broad biomedical applied modeling |
| 6 | Frontiers in Oncology, Hematologic Malignancies | Moderate | Needs stronger hematology discussion | Clinical CML monitoring application |
| 7 | BMC Cancer | Weaker | Risky without external validation | Oncology biomarker/modeling study |
| 8 | Cancers | Optional | Possible, but not best strategic fit | Oncology applied model with reproducible details |
| 9 | Statistical Methods in Medical Research | Stretch | Not ready | Requires simulation and stronger method innovation |
| 10 | Statistics in Medicine | Stretch | Not ready | Requires simulation, method comparison, and general method contribution |

## First-Choice Target: BMC Medical Research Methodology

This is the best fit because the paper's main contribution is methodological for healthcare research: how to model repeated molecular monitoring with irregular visits, assay-floor values, and interval-observed DMR onset.

Recommended title:

> A Joint Longitudinal-Interval Model for Real-World Molecular Monitoring and Deep Molecular Response in Chronic Myeloid Leukemia

Best contribution statement:

> We develop and apply an assay-floor-aware joint longitudinal-interval modeling framework for irregular biomarker monitoring, using real-world CML molecular response data as the motivating clinical application.

Minimum additions before submission:

- posterior predictive check for longitudinal MRD;
- calibration or descriptive validation of interval DMR probabilities;
- sensitivity analysis for assay-floor handling;
- sensitivity analysis for priors or visit-gap term;
- clear data/code availability statement.

## Practical Fallback: PLOS ONE

PLOS ONE is appropriate if the goal is reliable peer-reviewed publication rather than maximizing methodological prestige. The results meet a scientifically valid applied-study standard if the manuscript avoids overclaiming and fully reports the workflow.

Recommended title:

> Reproducible Joint Modeling of Interval-Observed Deep Molecular Response in Real-World CML Monitoring

Best contribution statement:

> This study provides a reproducible clinical-statistical workflow showing that serial MRD trajectories are strongly associated with interval-observed DMR probability in real-world CML monitoring.

## Computational Alternative: BMC Bioinformatics

BMC Bioinformatics becomes attractive if the paper emphasizes the reusable workflow, data processing, Stan model, benchmark models, and reproducible outputs.

Recommended title:

> A Reproducible Joint Modeling Workflow for Longitudinal BCR-ABL1 Monitoring and Deep Molecular Response in CML

To fit this journal better:

- make code availability prominent;
- present model inputs/outputs as a reusable workflow;
- add benchmark comparison against exact-time or baseline-only models;
- include enough implementation detail for reuse.

## Why Not BMC Cancer as First Choice

BMC Cancer is not the best first target because the paper is closer to applied methodology than a validated cancer biomarker study. If framed as a prognostic biomarker paper, reviewers may expect independent validation. The current single-cohort analysis is better suited to methodology or reproducible modeling journals.

## Why Not Statistics in Medicine Now

The current model is clinically meaningful and publishable, but it is not yet a full statistical-methods paper. For Statistics in Medicine or Statistical Methods in Medical Research, the manuscript would need:

- simulation studies;
- bias/coverage/calibration evaluation;
- comparison with exact-time, interval-censored, two-stage, and standard joint models;
- a general method contribution beyond the CML application.

## Recommended Submission Strategy

Submit first to **BMC Medical Research Methodology** after adding a compact sensitivity/calibration section. If the authors want a lower-risk route with broader acceptance criteria, submit to **PLOS ONE**. If the manuscript is rewritten to emphasize reusable computational workflow and model implementation, choose **BMC Bioinformatics**.

## Final Go/No-Go

| Question | Decision |
|---|---|
| Is there enough for a peer-reviewed journal? | Yes |
| Is it enough for a clinical decision-tool paper? | No |
| Is it enough for Statistics in Medicine? | Not yet |
| Best first target | BMC Medical Research Methodology |
| Best pragmatic fallback | PLOS ONE |
| Best computational target | BMC Bioinformatics |

