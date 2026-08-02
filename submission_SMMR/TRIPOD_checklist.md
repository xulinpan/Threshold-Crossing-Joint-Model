# TRIPOD Checklist — Prediction Model Development

**Study type:** Development of a prediction model, with internal (apparent)
evaluation only (no external validation). TRIPOD type: **Development only (D).**

**Manuscript:** *A Bayesian Joint Longitudinal–Interval and Multi-State
Threshold-Crossing Framework for Irregular Molecular Monitoring in Chronic
Myeloid Leukemia.*

Section references are to the current manuscript; "Suppl." = Supplementary
Material.

| Item | Section / Topic | Checklist item | Reported? | Location |
|---|---|---|---|---|
| 1 | Title | Identify study as developing a prediction model; target population; outcome | Yes | Title; Abstract |
| 2 | Abstract | Structured summary of objectives, setting, participants, outcome, methods, results, conclusions | Yes | Abstract |
| 3a | Introduction – Background | Rationale, including references to existing models | Yes | §1 |
| 3b | Introduction – Objectives | Specify objectives, including whether development and/or validation | Yes | §1 (contributions) |
| 4a | Methods – Source of data | Study design / data source (e.g., cohort) | Yes | §2 (retrospective monitoring cohort) |
| 4b | Methods – Source of data | Key study dates / follow-up | Yes | §2 (median follow-up 39.0 mo) |
| 5a | Methods – Participants | Study setting, locations | Yes | §2 (single tertiary haematology centre) |
| 5b | Methods – Participants | Eligibility criteria | Yes | §2, §4.1 |
| 5c | Methods – Participants | Details of treatments received | Yes | §1–§2 (TKI therapy) |
| 6a | Methods – Outcome | Definition of outcome predicted | Yes | §2, §3.2 (first documented DMR, log-MRD ≤ −4.5) |
| 6b | Methods – Outcome | Blinded outcome assessment | N/A | Routine assay ascertainment |
| 7a | Methods – Predictors | Predictors used, definitions, measurement | Yes | §3.1 (latent log-MRD trajectory; sample source) |
| 7b | Methods – Predictors | Blinded predictor assessment | N/A | — |
| 8 | Methods – Sample size | How sample size was arrived at | Yes | §4.1 (n=87); §3.7/§6 (simulation of operating characteristics) |
| 9 | Methods – Missing data | How missing data were handled | Yes | §2 (cleaning audit); §4.1 (complete covariates in 62); assay-floor left-censoring §3.1 |
| 10a | Methods – Statistical analysis | How predictors were handled | Yes | §3.1 (log(1+t), left-censoring) |
| 10b | Methods – Statistical analysis | Type of model, building, internal validation | Yes | §3 (Bayesian joint & multi-state HMC); §3.6 (patient-cluster bootstrap, optimism correction) |
| 10c | Methods – Statistical analysis | Performance measures | Yes | §3.6 (calibration-in-the-large, calibration slope, Brier) |
| 10d | Methods – Statistical analysis | Model updating (e.g., recalibration) | Yes | §3.4 (horizon-specific logistic recalibration) |
| 11 | Methods – Risk groups | How risk groups were created | N/A | Continuous predictions; no risk grouping |
| 12 | Methods – Development vs validation | For validation, how it relates to development | N/A | Development only |
| 13a | Results – Participants | Flow of participants; numbers | Yes | §2 (504 records/89 → 495/87); §4.1 |
| 13b | Results – Participants | Characteristics; predictors, outcome | Yes | §4.1; Table 1; Table 2 |
| 13c | Results – Participants | (Validation) comparison to development data | N/A | — |
| 14a | Results – Model development | Numbers of participants and outcome events | Yes | §4.1 (68 DMR, 19 censored) |
| 14b | Results – Model development | Unadjusted association of predictors/outcome | N/A | Model-based analysis |
| 15a | Results – Model specification | Full model (regression coefficients, intercept) | Yes | §5.2; Table 5 (posterior summaries) |
| 15b | Results – Model specification | How to use the model | Yes | §3.4, §5.5 (dynamic prediction) |
| 16 | Results – Model performance | Performance measures (with CIs) | Yes | §5.3–§5.5 (calibration, Brier, comparison); §6 (simulation coverage) |
| 17 | Results – Model updating | Results of any model updating | Yes | §5.5 (recalibrated Brier; Table 9) |
| 18 | Discussion – Limitations | Limitations (e.g., non-representative sample) | Yes | §7 Limitations |
| 19a | Discussion – Interpretation | (Validation) interpretation vs development | N/A | Development only |
| 19b | Discussion – Interpretation | Overall interpretation, considering objectives, other evidence | Yes | §7 (Clinical interpretation; Principal findings) |
| 20 | Discussion – Implications | Potential clinical use and implications | Yes | §7 (monitoring-support, not treatment-decision rule) |
| 21 | Other – Supplementary | Supplementary resources (code, data) | Yes | §Data and Code Availability; Suppl. |
| 22 | Other – Funding | Funding source and role | Yes | Title Page (Funding) |

**Notes.**
- The model is developed and evaluated by *internal (apparent)* performance;
  dynamic predictions are reported as apparent accuracy, and external/temporal
  validation is identified as required future work (§7). This is stated
  explicitly to avoid overstating item 16/19.
- Assay-floor left-censoring (item 9) and interval ascertainment of the outcome
  (item 6a) are handled in the likelihood rather than by imputation; this is a
  methodological feature of the model, described in §3.
