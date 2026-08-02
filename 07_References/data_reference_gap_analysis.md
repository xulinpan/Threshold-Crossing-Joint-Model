# Data-Reference Gap Analysis for the GLW CML Study

Prepared: 2026-07-09

Purpose: combine the cleaned GLW dataset and selected peer-reviewed references to identify manuscript-ready research gaps and the strongest contribution of this study.

## Dataset Evidence

- The model-ready GLW dataset contains 87 patients, 495 longitudinal molecular observations, and 275 at-risk survival intervals.
- DMR was observed in 68 patients (78.2%); 19 patients were censored.
- Median follow-up was 39.0 months (IQR 15.8, 65.8).
- Median visits per patient was 5 (IQR 4, 7).
- Visit timing was irregular: median visit gap was 6.0 months (IQR 3.5, 10.9), 50.3% of longitudinal gaps were longer than 6 months, and 20.0% were longer than 12 months.
- Among patients with DMR, the first observed DMR time had median 14.8 months (IQR 9.0, 27.0); 28 patients reached observed DMR by 12 months.
- Complete core covariates were available for 62 of 87 patients (71.3%).
- The cleaned public modeling data use coded patient identifiers only, with direct identifiers isolated in the private patient key.

## Main Research Gap

Current CML guidelines and recent clinical papers clearly define molecular-response milestones, deep molecular response, MR4.5, treatment-free remission, and TKI management. However, the peer-reviewed literature is less directly focused on privacy-preserving real-world datasets where BCR-ABL1/log-MRD trajectories are observed at irregular clinical visit times and the first DMR event is only known to occur within a monitoring interval.

The GLW dataset therefore supports a focused gap: **how to model time to deep molecular response using repeated molecular measurements from irregular real-world monitoring, while accounting for interval-observed event timing, visit gaps, assay-floor behavior, and incomplete baseline covariates.**

## Gap Matrix

| Data finding | Relevant references | What the literature covers | Gap supported by GLW data | Manuscript implication |
|---|---|---|---|---|
| DMR is common but first observed DMR is interval-timed: 68/87 events, median event time 14.8 months among events. | Hochhaus2020ELN; Hughes2014EarlyMolecularResponse; Rasekh2026TFRPredictor; Ono2026MR45 | Response milestones and early molecular response are clinically meaningful; DMR/MR4.5 is linked to TFR eligibility. | Many analyses treat response time as a measured date, but routine monitoring only locates DMR between visits. | Use interval-survival or joint longitudinal-survival modeling rather than exact-time survival alone. |
| Monitoring is irregular: median gap 6.0 months; 50.3% of gaps >6 months; 20.0% >12 months. | Hochhaus2020ELN; Senapati2023CMLManagement; Verweij2026RemoteMonitoring; Parihar2026TFR | Guidelines and real-world studies discuss monitoring schedules and feasibility. | Few studies directly quantify how irregular monitoring affects estimated DMR timing in small real-world cohorts. | Highlight irregular follow-up as a design feature requiring gap-aware modeling. |
| The dataset includes repeated log-MRD/BCR-ABL1 measurements and derived DMR/CMR states. | Kim2026BCRABLqPCR; Hughes2014EarlyMolecularResponse; Atallah2021TFR | Molecular monitoring and response thresholds are well established. | Less attention is paid to modeling the continuous log-MRD trajectory jointly with DMR onset under measurement error or detection-floor behavior. | Treat log-MRD as a longitudinal biomarker that predicts transition into DMR, not only as categorical milestones. |
| Many observations sit at the deep-response floor, with log-MRD values reaching -5. | Kim2026BCRABLqPCR; Ono2026MR45; Atallah2021TFR | Assay performance and MR4.5/DMR concepts are described. | Modeling often treats floor values as exact, even though deep molecular response may be constrained by assay sensitivity or reporting limits. | Discuss possible lower-limit/floor effects and consider sensitivity analyses. |
| Cytogenetic PH+ percentages are available at baseline and landmark time points but are incomplete. | Tokac2026CytogeneticDiversity; Wang2026AdditionalChromosomal; Takahashi2024JSHGuidelines | Cytogenetic abnormalities and Philadelphia chromosome status remain clinically relevant. | There is limited integration of sparse serial cytogenetic response with molecular trajectories in routine datasets. | Use PH+ data descriptively or as secondary covariates, and emphasize sparse-covariate limitations. |
| Complete core covariates are available for 62/87 patients. | Jabbour2025JAMAReview; Senapati2023CMLManagement; Oehler2026TKIManagement | Clinical reviews describe important baseline and treatment factors. | Real-world datasets often have incomplete covariate capture, limiting standard multivariable modeling. | Report complete-case results transparently and plan sensitivity analyses for missing covariates. |
| The data are de-identified into public model-ready files, with a private patient key kept separate. | Verweij2026RemoteMonitoring; Uwizeyimana2026Rwanda | Real-world and remote-monitoring studies show the value of routine CML data. | Privacy-preserving, reproducible model-ready CML molecular datasets are not commonly available in the literature. | Position the cleaned GLW workflow as a reproducible data-engineering contribution. |
| The dataset is small but longitudinally rich: 87 patients and 495 observations. | Wulfsohn1997JointModel; Rizopoulos2016JMbayes; Lovblom2026MixedObservationJointModel; Niu2026QuantileJointModel | Joint models support longitudinal biomarkers and event outcomes. | Recent methods exist, but applied CML examples combining interval-observed DMR and irregular log-MRD monitoring are scarce. | The modeling contribution can be framed as an applied joint/interval method for CML molecular response. |

## Proposed Manuscript Gap Paragraph

Although modern CML guidelines and recent studies emphasize early molecular response, DMR/MR4.5, and treatment-free remission, routine clinical datasets pose a different analytic problem: molecular response is observed through irregular BCR-ABL1 monitoring, and the first attainment of DMR is often only known to occur between two visits. In the GLW cohort, 87 patients contributed 495 longitudinal molecular observations, with a median of 5 visits per patient and a median visit gap of 6.0 months; half of visit gaps exceeded 6 months. DMR was observed in 68 patients, but the timing of first DMR is interval-observed rather than exact. This creates a methodological and clinical gap for models that jointly use continuous log-MRD trajectories and interval-observed DMR onset while handling irregular monitoring, assay-floor behavior, and incomplete covariates. Addressing this gap can improve estimation of molecular-response dynamics in privacy-preserving real-world CML data.

## Recommended Contribution Statement

This study contributes a privacy-preserving, model-ready real-world CML dataset and an analysis framework for deep molecular response that links repeated log-MRD measurements with interval-observed DMR timing under irregular clinical monitoring. The work complements guideline and TFR literature by focusing on the data structure clinicians actually generate in routine follow-up rather than assuming exact response times or fully regular visit schedules.

## Priority References for Gap Framing

- `Hochhaus2020ELN`: guideline anchor for CML treatment and response monitoring.
- `Jabbour2025JAMAReview`: modern clinical overview.
- `Haddad2026Algorithms`: recent treatment-algorithm update.
- `Hughes2014EarlyMolecularResponse`: early molecular-response prognostic anchor.
- `Rasekh2026TFRPredictor`: early molecular response and TFR eligibility.
- `Ono2026MR45`: MR4.5 and TKI discontinuation.
- `Kim2026BCRABLqPCR`: assay-performance context for molecular monitoring.
- `Verweij2026RemoteMonitoring`: real-world/remote monitoring context.
- `Tokac2026CytogeneticDiversity`: Philadelphia chromosome/cytogenetic framing.
- `Wulfsohn1997JointModel`: foundational joint longitudinal-survival modeling.
- `Rizopoulos2016JMbayes`: Bayesian joint-model implementation.
- `Lovblom2026MixedObservationJointModel`: recent mixed-observation joint-model method.

