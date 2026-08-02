# Major Methodological Revision Package

Manuscript topic: Bayesian modeling for irregular molecular monitoring and interval-observed deep molecular response in chronic myeloid leukemia.

Prepared: 2026-07-09

## Executive Editorial Recommendation

The manuscript should be revised as a three-level modeling paper rather than as a single clinical prediction model. The revised scientific claim should be:

> This study develops and audits a Bayesian framework for CML molecular monitoring under irregular visits, assay-floor measurements, and interval-observed DMR. The current fitted model is a repaired monitoring-oriented development model. The most coherent methodological upgrade is a latent threshold-crossing formulation, with informative monitoring and dynamic prediction treated as advanced exploratory extensions. External validation remains required before clinical decision use.

This framing resolves the main reviewer risks:

- endpoint circularity is addressed by redefining DMR onset as latent threshold crossing;
- DMR/CMR identity is explained as an assay-floor feature of the retained dataset, not a coding error;
- time scale is standardized to years in the likelihood and months in clinical reporting;
- patient-level underprediction is handled through dynamic prediction and recalibration;
- model comparison and validation are added without overclaiming treatment-decision utility.

## A. Revised Model Hierarchy

### Level 1: Essential Repair

Purpose: maintain a stable primary fitted model suitable for current manuscript reporting.

Model:

- Bayesian joint longitudinal--interval model.
- Longitudinal repeated log-MRD with patient random intercept and random log-time slope.
- Independent random effects, not a correlated random-effect matrix.
- Assay-floor observations left-censored at log-MRD \(\le -5.0\).
- First observed DMR modeled as interval-observed.

Verified current facts:

- 87 patients.
- 495 longitudinal log-MRD observations.
- 275 at-risk DMR intervals.
- 68 patients reached DMR and 19 were censored.
- 239 retained longitudinal observations were floor-coded at \(-5.0\).
- DMR and CMR counts are identical in the processed dataset because there are no retained observations in \((-5.0,-4.5]\).
- The Stan likelihood uses years after converting month variables by dividing by 12.

Status:

- Current repaired independent-random-effects model has been fitted.
- No divergences, maximum R-hat 1.003, minimum bulk ESS 1359, minimum tail ESS 1395, minimum E-BFMI 0.610.
- This level can support the main manuscript if the paper remains an applied monitoring-methods development article.

### Level 2: Strong Methodological Upgrade

Purpose: solve endpoint circularity and make DMR definition model-coherent.

Recommended upgraded primary model:

- Latent threshold-crossing model.
- DMR onset is defined as:
  \[
  T_i^{\mathrm{DMR}}=\inf\{t\ge 0:m_i(t)\le -4.5\}.
  \]
- Observed DMR is interval-observed because monitoring occurs at visits.
- Longitudinal assay values remain exact, left-censored, or optionally interval-censored depending on reporting.
- Compare deterministic and probabilistic threshold-crossing versions.

Status:

- Requires reanalysis.
- A Stan skeleton/design exists, but posterior estimates and validation metrics must not be reported until fitted.

### Level 3: Advanced Extension

Purpose: evaluate whether the real-world observation process affects inference and whether the model can support dynamic monitoring.

Extensions:

- Exploratory informative visit process:
  \[
  \lambda_i^V(t)=\exp\{\delta_0+\delta_1[m_i(t)+4.5]+\delta_2\log(1+t)+u_i\}.
  \]
- Dynamic prediction:
  \[
  \Pr(T_i^{\mathrm{DMR}}\le s+H\mid T_i^{\mathrm{DMR}}>s,\mathcal H_i(s)).
  \]
- Recalibration:
  \[
  \mathrm{logit}(p_{\mathrm{cal}})=a_H+b_H\mathrm{logit}(p_{\mathrm{orig}}).
  \]
- Patient-cluster bootstrap or leave-one-patient-out internal validation.

Status:

- Informative visit-process model requires reanalysis.
- Apparent dynamic prediction and recalibration outputs exist using the fitted posterior, but strict prospective landmark prediction requires reanalysis with landmark-specific posterior updating or leave-one-patient-out refitting.
- Decision-curve analysis should remain exploratory unless a clinically justified action threshold is prespecified.

## B. Full Mathematical Notation

### Time Scale

Let \(t\) denote years since imatinib initiation. The raw and processed clinical files may store months, but all likelihood components use years:
\[
t_{\mathrm{years}}=t_{\mathrm{months}}/12.
\]

Tables and figures may show months for clinical readability.

### Latent Biological MRD Trajectory

Let \(m_i(t)\) denote the latent biological log-MRD trajectory for patient \(i\):
\[
m_i(t)=\beta_0+f(t)+b_{0i}+b_{1i}g(t),
\]
with the current stable choice
\[
f(t)=\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2,
\qquad
g(t)=\log(1+t).
\]

Independent random effects are recommended:
\[
b_{0i}\sim N(0,\tau_0^2), \qquad b_{1i}\sim N(0,\tau_1^2).
\]

The observation-scale mean includes sample source:
\[
\mu_{ij}=m_i(t_{ij})+\beta_{\mathrm{BM}}I(\mathrm{bone\ marrow}_{ij}).
\]

The biological DMR crossing should use \(m_i(t)\), not \(\mu_{ij}\), because sample source is an observation-level measurement shift.

### Longitudinal Observation Model

For exact non-floor observations:
\[
y_{ij}\mid m_i(t_{ij})\sim N(\mu_{ij},\sigma_y^2).
\]

For floor observations at \(c_F=-5.0\):
\[
\Pr(y_{ij}^{\mathrm{obs}}=c_F\mid m_i(t_{ij}))
=
\Pr(y_{ij}^{\ast}\le c_F)
=
\Phi\left\{\frac{c_F-\mu_{ij}}{\sigma_y}\right\}.
\]

If threshold categories are reported, use interval censoring:
\[
\Pr(a<y_{ij}^{\ast}\le b)=
\Phi\left\{\frac{b-\mu_{ij}}{\sigma_y}\right\}
-
\Phi\left\{\frac{a-\mu_{ij}}{\sigma_y}\right\}.
\]

### Level 1 Interval-Hazard DMR Model

For at-risk interval \(k\), define \(L_{ik}\), \(R_{ik}\), \(\Delta_{ik}=R_{ik}-L_{ik}\), and \(t_{ik}^{mid}=(L_{ik}+R_{ik})/2\). The repaired interval-hazard model is:
\[
h_{ik}^{\mathrm{DMR}}=
\exp\{\gamma_0+\gamma_1\log(1+t_{ik}^{mid})
+\gamma_2\log(1+\Delta_{ik})
+\alpha m_i(t_{ik}^{mid})\}.
\]

\[
p_{ik}^{\mathrm{DMR}}=
1-\exp\{-h_{ik}^{\mathrm{DMR}}\Delta_{ik}\}.
\]

This model is useful as a fitted repaired baseline but has potential endpoint circularity because DMR is defined by log-MRD.

### Level 2 Deterministic Threshold-Crossing Model

\[
T_i^{\mathrm{DMR}}=\inf\{t\ge 0:m_i(t)\le c_{\mathrm{DMR}}\},
\qquad c_{\mathrm{DMR}}=-4.5.
\]

For a patient with observed first DMR in \((L_i,R_i]\):
\[
\mathcal L_i^{\mathrm{DMR}}
=
I\{M_i(L_i)>c_{\mathrm{DMR}},\ M_i(R_i)\le c_{\mathrm{DMR}}\},
\]
where
\[
M_i(t)=\min_{0\le u\le t}m_i(u).
\]

For censoring at \(C_i\):
\[
\mathcal L_i^{\mathrm{cens}}=I\{M_i(C_i)>c_{\mathrm{DMR}}\}.
\]

Direct indicator likelihoods are difficult for HMC. A smooth approximation is recommended:
\[
\log \mathcal L_i^{\mathrm{DMR}}
\approx
\log\mathrm{logit}^{-1}\left(\frac{M_i(L_i)-c_{\mathrm{DMR}}}{\epsilon}\right)
+
\log\mathrm{logit}^{-1}\left(\frac{c_{\mathrm{DMR}}-M_i(R_i)}{\epsilon}\right),
\]
where \(\epsilon=0.03\) to \(0.05\) log units.

Status: requires reanalysis.

### Level 2 Probabilistic Threshold-Crossing Model

Let the operational DMR threshold have uncertainty:
\[
C_i\sim N(c_{\mathrm{DMR}},\sigma_{\mathrm{thr}}^2).
\]

For observed DMR in \((L_i,R_i]\):
\[
\mathcal L_i^{\mathrm{DMR}}
=
\Pr\{M_i(R_i)\le C_i<M_i(L_i)\}
=
\Phi\left(\frac{M_i(L_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{thr}}}\right)
-
\Phi\left(\frac{M_i(R_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{thr}}}\right).
\]

For censoring:
\[
\mathcal L_i^{\mathrm{cens}}
=
\Pr\{C_i<M_i(C_i)\}
=
\Phi\left(\frac{M_i(C_i)-c_{\mathrm{DMR}}}{\sigma_{\mathrm{thr}}}\right).
\]

Status: requires reanalysis.

### Level 3 Informative Visit Process

Let \(v_{i1},\dots,v_{iM_i}\) be visit times in \([E_i,C_i]\):
\[
\lambda_i^V(t)
=
\exp\{\delta_0+\delta_1[m_i(t)+4.5]+\delta_2\log(1+t)+u_i\},
\qquad
u_i\sim N(0,\sigma_u^2).
\]

The visit-process likelihood is:
\[
L_i^V=
\left\{\prod_{\ell=1}^{M_i}\lambda_i^V(v_{i\ell})\right\}
\exp\left\{-\int_{E_i}^{C_i}\lambda_i^V(s)\,ds\right\}.
\]

The integral should be evaluated by quadrature:
\[
\int_{E_i}^{C_i}\lambda_i^V(s)\,ds
\approx
\sum_q w_{iq}\lambda_i^V(s_{iq}).
\]

Status: requires reanalysis.

### Dynamic Prediction and Recalibration

At landmark \(s\), horizon \(H\):
\[
\pi_i(s,H)=
\Pr(T_i^{\mathrm{DMR}}\le s+H\mid T_i^{\mathrm{DMR}}>s,\mathcal H_i(s)).
\]

For the interval-hazard approximation:
\[
\pi_i(s,H\mid\theta,b_i)
=
1-\exp\left\{-\int_s^{s+H}h_i(u\mid\theta,b_i)du\right\}.
\]

Recalibration:
\[
\mathrm{logit}\{\pi_{i,\mathrm{cal}}(s,H)\}
=
a_H+b_H\mathrm{logit}\{\pi_{i,\mathrm{orig}}(s,H)\}.
\]

Strict prospective dynamic prediction requires posterior updating of \(b_i\) using only \(\mathcal H_i(s)\). Existing apparent dynamic prediction results use the fitted posterior and are therefore development/prototype estimates.

## C. Recommended Priors

### Longitudinal Model

\[
\beta_0\sim N(-2.5,2^2),\quad
\beta_1\sim N(-1,1^2),\quad
\beta_2\sim N(0,0.5^2).
\]

\[
\beta_{\mathrm{BM}}\sim N(0,1^2),\qquad
\sigma_y\sim \mathrm{Exponential}(1).
\]

\[
\tau_0,\tau_1\sim \mathrm{Exponential}(1).
\]

If a spline time trend is used:
\[
\Delta^2\theta_k\sim N(0,\sigma_f^2),\qquad
\sigma_f\sim \mathrm{Exponential}(2).
\]

### Interval-Hazard Baseline Model

\[
\gamma_0\sim N(-2,2^2),\qquad
\gamma_1,\gamma_2\sim N(0,1^2),
\]
\[
\alpha\sim N(-0.5,0.75^2).
\]

### Threshold-Crossing Model

For deterministic crossing:

- \(\epsilon\) fixed at 0.03 to 0.05 for smooth approximation.
- Do not estimate \(\epsilon\) in the small cohort.

For probabilistic crossing:
\[
\sigma_{\mathrm{thr}}\sim \mathrm{HalfNormal}(0,0.1)
\]
or fixed sensitivity values such as 0.03, 0.05, and 0.10 log units.

Status: requires reanalysis.

### Informative Visit Process

\[
\delta_0\sim N(\log 2,1^2),
\quad
\delta_1\sim N(0,0.5^2),
\quad
\delta_2\sim N(0,0.5^2),
\quad
\sigma_u\sim \mathrm{Exponential}(2).
\]

Avoid estimating correlations among \(u_i,b_{0i},b_{1i}\) in the primary visit-process sensitivity model.

Status: requires reanalysis.

## D. Stan Implementation Plan

### Immediate Fitted Model

Use:

- `04_Code/Stan/glw_joint_interval_dmr_independent.stan`
- `04_Code/R/06_fit_stan_joint_model_independent.R`
- `04_Code/R/08_generate_renewed_model_reporting.R`

This is the repaired Level 1 model.

### Threshold-Crossing Upgrade

Use or extend:

- `04_Code/Stan/glw_latent_threshold_crossing_skeleton.stan`

Implementation details:

1. Separate biological trajectory \(m_i(t)\) from observation mean \(\mu_{ij}\).
2. Keep floor observations as left-censored.
3. Compute \(M_i(t)\) on a dense patient-specific grid or under a monotone/spline trajectory.
4. Implement deterministic crossing using smooth logistic approximations.
5. Implement probabilistic crossing using integrated threshold uncertainty.
6. Generate posterior crossing-time summaries and interval probabilities.

Status: requires reanalysis.

### Visit-Process Extension

Use:

- `04_Code/Stan/glw_joint_interval_dmr_visit_process_skeleton.stan`
- `04_Code/R/11_prepare_visit_process_data.R`

Implementation details:

1. Exclude deterministic baseline visits.
2. Use treatment initiation to first DMR/censoring as the default risk window.
3. Pass `id_visit`, `t_visit`, `id_quad`, `t_quad`, and `w_quad`.
4. Keep \(u_i\) independent of longitudinal random effects.
5. Report \(\delta_1\), \(\sigma_u\), visit-count PPC, and changes in DMR predictions.

Status: requires reanalysis.

### Dynamic Prediction

Use:

- `04_Code/R/12_dynamic_prediction_recalibration.R`

For strict dynamic prediction, add a second-stage posterior-updating script that re-estimates or importance-updates \(b_i\) using only MRD history up to \(s\).

Status:

- apparent development analysis exists;
- strict prospective dynamic prediction requires reanalysis.

## E. R Analysis Workflow

Recommended workflow:

1. Clean and verify raw data:
   - `04_Code/R/01_clean_generate_model_data.R`
   - Verify DMR/CMR definitions.
   - Adjudicate the excluded raw record with IS = 0.002 and missing essential fields.
2. Generate figures/tables:
   - `04_Code/R/02_eda_figures_tables_latex.R`
3. Build Stan data:
   - `04_Code/R/04_prepare_joint_model_data.R`
4. Fit repaired Level 1 primary model:
   - `04_Code/R/06_fit_stan_joint_model_independent.R`
5. Generate repaired model reporting:
   - `04_Code/R/08_generate_renewed_model_reporting.R`
6. Run sensitivity/PPC/calibration:
   - `04_Code/R/07_sensitivity_calibration_ppc.R`
7. Compare models:
   - `04_Code/R/09_model_comparison_validation.R`
8. Compare time functions:
   - `04_Code/R/10_time_function_comparison.R`
9. Prepare informative visit-process data:
   - `04_Code/R/11_prepare_visit_process_data.R`
10. Dynamic prediction/recalibration:
   - `04_Code/R/12_dynamic_prediction_recalibration.R`
11. Fit threshold-crossing model:
   - Requires reanalysis.
12. Run simulation study:
   - Requires new simulation script and reanalysis.

All scripts should keep the RStudio-safe project-root finder already used in this repository.

## F. Model-Comparison Framework

Compare:

1. Descriptive Kaplan--Meier:
   - first observed DMR visit time;
   - descriptive only because true onset is interval-observed.
2. Interval timing + visit-gap model:
   - no longitudinal MRD;
   - tests interval-aware structure alone.
3. Landmark MRD model:
   - last observed MRD before 6, 12, or 18 months;
   - clinically familiar but sample-limited.
4. Joint exact-floor model:
   - treats floor values as exact \(-5.0\);
   - assay-floor sensitivity.
5. Repaired joint left-censored model:
   - current Level 1 primary.
6. Latent threshold-crossing model:
   - methodological upgrade;
   - requires reanalysis.

Existing verified model-comparison results:

- Primary joint model interval Brier: 0.082.
- Primary joint model patient Brier: 0.116.
- Joint exact-floor interval Brier: 0.094.
- Interval timing + visit-gap patient Brier: 0.290.
- Descriptive Kaplan--Meier patient Brier: 0.362.
- Landmark MRD patient Brier: 0.267.

Threshold-crossing model comparison metrics require reanalysis.

## G. Calibration and Validation Framework

### Required Metrics

For interval-level, patient-level, and dynamic landmark predictions:

- Brier score;
- calibration intercept;
- calibration slope;
- mean predicted versus observed rate;
- grouped calibration plot;
- bootstrap optimism correction or leave-one-patient-out validation.

### Current Verified Calibration Findings

Level 1 repaired model:

- Interval observed event rate: 0.247.
- Interval mean predicted probability: 0.247.
- Interval Brier score: 0.082.
- Patient observed DMR rate: 0.782.
- Patient mean predicted DMR probability: 0.581 in the renewed check; 0.597 in model-comparison fixed-posterior summary.
- Patient-level underprediction remains clinically important.

Dynamic prediction:

- Apparent/recalibrated analysis exists.
- 6-month horizon: original mean 0.506 versus observed 0.361; recalibrated Brier 0.157; optimism-corrected 0.167.
- 12-month horizon: original mean 0.626 versus observed 0.576; recalibrated Brier 0.161; optimism-corrected 0.172.
- 24-month horizon: original mean 0.721 versus observed 0.828; recalibrated Brier 0.055; optimism-corrected 0.063.

Strict external or temporal validation: requires reanalysis/new data.

Decision-curve analysis:

- Treat as exploratory unless a clinical action threshold is prespecified.

## H. Simulation Study Design

Purpose: determine whether the proposed models recover DMR onset and calibrated DMR probabilities under irregular visits, floor censoring, and informative monitoring.

### Data-Generating Scenarios

1. Non-informative visits:
   - visit times independent of \(m_i(t)\).
2. Informative visits:
   - \(\lambda_i^V(t)\) depends on \(m_i(t)\).
3. Exact values only:
   - no floor censoring.
4. Assay-floor censoring:
   - observations below \(-5.0\) reported as floor.
5. Threshold uncertainty:
   - probabilistic DMR classification around \(-4.5\).
6. Sparse monitoring:
   - larger visit gaps and fewer visits per patient.

### Sample Sizes

Use:

- \(N=87\), matching current cohort;
- \(N=150\), moderate external validation-like cohort;
- \(N=300\), asymptotic benchmark.

### Simulation Steps

1. Simulate \(b_{0i},b_{1i}\) and \(m_i(t)\).
2. Simulate visit process.
3. Simulate noisy log-MRD with floor censoring.
4. Define true DMR onset by threshold crossing.
5. Record observed DMR interval from visit times.
6. Fit:
   - interval-hazard joint model;
   - deterministic threshold-crossing model;
   - probabilistic threshold-crossing model;
   - no-floor-censoring sensitivity;
   - optional informative visit-process model.
7. Evaluate:
   - bias in DMR onset;
   - coverage of parameters;
   - calibration of DMR probability;
   - Brier score;
   - effect of ignoring informative monitoring.

Status: requires reanalysis.

## I. Revised Manuscript Title

Recommended title:

**Bayesian Threshold-Crossing Joint Modeling for Irregular Molecular Monitoring and Interval-Observed Deep Molecular Response in Chronic Myeloid Leukemia**

Alternative conservative title if the threshold-crossing model remains proposed but not fitted:

**A Bayesian Joint Longitudinal--Interval Framework for Irregular Molecular Monitoring and Deep Molecular Response in Chronic Myeloid Leukemia**

## J. Revised Abstract

### Abstract Template

**Background:** Deep molecular response (DMR) is central to chronic myeloid leukemia monitoring, but routine molecular data are often irregularly observed, floor-limited at deep response levels, and linked to response endpoints detected only at clinic visits. These features complicate standard survival, landmark, and clinical prediction analyses.

**Methods:** We analyzed a real-world CML monitoring cohort of 87 patients with 495 longitudinal log-MRD observations and 275 at-risk intervals for first documented DMR. DMR was defined as log-MRD \(\le -4.5\), and complete molecular response was defined as log-MRD \(\le -5.0\). Treatment time was modeled in years, with clinical summaries reported in months. The repaired primary model was a Bayesian joint longitudinal--interval model with independent patient-specific random intercept and random log-time slope, left-censoring for assay-floor observations, and an interval likelihood for first documented DMR. To address endpoint circularity, we designed a latent threshold-crossing upgrade in which DMR onset is defined as the first time the latent biological MRD trajectory crosses \(-4.5\). Model adequacy was evaluated using posterior predictive checks, calibration summaries, model comparison against simpler alternatives, and dynamic prediction with recalibration.

**Results:** In the processed analysis data, DMR and CMR counts were identical because no retained observations had log-MRD in \((-5.0,-4.5]\); all retained DMR-positive measurements were floor-coded at \(-5.0\). The repaired joint model had acceptable computation, with no divergent transitions, maximum R-hat 1.003, and minimum bulk effective sample size 1359. Interval-level calibration was close in aggregate, whereas patient-level cumulative DMR probability was underpredicted. The primary joint model outperformed Kaplan--Meier, interval-only, landmark, and exact-floor alternatives by Brier score in the current internal comparison. Results for the latent threshold-crossing primary upgrade, informative visit-process extension, and simulation study require reanalysis.

**Conclusions:** The revised framework aligns CML molecular-response modeling with the clinical observation process: repeated MRD measurements, assay-floor reporting, irregular visits, and interval-observed DMR. The repaired model supports monitoring-oriented inference but should not be used as a treatment-decision rule without external validation. A latent threshold-crossing formulation is the preferred methodological upgrade to avoid endpoint circularity.

## K. Revised Methods Section

Use the manuscript-ready LaTeX companion:

`08_Model/major_methodological_revision_sections.tex`

Core method changes:

1. Endpoint audit and time-scale standardization.
2. Assay-floor left-censoring.
3. Repaired Level 1 Bayesian joint model.
4. Level 2 latent threshold-crossing upgrade.
5. Model comparison.
6. Calibration and dynamic prediction.
7. Informative visit-process sensitivity.
8. Simulation plan.

## L. Revised Results Section Template

The revised Results must separate verified current results from items requiring reanalysis.

Verified current results:

- cohort assembly;
- endpoint audit;
- time-scale confirmation;
- assay-floor frequency;
- repaired Level 1 model diagnostics;
- repaired Level 1 posterior summaries;
- model comparison table;
- current dynamic prediction/recalibration prototype.

Requires reanalysis:

- deterministic threshold-crossing posterior results;
- probabilistic threshold-crossing posterior results;
- threshold-crossing model comparison metrics;
- informative visit-process posterior results;
- strict prospective dynamic prediction with landmark-specific updating;
- simulation results;
- external validation.

## M. Revised Discussion Section

Main discussion points:

1. The central methodological contribution is modeling the actual CML monitoring data structure rather than forcing exact event times.
2. The endpoint-circularity concern is real and is best resolved by a latent threshold-crossing model.
3. Identical DMR/CMR descriptive counts are a floor-reporting feature of this dataset.
4. The repaired interval-hazard model is useful as a stable development model but should be framed conservatively.
5. Model comparison supports the value of joint longitudinal--interval modeling over simpler alternatives in this cohort.
6. Dynamic prediction requires recalibration and external validation.
7. Informative monitoring is plausible but weakly identifiable in 87 patients and should remain exploratory.
8. The model is not a treatment-free remission or treatment-decision rule.

## N. Reviewer-Response Style Explanation

**Reviewer concern: The model is circular because DMR is defined by log-MRD but latent log-MRD predicts DMR.**

Response: We agree and revised the framework. The repaired interval-hazard model is retained as a stable development benchmark, but the preferred methodological upgrade defines DMR onset as latent threshold crossing, \(T_i^{\mathrm{DMR}}=\inf\{t:m_i(t)\le -4.5\}\). In this formulation, DMR is no longer a separate event predicted by MRD; it is a functional of the latent MRD trajectory observed only through irregular visits. Fitted threshold-crossing results require reanalysis.

**Reviewer concern: DMR and CMR counts are identical, suggesting a coding error.**

Response: We audited the processed and raw molecular files. DMR and CMR flags match direct threshold recomputation. The retained analysis dataset contains no observations in \((-5.0,-4.5]\); all retained DMR-positive values are floor-coded at \(-5.0\). Thus the identical counts reflect assay reporting and retained-data structure, not an endpoint-coding error. One excluded raw record with IS = 0.002 and missing essential fields should be adjudicated before final submission.

**Reviewer concern: The manuscript mixes months and years.**

Response: We standardized notation. Clinical summaries report months, but all likelihood time variables use years after dividing month variables by 12. Model coefficients are interpreted on the \(\log(1+\mathrm{years})\) scale.

**Reviewer concern: Floor values should not be treated as exact.**

Response: The repaired primary model treats floor observations as left-censored at log-MRD \(\le -5.0\), using a normal CDF likelihood. Exact-floor treatment is now only a sensitivity model.

**Reviewer concern: The model is overparameterized.**

Response: We simplified the primary random-effects structure by using independent patient-specific random intercept and random slope terms, removing the weakly estimated correlation parameter. More complex structures, including informative visit-process random effects, are labeled exploratory.

**Reviewer concern: Patient-level predictions are undercalibrated.**

Response: We added patient-level calibration, dynamic landmark prediction, and logistic recalibration. The current recalibration analysis is internal and developmental; strict prospective dynamic prediction and external validation remain required.

**Reviewer concern: The complex joint model lacks benchmarks.**

Response: We added Kaplan--Meier, interval-only, landmark, exact-floor, and primary joint-model comparisons. In the current internal comparison, the primary joint model had the best interval and patient Brier scores. Threshold-crossing comparison requires reanalysis.

**Reviewer concern: Clinical utility is overclaimed.**

Response: We revised the interpretation. The model is presented as a monitoring-oriented development framework, not as a treatment-decision rule, treatment-free remission eligibility tool, or externally validated clinical prediction model.
