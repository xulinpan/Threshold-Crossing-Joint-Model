# Informative Visit-Process Sensitivity Model

## Purpose

The current CML molecular-monitoring data are irregular: patients with more frequent visits have more opportunities to document DMR, and monitoring frequency may itself depend on latent molecular burden, physician concern, adherence, or clinical workflow. This document specifies an exploratory joint model that adds a recurrent monitoring-visit process to the existing longitudinal and interval-DMR model.

This extension should be treated as an advanced sensitivity analysis, not as the primary model. The cohort has 87 patients, 495 molecular observations, and 275 at-risk intervals; therefore, identifiability and overfitting are central concerns.

## A. Joint Mathematical Model

Let \(t\) denote treatment time in years. The processed data may report time in months for clinical interpretation, but the likelihood should use years, consistent with the current Stan model.

### Latent Longitudinal MRD Trajectory

For patient \(i\), let \(\eta_i(t)\) be the latent true log-MRD trajectory:
\[
\eta_i(t)=
\beta_0+\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2
+\beta_{\mathrm{BM}} I(\mathrm{bone\ marrow}_{ij})
+b_{0i}+b_{1i}\log(1+t).
\]

For the DMR interval and visit-process submodels, the biological latent trajectory should exclude sample-source adjustment:
\[
\eta_i^{\ast}(t)=
\beta_0+\beta_1\log(1+t)+\beta_2\{\log(1+t)\}^2
+b_{0i}+b_{1i}\log(1+t).
\]

Independent random effects are recommended for the exploratory version:
\[
b_{0i}\sim N(0,\tau_0^2), \qquad
b_{1i}\sim N(0,\tau_1^2).
\]

For non-floor observations:
\[
y_{ij}\mid \eta_i(t_{ij}),\sigma_y \sim
N\{\eta_i(t_{ij}),\sigma_y^2\}.
\]

For assay-floor observations at \(c=-5.0\):
\[
\Pr(y_{ij}^{\mathrm{obs}}=c)
=
\Pr(y_{ij}^{\ast}\le c)
=
\Phi\left\{\frac{c-\eta_i(t_{ij})}{\sigma_y}\right\}.
\]

### Interval-Observed First DMR

Let \(L_{ik}\) and \(R_{ik}\) be the start and end of at-risk interval \(k\), with \(\Delta_{ik}=R_{ik}-L_{ik}\), and \(t_{ik}^{mid}=(L_{ik}+R_{ik})/2\). The current interval hazard model can be retained:
\[
h_{ik}^{\mathrm{DMR}}=
\exp\{\gamma_0+\gamma_1\log(1+t_{ik}^{mid})
+\gamma_2\log(1+\Delta_{ik})
+\alpha \eta_i^{\ast}(t_{ik}^{mid})\}.
\]

The probability of first DMR during the interval is:
\[
p_{ik}^{\mathrm{DMR}}=
1-\exp\{-h_{ik}^{\mathrm{DMR}}\Delta_{ik}\}.
\]

The interval contribution is:
\[
L_i^{\mathrm{DMR}}=
\prod_{k=1}^{K_i}
\left(p_{ik}^{\mathrm{DMR}}\right)^{d_{ik}}
\left(1-p_{ik}^{\mathrm{DMR}}\right)^{1-d_{ik}},
\]
where \(d_{ik}=1\) only for the first interval ending in documented DMR.

### Informative Visit Process

Let \(v_{i1},\ldots,v_{iM_i}\) be observed monitoring visit times for patient \(i\) within the visit-process risk window \([E_i,C_i]\). For the most parsimonious sensitivity analysis, use \(E_i=0\) only if monitoring is known to be captured from treatment initiation. If early monitoring history is incomplete, use the first reliably observed monitoring date as \(E_i\) and condition on entry.

The monitoring-visit intensity is:
\[
\lambda_i^{V}(t)=
\exp\{\delta_0+\delta_1[\eta_i^{\ast}(t)+4.5]
+\delta_2\log(1+t)+u_i\},
\]
where \(u_i\) is a patient-level monitoring-intensity random effect:
\[
u_i\sim N(0,\sigma_u^2).
\]

The centering term \(\eta_i^{\ast}(t)+4.5\) makes \(\delta_0\) the log visit rate at the DMR threshold. Under this coding, \(\delta_1>0\) means that patients with higher log-MRD, i.e. worse molecular burden, tend to be monitored more frequently. \(\delta_1<0\) means deeper responders are monitored more frequently.

For the small-cohort exploratory model, \(u_i\) should initially be independent of \(b_{0i}\) and \(b_{1i}\). Correlation between monitoring intensity and the longitudinal random effects is scientifically plausible, but estimating a full random-effect covariance matrix with 87 patients risks non-identifiability.

## B. Likelihood Contribution for Observed Visit Times

Conditional on the latent trajectory and \(u_i\), the visit times are modeled as an inhomogeneous recurrent-event point process. Excluding deterministic baseline measurements, the visit-process likelihood is:
\[
L_i^{V}
=
\left\{
\prod_{\ell=1}^{M_i}
\lambda_i^{V}(v_{i\ell})
\right\}
\exp\left\{
-\int_{E_i}^{C_i}\lambda_i^{V}(s)\,ds
\right\}.
\]

The log-likelihood is:
\[
\ell_i^{V}
=
\sum_{\ell=1}^{M_i}\log \lambda_i^{V}(v_{i\ell})
-
\int_{E_i}^{C_i}\lambda_i^{V}(s)\,ds.
\]

The full joint likelihood is:
\[
L_i =
L_i^{Y}
\times L_i^{\mathrm{DMR}}
\times L_i^{V},
\]
with conditional independence of the observed MRD values, DMR interval process, and monitoring visit process given the latent trajectory and random effects.

In practice, the integral has no convenient closed form because \(\eta_i^{\ast}(t)\) is nonlinear. Use quadrature:
\[
\int_{E_i}^{C_i}\lambda_i^{V}(s)\,ds
\approx
\sum_{q=1}^{Q_i} w_{iq}
\lambda_i^{V}(s_{iq}),
\]
where \(s_{iq}\) and \(w_{iq}\) are quadrature nodes and weights over the observed monitoring-risk window.

## C. Priors for a Small Dataset

Use the current primary model priors for the longitudinal and DMR components:
\[
\beta_0\sim N(-2.5,2^2),\quad
\beta_1\sim N(-1,1^2),\quad
\beta_2\sim N(0,0.5^2),
\]
\[
\beta_{\mathrm{BM}}\sim N(0,1^2),\quad
\sigma_y\sim \mathrm{Exponential}(1),
\]
\[
\gamma_0\sim N(-2,2^2),\quad
\gamma_1,\gamma_2\sim N(0,1^2),\quad
\alpha\sim N(-0.5,0.75^2),
\]
\[
\tau_0,\tau_1\sim \mathrm{Exponential}(1).
\]

For the visit process:
\[
\delta_0\sim N(\log 2,1^2),
\]
where 2 visits/year is a weak prior center consistent with roughly semiannual monitoring. If a different empirical average visit rate is preferred, center \(\delta_0\) at that rate and keep the prior broad.

\[
\delta_1\sim N(0,0.5^2),
\qquad
\delta_2\sim N(0,0.5^2),
\qquad
\sigma_u\sim \mathrm{Exponential}(2).
\]

These priors imply shrinkage toward a non-informative visit process while allowing clinically meaningful departures. Avoid a diffuse prior on \(\delta_1\), because \(\delta_1\), \(u_i\), and the longitudinal random effects can otherwise trade off against each other.

## D. Practical Stan Implementation Strategy

1. Keep the current primary model unchanged for the main manuscript.
2. Build a separate sensitivity Stan program with the same longitudinal and DMR components plus a visit-process block.
3. Exclude deterministic baseline visits from the visit-process likelihood. Count only follow-up monitoring opportunities whose timing plausibly reflects clinical monitoring behavior.
4. Define each patient's monitoring-risk window carefully:
   - If monitoring is complete from treatment initiation, set \(E_i=0\).
   - If early monitoring is incomplete, condition on entry at the first reliably observed visit and integrate from \(E_i\) to \(C_i\).
   - For DMR ascertainment bias, consider using visits up to first documented DMR or censoring as the primary visit-process risk window.
5. Construct quadrature nodes in R and pass them to Stan as `id_quad`, `t_quad`, and `w_quad`.
6. Pass observed visit times as `id_visit` and `t_visit`.
7. Compute:
   \[
   \sum_{\ell}\log\lambda_i^V(v_{i\ell})
   -
   \sum_q w_{iq}\lambda_i^V(s_{iq})
   \]
   in the model block.
8. Generate patient-level visit-process log-likelihood values in `generated quantities` to support patient-level LOO or bootstrap comparison.
9. Do not estimate correlations among \(u_i\), \(b_{0i}\), and \(b_{1i}\) in the first implementation. If the independent-\(u_i\) model is stable and \(\delta_1\) is clearly supported, a correlated random-effect sensitivity model can be explored later.

The companion Stan skeleton is:

`04_Code/Stan/glw_joint_interval_dmr_visit_process_skeleton.stan`

The companion data-preparation script is:

`04_Code/R/11_prepare_visit_process_data.R`

The current prepared visit-process dataset uses the risk window from treatment initiation to first documented DMR or censoring, excludes deterministic \(t=0\) baseline visits, and creates 32 midpoint quadrature nodes per patient. In the current data this produced 224 observed follow-up monitoring visits and 2,784 quadrature nodes. The median number of modeled pre-DMR/censoring visits per patient was 2, which reinforces that this model should be treated as exploratory.

## E. Diagnostics for Informative Monitoring

### Evidence That Visit Timing Is Informative

Do not use a simple significance-test mindset. Instead assess whether:

- the posterior for \(\delta_1\) is concentrated away from zero;
- the posterior for \(\sigma_u\) supports meaningful between-patient monitoring heterogeneity;
- posterior predictive checks reproduce visit counts per patient, median gap, gaps longer than 6 months, and gaps longer than 12 months;
- patient-level time-rescaled residuals for the visit process show no strong systematic lack of fit;
- the visit-process model improves patient-level expected log predictive density for visit timing without degrading longitudinal or DMR calibration.

### Bias Assessment for DMR Probabilities

Compare the current primary model with the visit-process sensitivity model:

- posterior mean and interval for \(\alpha\);
- interval-level DMR Brier score and calibration slope;
- patient-level cumulative DMR Brier score and calibration slope;
- mean absolute change in patient-level predicted DMR probability;
- proportion of patients whose predicted DMR probability changes by more than 0.05 or 0.10;
- calibration stratified by visit frequency or median visit gap.

If \(\delta_1\) is near zero and DMR probabilities change minimally, this supports the practical robustness of the primary model. If \(\delta_1\) is nonzero and DMR probabilities shift materially, report the primary estimates as potentially sensitive to informative monitoring.

### HMC and Identifiability Diagnostics

Check:

- divergent transitions;
- maximum treedepth;
- E-BFMI;
- R-hat and bulk/tail ESS for \(\delta_1\), \(\delta_2\), \(\sigma_u\), \(\alpha\), and random-effect scales;
- posterior pair plots for \(\delta_1\), \(\sigma_u\), \(\tau_0\), \(\tau_1\), and \(\alpha\);
- prior-to-posterior learning for \(\delta_1\);
- sensitivity to the prior scale on \(\delta_1\) and \(\sigma_u\).

## F. Manuscript-Ready Sensitivity-Analysis Section

### Informative Monitoring Sensitivity Analysis

Because molecular monitoring was irregular, patients with more frequent visits had more opportunities to document DMR. In addition, visit frequency may depend on latent molecular response, clinician concern, adherence, or local practice patterns. To evaluate whether the observation process could influence estimated DMR probabilities, we prespecified an exploratory informative-monitoring sensitivity model.

The sensitivity model extended the joint longitudinal--interval framework by adding a recurrent monitoring-visit process. Let \(\eta_i(t)\) denote the latent log-MRD trajectory for patient \(i\). Conditional on \(\eta_i(t)\), monitoring visits were modeled as an inhomogeneous point process with intensity
\[
\lambda_i^{V}(t)=
\exp\{\delta_0+\delta_1[\eta_i(t)+4.5]
+\delta_2\log(1+t)+u_i\},
\]
where \(u_i\sim N(0,\sigma_u^2)\) is a patient-level monitoring-intensity random effect. The latent trajectory was centered at the DMR threshold so that \(\delta_1\) represents the log visit-rate ratio associated with a one-unit higher log-MRD value relative to the DMR threshold. Positive \(\delta_1\) would indicate more frequent monitoring among patients with higher molecular burden, whereas negative \(\delta_1\) would indicate more frequent monitoring among deeper responders.

For observed follow-up visit times \(v_{i1},\ldots,v_{iM_i}\) in the monitoring-risk window \([E_i,C_i]\), the visit-process likelihood contribution was
\[
L_i^{V}=
\left\{\prod_{\ell=1}^{M_i}\lambda_i^{V}(v_{i\ell})\right\}
\exp\left\{-\int_{E_i}^{C_i}\lambda_i^{V}(s)\,ds\right\}.
\]
The cumulative intensity was evaluated by numerical quadrature. Deterministic baseline measurements were excluded from the visit-process likelihood. The longitudinal assay-floor likelihood and the interval-observed DMR likelihood were otherwise unchanged from the primary model.

This model was treated as exploratory because the cohort included only 87 patients and the visit process introduces additional parameters that may be weakly identified. The primary assessment focused on whether \(\delta_1\) was supported by the posterior distribution, whether posterior predictive checks reproduced the observed visit-count and visit-gap distributions, and whether interval-level or patient-level DMR calibration materially changed after accounting for visit intensity. The sensitivity model was not used as the primary clinical prediction model.

## G. Limitations and Overfitting Risk

This extension is scientifically attractive but statistically fragile. The same latent MRD trajectory explains the longitudinal measurements, DMR probability, and visit intensity. As a result, \(\delta_1\) may be confounded with patient random effects, DMR association \(\alpha\), and the monitoring-intensity frailty \(u_i\). With only 87 patients, estimating a fully correlated random-effect structure would likely be overparameterized.

The visit process is also only interpretable if the recorded visit history represents the actual monitoring process. If early visits are missing from the dataset, treating the absence of recorded visits as absence of monitoring would bias the intensity model. In that case, the model should condition on a reliable entry time or be omitted.

Because DMR is only documented at visits, the visit process and the DMR observation process are structurally linked. This sensitivity analysis can evaluate whether primary estimates are robust to informative monitoring, but it cannot fully recover unobserved biological DMR onset without stronger assumptions or external monitoring-schedule information. Therefore, results from this model should be reported as exploratory evidence about observation-process sensitivity rather than as validated clinical decision support.
