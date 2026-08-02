# GLW Joint Interval-DMR Model

Generated: 2026-07-09

## Model Choice

The primary model should be a joint longitudinal-interval survival model. This follows the data-reference gap analysis: the GLW cohort has repeated log-MRD/BCR-ABL1 measurements, irregular visit gaps, and DMR onset that is observed within monitoring intervals rather than at exact event times.

The model is grounded in the selected references:

- CML response and monitoring guidelines: `Hochhaus2020ELN`, `Jabbour2025JAMAReview`, `Haddad2026Algorithms`, `Senapati2023CMLManagement`.
- Early molecular response, DMR/MR4.5, and TFR: `Hughes2014EarlyMolecularResponse`, `Rasekh2026TFRPredictor`, `Ono2026MR45`, `Atallah2021TFR`.
- Molecular assay and real-world monitoring: `Kim2026BCRABLqPCR`, `Verweij2026RemoteMonitoring`, `Parihar2026TFR`.
- Joint longitudinal-survival methods: `Wulfsohn1997JointModel`, `Rizopoulos2016JMbayes`, `Lovblom2026MixedObservationJointModel`, `Niu2026QuantileJointModel`.

## Data Inputs

Primary processed files:

- `03_Data/Processed/real_longitudinal_analysis.csv`
- `03_Data/Processed/real_interval_survival_analysis.csv`
- `03_Data/Processed/real_patient_level_analysis.csv`

Prepared Stan data:

- `03_Data/Processed/stan_data_real_joint_interval_dmr.rds`

## Longitudinal Submodel

For patient \(i\) at visit \(j\), let \(y_{ij}\) be observed log-MRD and \(t_{ij}\) be years from imatinib start. The latent molecular trajectory is

\[
\mu_{ij} =
\beta_0 + \beta_1 \log(1+t_{ij}) + \beta_2 \{\log(1+t_{ij})\}^2
+ \beta_3 I(\mathrm{bone\ marrow}_{ij}) + b_{0i} + b_{1i}\log(1+t_{ij}).
\]

Non-floor observations use

\[
y_{ij} \sim N(\mu_{ij}, \sigma_y^2).
\]

Because many deep-response values are reported at the assay floor, observations with \(y_{ij}\leq -5\) are treated as left-censored:

\[
P(y_{ij}\leq -5) = \Phi\{(-5 - \mu_{ij})/\sigma_y\}.
\]

## Interval DMR Submodel

For interval \(k\) of patient \(i\), define the interval start \(L_{ik}\), end \(R_{ik}\), and gap \(\Delta_{ik}=R_{ik}-L_{ik}\). DMR onset is modeled with an interval hazard:

\[
h_{ik} =
\exp\{\gamma_0 + \gamma_1 \log(1+t_{ik}^{mid})
+ \gamma_2 \log(1+\Delta_{ik}) + \alpha \mu_i(t_{ik}^{mid})\}.
\]

The probability that first DMR occurs within the interval is

\[
P(d_{ik}=1 \mid \mathrm{at\ risk}) =
1 - \exp(-h_{ik}\Delta_{ik}).
\]

The association parameter \(\alpha\) links the latent log-MRD trajectory to DMR risk. Since lower log-MRD indicates deeper response, a negative \(\alpha\) means lower molecular burden is associated with higher interval probability of DMR.

## Generated Code

- `04_Code/Stan/glw_joint_interval_dmr.stan`: Bayesian joint model.
- `04_Code/R/04_prepare_joint_model_data.R`: creates the Stan data list with raw interval lengths.
- `04_Code/R/05_fit_benchmark_models.R`: fits installed-package benchmark models using `survival` and `nlme`.

## Benchmark Models

Because `rstan` and `cmdstanr` are not installed in this workspace, the generated workflow includes two benchmark models that run now:

1. Interval-censored Weibull survival model for first DMR.
2. Longitudinal mixed model for log-MRD over time with patient-level random effects.

These benchmark models are not the final contribution. They are quality-control models that verify the processed data support the proposed joint model.

