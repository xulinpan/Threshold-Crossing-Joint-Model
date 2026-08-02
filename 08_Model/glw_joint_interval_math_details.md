# Mathematical Details of the GLW Joint Interval-DMR Model

Generated: 2026-07-09

This file summarizes the mathematical model implemented in `04_Code/Stan/glw_joint_interval_dmr.stan`.

## Data Structure

Patients are indexed by \(i=1,\ldots,N\). Patient \(i\) has repeated molecular measurements \(y_{ij}\) and at-risk intervals \((L_{ik},R_{ik}]\) for first DMR.

- \(y_{ij}\): observed log-MRD.
- \(t_{ij}\): years from imatinib start.
- \(x_{ij}=\log(1+t_{ij})\): transformed monitoring time.
- \(z_{ij}\): bone-marrow sample indicator.
- \(\Delta_{ik}=R_{ik}-L_{ik}\): monitoring interval length.
- \(d_{ik}=1\): first DMR occurs in interval \(k\).

## Longitudinal Submodel

\[
\mu_{ij}
=\beta_0+\beta_1x_{ij}+\beta_2x_{ij}^2+\beta_3z_{ij}
+b_{0i}+b_{1i}x_{ij}.
\]

The random effects are

\[
(b_{0i},b_{1i})^T\sim N_2(0,\Sigma_b).
\]

For non-floor observations:

\[
y_{ij}\sim N(\mu_{ij},\sigma_y^2).
\]

For floor observations at \(c=-5\):

\[
P(y_{ij}\leq c)=\Phi((c-\mu_{ij})/\sigma_y).
\]

## Interval-DMR Submodel

At the interval midpoint,

\[
\mu_i(t^{mid}_{ik})
=\beta_0+\beta_1x^{mid}_{ik}+\beta_2(x^{mid}_{ik})^2
+b_{0i}+b_{1i}x^{mid}_{ik}.
\]

The log hazard for first DMR in interval \(k\) is

\[
\eta_{ik}
=\gamma_0+\gamma_1x^{mid}_{ik}
+\gamma_2\log(1+\Delta_{ik})
+\alpha\mu_i(t^{mid}_{ik}).
\]

The cumulative hazard is

\[
H_{ik}=\exp(\eta_{ik})\Delta_{ik}.
\]

The event probability is

\[
P(d_{ik}=1)=1-\exp(-H_{ik}).
\]

The interval likelihood contribution is

\[
\ell^D_{ik}
=d_{ik}\log(1-\exp(-H_{ik}))-(1-d_{ik})H_{ik}.
\]

## Joint Likelihood

\[
L(\theta)
=
\prod_{i=1}^{N}
\left[
\prod_{j=1}^{n_i}L^Y_{ij}(\theta)
\prod_{k=1}^{K_i}L^D_{ik}(\theta)
\right].
\]

The longitudinal and interval-DMR processes are linked through the shared latent trajectory and random effects.

## Priors

\[
\beta_0\sim N(-2.5,2^2),\quad
\beta_1\sim N(-1,1^2),\quad
\beta_2\sim N(0,0.5^2),\quad
\beta_3\sim N(0,1^2),
\]

\[
\sigma_y\sim \mathrm{Exponential}(1),
\]

\[
\gamma_0\sim N(-2,2^2),\quad
\gamma_1\sim N(0,1^2),\quad
\gamma_2\sim N(0,1^2),\quad
\alpha\sim N(-0.5,0.75^2).
\]

Random-effect standard deviations use exponential priors and the correlation matrix uses an LKJ prior.

## Key Interpretation

- \(\alpha<0\): lower log-MRD is associated with higher interval probability of first DMR.
- \(\gamma_2\): evaluates whether longer visit gaps are associated with DMR detection after accounting for interval length.
- \(\tau_0,\tau_1\): quantify patient heterogeneity in baseline molecular burden and response slope.

The full LaTeX version is in `08_Model/glw_joint_interval_math_details.tex`.

