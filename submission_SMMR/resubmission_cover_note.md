# Resubmission cover-letter paragraph (Point 4 overhaul)

*Drop-in paragraph for the resubmission cover letter to the Editor.*

---

We are grateful to the reviewer for the detailed critique of our simulation
study, which we have taken seriously and addressed by **redesigning the study
from the ground up**. In the original submission the simulation evaluated a
marginal maximum-*a-posteriori* estimator with Wald intervals — not the
Hamiltonian Monte Carlo (HMC) estimator used in the application — and reported
coverages (e.g., 0.50 for the time slope, 0.12 for a variance component) that,
as the reviewer rightly noted, signalled inferential problems. The revised study
now fits the **same HMC estimator used in the application**, at **200 replicates**
for the central scenario, and reports bias, Bayesian credible-interval coverage,
interval width, and computational-failure rates separately, each with Monte Carlo
standard errors. Under the actual Bayesian estimator the variance components,
residual standard deviation, and the biomarker–event association attain **near-nominal
coverage** (τ_b0 0.94, τ_b1 0.92, α 0.96 at n = 87, improving at n = 150); the
earlier alarming values were an artefact of the Wald/MAP approximation, which we
now state explicitly. We have removed the unqualified word "recovered" from the
abstract and contributions, and we characterise the one genuine finite-sample
effect that remains — attenuation of the latent time slope under heavy
assay-floor censoring — as a limitation that shrinks with sample size rather than
concealing it. We added **simulation-based calibration** (150 prior-predictive
draws), which confirms the posterior is well calibrated (rank statistics uniform
for 10 of 11 parameters), and a **misspecification grid** (wrong trajectory,
heterogeneous floors, heavy-tailed errors, informative monitoring) showing that
the association parameter is robust throughout while the variance components are
sensitive only to gross violations of their assumptions. We believe these changes
convert the simulation from a weakness into support for the paper's claims, and
we thank the reviewer for prompting them.
