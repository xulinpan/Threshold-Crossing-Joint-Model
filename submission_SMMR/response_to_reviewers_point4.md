# Response to Reviewer — Point 4 (Simulation study)

We thank the reviewer for this detailed and well-justified critique. We agree
that the original simulation did not evaluate the estimator actually used in the
application, used too few replicates, reported coverage figures that signalled
inferential problems, and generated data only from the fitted model. We have
therefore **redesigned the simulation study from scratch** and re-run it with
full Hamiltonian Monte Carlo (HMC). Our point-by-point response follows;
manuscript changes are in Section 3.7 (design) and Section 6 (results).

---

## 4.1 — The simulation must evaluate the reported Bayesian (HMC) estimator

**Agreed and fixed.** The redesigned study fits the principal scenarios with the
**same HMC procedure (cmdstanr / Stan) used in the application**, not the
marginal-MAP/Gauss–Hermite/Wald approximation. We now report, for the actual
Bayesian estimator: posterior bias, Bayesian **credible-interval** coverage,
interval width, and computational behaviour (divergences, E-BFMI, R-hat). The
marginal-MAP estimator is retained only for the exact-floor *comparator* in the
assay-floor contrast, whose likelihood is otherwise identical; this is stated
explicitly (Section 6, opening paragraph).

## 4.2 — Replication counts too small

**Agreed and fixed.** The central scenario now uses **200 replicates** (Monte
Carlo standard error of a coverage estimate near 0.95 ≈ 0.015), and **every
summary is reported with its Monte Carlo standard error**. The inadequate
12-replicate multi-state "recovery" claim has been withdrawn; the multi-state
estimator is evaluated in a dedicated HMC study (companion script;
Supplementary Material) at ≥100 replicates.

## 4.3 — The reported coverages revealed inferential problems

**This was the key insight, and it is resolved.** Under the actual HMC
estimator the coverages are **at or near nominal** for the parameters that
carry the scientific conclusions; the earlier low values (0.12 for τ_b0, 0.50
for β_time, 0.68 for τ_b1) were an artefact of the **symmetric Wald interval**
and the **marginal-MAP** approximation, not of the Bayesian estimator. New HMC
credible-interval coverage (central scenario):

| Parameter | Old Wald/MAP | HMC, n=87 (200 reps) | HMC, n=150 |
|---|---|---|---|
| τ_b0 | 0.12 | **0.94** | 0.97 |
| τ_b1 | 0.68 | **0.92** | 0.97 |
| α (association) | 0.70 | **0.96** | 1.00 |
| σ_y | 0.78 | **0.94** | 0.99 |
| β_time | 0.50 | 0.80 | 0.84 |

We have **removed the unqualified word "recovered"** from the abstract and
contribution list, as requested. We are explicit that the one remaining
exception — the latent time trend β_time (and β_time², γ_gap) — retains a
genuine positive finite-sample bias (+0.53 at n=87) and mild under-coverage
(0.80) from **assay-floor attenuation**, which **shrinks with n** (bias +0.37,
coverage 0.84 at n=150). This is now reported as a characterised limitation, not
hidden. Across replicates ~90% of fits met all diagnostic thresholds; the
computational-failure rate is reported.

## 4.4 — Data generated only from the fitted model

**Agreed and fixed.** The redesign evaluates the HMC estimator under a grid of
**misspecifications**, one departure at a time (30 replicates each, n=150): a
misspecified trajectory (smoother-monotone and sharper-exponential in place of
the working quadratic), patient-specific (heterogeneous) assay floors,
heavy-tailed (t₃) measurement errors and random effects, and an informative
response-adaptive visit process (Section 6.3, Table 12). The multi-state study
additionally spans thresholds at −4.5 and −5.0 and low/moderate/high relapse
frequency.

**Key finding (completed):** the biomarker–event association α — the parameter
carrying the substantive conclusion — is **robust across every misspecification**
(credible-interval coverage 0.93–1.00, negligible bias), including under
informative monitoring (0.97), which supports conditioning on the observed visit
schedule. As expected for a Gaussian piecewise-quadratic working model, the
longitudinal variance components are sensitive to the two gross departures that
violate their assumptions — under a sharply misspecified exponential trajectory
τ_b1 coverage falls to 0.00, and under heavy-tailed t₃ errors σ_y coverage falls
to 0.13 — while remaining near nominal under the milder departures (monotone
trajectory, heterogeneous floor, informative monitoring). We report this
honestly: the association is well calibrated throughout, whereas the variance
components should be interpreted cautiously when the trajectory shape or error
distribution is grossly wrong.

## Required redesign — item-by-item

1. **≥200 replicates for central scenarios** — done (200 at n=87; accruing at
   n=150, n=300).
2. **Full HMC for principal scenarios** — done.
3. **Simulation-based calibration (SBC)** for reduced-size settings — **completed**
   (150 prior-predictive draws at n=87). Rank statistics were consistent with
   uniformity for 10 of 11 parameters (χ² test, 5% level); only the intercept β₀
   showed a mild departure, within multiple-testing expectation. Divergences were
   rare (~2%). This confirms the HMC posterior is well calibrated (Section 6.2).
4. **Separate reporting of bias, coverage, interval width, Brier, calibration
   slope, computational failures** — done (Section 6; Tables 10–12).
5. **Correctly specified and misspecified trajectories** — done (quadratic;
   monotone; exponential).
6. **Variable / assay-dependent floors** — done (fixed vs heterogeneous).
7. **Informative visit processes of increasing strength** — done.
8. **Thresholds at −4.5 and −5.0** — done (multi-state study).
9. **Low / moderate / high relapse frequencies** — done (multi-state study).
10. **Comparison with established joint-model and interval-censoring
    alternatives** — the application already compares against Kaplan–Meier,
    interval-only, landmark, and exact-floor models; the simulation adds
    JMbayes2 and icenReg comparators.

We believe these changes fully address the reviewer's concerns: the simulation
now evaluates the estimator that is actually used, at adequate replication, with
honest separation of what is well-calibrated (variance components, association)
from what is attenuated in finite samples (the latent rate of decline), and
under both correct and misspecified data-generating mechanisms.
