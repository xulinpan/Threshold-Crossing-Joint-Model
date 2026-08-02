# Editorial Audit — *Statistics in Medicine* style (rev8, real source)

**Manuscript:** A Bayesian Joint Longitudinal–Interval and Multi-State Threshold-Crossing Framework for Irregular Molecular Monitoring in CML (`glw_full_manuscript_smmr_rev8.tex`, 26 pp.)

**Audit date:** 26 July 2026 · **Basis:** the real source manuscript with the c_F=−4.5 and serial-correlation sensitivity results now wired in.

**Recommendation: Minor-to-Major Revision.** This is a genuine, well-executed methods paper on a real cohort: reproducible pipeline, exact numbers from the fitted models, real vector figures, a coherent methodological story, and — new in rev8 — two honest sensitivity analyses that strengthen it. The remaining issues are specific and mostly fixable in text, but three are substantive: an internal tension the new serial-correlation result creates with the σ_thr identifiability argument, an estimator mismatch in the prior-sensitivity table, and an unresolved discrepancy in the multi-state threshold width.

---

## Strengths (worth stating for the record)

- The new censoring-point analysis is exactly right: the association α is censoring-robust (−1.27 vs −1.30) while variance/shape parameters move moderately, cleanly separated and honestly reported.
- The serial-correlation analysis is a candid **null result** (σ_total 1.82 ≈ primary 1.81) rather than a claimed fix — good scientific practice.
- Figures and tables are the authoritative vector/generated outputs; the build is clean (26 pp, 0 undefined references).

---

## A. Internal-consistency issues (address first)

### A1. σ_y is now argued to be "measurement error" in §3.3 and "not measurement error" in §4.4
§3.3 rests the identifiability of the threshold width on *"separating the threshold width σ_thr (threshold ambiguity) from the **measurement error** σ_y."* But the new §4.4 shows σ_y ≈ 1.8 is **not** assay measurement error — a serial-correlation decomposition leaves it essentially unchanged, so it is largely unexplained biological/model variance. These two passages now contradict each other. Since the σ_thr↔σ_y separation is one of the paper's two headline identifiability claims, this must be reconciled: either recast σ_y in §3.3 as "residual (observation-level) variability" and argue the separation still holds, or temper the identifiability claim. As written, a referee will read §4.4 and then distrust §3.3.

### A2. The prior-sensitivity table is on a different estimator than the primary analysis
Table 6's caption states the values come from *"the marginal-likelihood re-fit used for this check,"* and its numbers differ from the HMC posterior: τ_b1 = 2.03 (vs 2.20), σ_y = 1.79 (vs 1.81), **α_MRD = −1.18 (vs −1.27)**. So the claim that "estimates were stable across prior families" is demonstrated on marginal-likelihood point estimates, not on the reported HMC posterior, and the α it reports is not the α the paper concludes with. Re-run the prior-sensitivity check under the same HMC estimator (the machinery clearly exists), or state prominently that Table 6 is a marginal-likelihood approximation and reconcile the α discrepancy.

### A3. The multi-state threshold width is unsettled (0.18 vs 0.61) and the HMC fit has 219 divergences
Table 9 reports σ_thr = 0.18 with **asymptotic** 95% intervals from marginal posterior modes. A full HMC threshold fit exists in the repository (`threshold_crossing_model/outputs/`) but gives σ_thr ≈ **0.61** (0.45–0.78) and logged **219 divergent transitions**. So (i) the two parameterizations disagree materially on the very quantity whose identifiability is a headline contribution, and (ii) the HMC fit that could supply proper intervals is not yet clean. The multi-state section currently rests on asymptotic intervals from one parameterization while a divergent HMC fit of a related model disagrees. Resolve which model is authoritative, fit it by HMC at adapt_delta ≥ 0.995 until divergences clear, and report HMC credible intervals in Table 9.

---

## B. Overclaims to temper

### B1. Abstract: "improved on simpler alternatives by Brier score"
The defensible, like-for-like improvement is over the **exact-floor** joint model (patient Brier 0.116 vs 0.123). Against Kaplan–Meier and interval-only the comparison is not clean — those comparators do not beat the intercept-only benchmark and have negative patient calibration slopes (see C1). Narrow the abstract claim to the exact-floor contrast, or qualify.

### B2. Abstract: "calibrated in the large" (unqualified)
Calibration-in-the-large holds at the **interval** level (0.247/0.247) but **not** at the patient level, where cumulative DMR is under-predicted (0.782 vs 0.581) and needs recalibration. The abstract should say "calibrated in the large at the interval level."

### B3. Abstract: "near-nominal credible-interval coverage" from the simulation
The central-scenario truth equals the fitted posterior means, so that coverage validates the sampler, not the model. This is fine if stated; the misspecification grid (§5.3) is the real test and shows σ_y coverage collapsing under t₃ errors — which, given A1, is the policy-relevant number and belongs in the abstract's caveat.

---

## C. Comparators and simulation

### C1. Kaplan–Meier and interval-only have negative patient calibration slopes
Table 7 shows patient cal. int./slope of 2.059/**−0.124** (KM) and 1.271/**−0.345** (interval-only). A negative slope means predictions anti-correlated with outcome — a sign of an implementation problem, not merely a weak model. The text discusses bootstrap and sample-size caveats but does not address the negative slopes. Either diagnose/repair the comparator implementations or state explicitly that they are descriptive references not to be read as competitive predictors.

### C2. Simulation design cannot detect model misspecification in the central scenario
Because truth = fitted posterior means (§5.1/§5 design), the central scenario characterises recovery under correct specification. The paper does run a misspecification grid (§5.3), which is the right instinct; make explicit in §5 that the central-scenario coverage is not evidence the model is correctly specified for the real data.

---

## D. Presentation and reporting

- **Terminology:** "minimal residual disease" → **measurable** residual disease (current nomenclature); use **BCR::ABL1** consistently. Define the control gene and IS conversion factor for log-MRD (reproducibility; interacts with the censoring point).
- **Colour-vision safety:** Figure 5 (dynamic calibration) and Figure 8 (multi-state) encode series by blue/red/green only; add redundant shape/linetype so they survive grayscale and red–green CVD. The figure-generating scripts (`04_Code/R/…`) are the place to fix this.
- **Table 9 caption** describes asymptotic intervals "reported for completeness"; once A3 is resolved, replace with HMC intervals.
- **Prior-sensitivity subsection** (now followed by the new censoring/serial subsection) reads well, but the τ_b1 "≤4%" stability statement is about the marginal-likelihood fit (see A2).
- **β_BM** is identified by few peripheral-blood samples and is correctly flagged as not interpretable; keep that caveat.

---

## E. Smaller points

- The new §4.4 states σ_u has "low effective sample size" — good honesty; note the ESS (~270) explicitly so readers know the decomposition is weakly identified, not that σ_u is precisely 0.43.
- §4.7 "negative log-likelihood fell by 253 for one additional variance parameter" compares nested hierarchical models by raw log-likelihood; report LOO-ELPD ± SE instead.
- Ethics/TRIPOD correctly reside in the submission package (`submission_SMMR/`), appropriate for the anonymized main text.
- Dynamic predictions are correctly labelled apparent, not prospective.

---

## F. Priority order

1. Reconcile the σ_y "measurement error" language in §3.3 with the §4.4 serial-correlation finding (A1).
2. Re-run prior sensitivity under HMC and reconcile α = −1.18 vs −1.27 (A2).
3. Settle the multi-state model (σ_thr 0.18 vs 0.61), refit by HMC without divergences, report HMC intervals (A3).
4. Narrow the abstract's Brier and calibration claims (B1, B2); surface the t₃ coverage caveat (B3).
5. Address the negative comparator slopes (C1).
6. Terminology, log-MRD definition, and colour-vision-safe figures (D).

### Bottom line
rev8 is close to a submittable methods paper and markedly stronger than the pre-refit version. The blocking items are now about *internal consistency* (A1–A3) rather than missing analyses — the serial-correlation null result, while excellent practice, exposes a tension with the identifiability narrative that must be resolved, and the multi-state threshold width needs a single clean HMC fit. Fix those and temper the abstract, and this is a Minor Revision.

*Prepared against Statistics in Medicine expectations: internal consistency, faithful reporting, comparator fairness, and calibrated claims.*
