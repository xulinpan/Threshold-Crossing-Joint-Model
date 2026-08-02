# Editorial Audit — *Statistics in Medicine* style (rev10)

**Manuscript:** A Bayesian Threshold-Crossing Joint Model for Irregular, Assay-Floor-Censored Molecular Monitoring in Chronic Myeloid Leukemia (`glw_full_manuscript_smmr_rev10.pdf`, 28 pp., ~12,300 words, 8 figures, 11 tables, 23 refs)

**Audit date:** 28 July 2026 · **Basis:** the current source manuscript.

**Recommendation: Major Revision — not sendable to referees in its present state.** One section contains an unfilled placeholder where results should be, and the paper's central empirical claim is asserted in the Abstract and Discussion on the strength of that missing evidence. The framework itself is in good shape — the restructuring around the threshold-crossing formulation works, the unified single-program section is a genuine advance, and the calibrated reporting of what the data do and do not support is better than most submissions. But the manuscript cannot go out until §5.3 is completed and three internal contradictions are resolved.

---

## A. Blocking

### A1. §5.3 "Robustness to misspecification" contains no results — only a placeholder
Page 20 reads, verbatim:

> *"[Results to be inserted from the completed runs: bias, credible-interval coverage, calibration slope, and computational-failure rate for each scenario, with Monte Carlo standard errors.]"*

The section describes five departures (misspecified trajectory ×2, threshold at −5.0, heterogeneous floors, heavy-tailed *t*₃ errors, informative monitoring) and then reports nothing. A referee reaching page 20 will stop there.

### A2. The Abstract and Conclusions assert a claim that A1 leaves unsupported
Three separate places state that the association is stable **"across the censoring point, prior families and all misspecifications examined"** (Abstract line 83; Introduction line 215; Conclusions line 1366). The misspecification evidence does not exist in the manuscript. As written this is an unsupported claim in the most visible position in the paper. Either complete §5.3 or delete the misspecification clause from all three locations.

### A3. §6.1 directly contradicts §4.6 on the Brier comparison
- **§4.6** (p. 14) now says the Kaplan–Meier and interval-only comparisons *"will not bear weight"*, notes their **negative** patient-level calibration slopes (−0.124, −0.345), and states: *"we... **do not claim predictive superiority over them**."*
- **§6.1** (p. 23) says the model *"**outperformed** Kaplan–Meier, interval-only, landmark, and exact-floor alternatives by Brier score."*

These cannot both stand. §6.1 must be brought into line with §4.6 — the defensible claim is the like-for-like exact-floor contrast (0.082 vs 0.094 interval; 0.116 vs 0.123 patient).

---

## B. Substantive

### B1. Replicate counts are inconsistent across the simulation
§3.8 and §5.2 specify **200** replicates for central scenarios; §5.3 states **50**; §5.5 uses **12**. The 12-replicate multi-state check is explicitly flagged as illustrative (good), but the 50-vs-200 discrepancy is unexplained, and since §5.3 has no results the reader cannot tell whether 50 was planned or achieved. State the achieved replicate count per scenario with Monte Carlo standard errors.

### B2. σ_thr now takes four values across the paper without a single reconciling statement
0.15 (simulation truth, Table 11), 0.18 (multi-state fit, §4.8), 0.27 (unified M3s), 0.78 (unified M2). §5.2's "conditionally identified" paragraph handles the M2/M3/M3s spread well, but Table 11's truth of 0.15 and §4.8's 0.18 sit outside that discussion. Add one sentence tying the simulation truth to the fitted values, or the reader will read the table as a fifth, contradictory estimate.

### B3. The spline result rests on nine events and is described inconsistently
§4.8 reports relapse calibration "improved" (predicted 0.15 vs observed 0.10) — but that is *over*-prediction replacing under-prediction, on 9 confirmed events, with σ_ζ = 4.88 even under a regularising prior. §5.5 correctly calls the multi-state recovery "illustrative rather than definitive." §4.8 should adopt the same register.

### B4. Patient-level calibration remains the weakest reported result
Observed 0.782 vs predicted 0.581, slope 2.82. The paper is honest about this and recalibrates, but the Abstract's clinical sentence ("landmark probabilities... are obtained as monitoring support") should not be read by a clinician as implying usable individual predictions. The existing "internal (apparent)" qualifier is doing a lot of work; consider stating the patient-level miscalibration in the Abstract itself.

---

## C. Compliance and presentation

### C1. Double-blind breach
§6.5 Limitations names *"a single centre in **Kunming, Yunnan Province, China**"* (line 1318), while the title block lists Yunnan University, Kunming. For a double-blind submission this is self-identifying. §3.1 was already genericised to "a single haematology centre in China"; make §6.5 match.

### C2. Ethics statement is a placeholder
Lines 1401–1407 state that the IRB approval "is [to be provided]" and refer to the Title Page. Journals require the approving body and protocol number in the manuscript. This must be completed before submission — it cannot be left to the cover material.

### C3. Colour-vision safety
Figures 5 and 8 encode series by blue/red/green alone. Add redundant shape/linetype so they survive greyscale printing and red–green colour-vision deficiency.

### C4. Figure 7 label collision
"floor n=87 / n=150 / n=300" overprint into an illegible cluster. Regenerate with repelled labels.

---

## D. What is working well (for the record)

- **The restructure succeeds.** Title, Abstract, Introduction and §5 now tell one story: an endpoint defined by a threshold on its own predictor, read as a crossing. The gap paragraph ("The gap is their intersection") is a model of how to position a methods contribution.
- **The unified single-program section (§5) is the strongest part of the paper.** That M1 reproduces the conventional joint model *exactly* is what licenses every other comparison, and it is properly foregrounded.
- **Calibrated honesty throughout** — conditional identifiability of σ_thr, curvature identified only with multi-state transitions, LOO scoped to the longitudinal block with the Pareto-k caveat, dynamic predictions labelled apparent. This is unusually disciplined.
- **Reproducibility** — synthetic data, full pipeline, archived release. Above the norm.

---

## E. Priority order

1. **Complete §5.3** (or remove the misspecification claim from Abstract, Intro, Conclusions) — A1/A2
2. **Reconcile §6.1 with §4.6** on Brier superiority — A3
3. Fix replicate-count inconsistency (50 vs 200) — B1
4. Add one sentence reconciling the four σ_thr values — B2
5. Match §4.8's register on relapse to §5.5's — B3
6. Genericise Kunming/Yunnan in §6.5 — C1
7. Complete the ethics statement — C2
8. Regenerate Figures 5, 7, 8 for colour-vision safety and label legibility — C3/C4

### Bottom line
The intellectual work is done and it is good. What remains is a completion problem, not a conception problem: one empty results section, one stale Discussion sentence that survived the §4.6 revision, and a handful of reporting items. Items 1 and 2 are the difference between a paper that gets reviewed and one that gets returned.

*Prepared against Statistics in Medicine expectations: completeness of reported results, internal consistency of claims, comparator fairness, and reporting compliance.*
