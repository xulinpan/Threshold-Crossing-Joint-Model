# Submission package — Statistical Methods in Medical Research (SAGE)

**Title:** A Bayesian Threshold-Crossing Joint Model for Irregular,
Assay-Floor-Censored Molecular Monitoring in Chronic Myeloid Leukemia

**Rebuilt:** 1 August 2026 from `02_LaTeX/glw_full_manuscript_smmr_rev10.tex`.

Review model: **double-anonymized**. The manuscript sent to reviewers must contain
no author-identifying information; all identifying content is on the Title Page,
which is *not* sent to reviewers.

## Files in this package

| Order | File | Upload designation (Sage Track) |
|------|------|----------------------------------|
| 1 | `00_cover_letter.pdf` | Cover Letter |
| 2 | `01_title_page.pdf` | Title Page (not sent to reviewers) |
| 3 | `02_manuscript_anonymized.pdf` | Main Document (anonymized) |
| 4 | `figures/Figure1.pdf … Figure8.pdf` | Figure (one file each, in order) |
| 5 | `04_supplementary_material.pdf` | Supplemental File for Review |
| 6 | `TRIPOD_checklist.md` | Supplemental File (reporting guideline) |
| — | `latex_source/` | Tex/LaTeX Suppl Files (upload after Main Document) |

`figures/FigureN.png` are 600-dpi raster alternatives; the `.pdf` figures are
vector and preferred.

## Figure map

Numbered by order of appearance in the manuscript.

| Submission | Source | Content |
|---|---|---|
| Figure1 | `05_Figures/figure_01_log_mrd_trajectories` | Observed log-MRD trajectories |
| Figure2 | `05_Figures/figure_02_time_to_dmr_km` | Time to documented DMR |
| Figure3 | `08_Model/figure_08_posterior_predictive_longitudinal` | Longitudinal posterior predictive check |
| Figure4 | `08_Model/figure_09_dmr_calibration` | DMR calibration |
| Figure5 | `08_Model/figure_13_dynamic_prediction_calibration` | Dynamic-prediction calibration |
| Figure6 | `05_Figures/figure_15_sim_floor_bias` | Assay-floor contrast, censored vs exact |
| Figure7 | `05_Figures/figure_16_sim_calibration` | Simulation calibration-in-the-large |
| Figure8 | `05_Figures/figure_17_multistate` | Multi-state transition probabilities |

## Rebuilding

The anonymized main file is **generated, not hand-edited**. Edit the source
manuscript, then regenerate and compile:

```powershell
cd submission_SMMR\latex_source
xelatex glw_manuscript_anonymized.tex
bibtex  glw_manuscript_anonymized
xelatex glw_manuscript_anonymized.tex
xelatex glw_manuscript_anonymized.tex
copy glw_manuscript_anonymized.pdf ..\02_manuscript_anonymized.pdf

cd ..
xelatex 00_cover_letter.tex
xelatex 01_title_page.tex
```

The generated file differs from the source manuscript in four ways only:
the author block is emptied; `natbib` is switched to superscript numbering per
SMMR house style; graphics and table paths are rebased; and the Statements and
Declarations block is replaced by a pointer to the Title Page, because the
approving ethics committee and the funding award numbers would identify the
authors and the study centre.

## What changed since the 16 July version

The simulation section was rebuilt and several claims corrected:

- **Table 10** regenerated from a paired HMC design (600 fits). The exact-floor
  bias in the random-slope SD is flat in *n* (−1.130, −1.125, −1.123 at
  *n* = 87, 150, 300), so the claim is asymptotic rather than finite-sample.
- **Table 11** extended to three sample sizes with Monte Carlo standard errors.
  Coverage is lower than previously reported because the simulated visit process
  now matches the cohort's spacing instead of being roughly twice as dense.
- **Table 13** reports the misspecification results previously described as
  pending. Coverage of the time slope under response-adaptive monitoring is
  0.44; this is reported and discussed rather than set aside.
- **Table 14** (new) prior sensitivity for the location-scale priors.
- **Table 15** (new) simulation-based calibration, with corrected binning.
- **Table 16** (new) Student-*t* observation model: the two specifications fit
  equally well by leave-one-out yet imply different variance decompositions.
- Figures 6 and 7 regenerated from the corresponding analyses.

## Notes

- Approximate word count: main text ~5,400; ~9,600 counting 8 figures and 16
  tables at 200 words each.
- Figures are legible in grayscale.
- **`response_to_reviewers_point4.md` and `resubmission_cover_note.md` predate
  this revision and must be rewritten before resubmission.**
