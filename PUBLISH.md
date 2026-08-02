# Publishing to GitHub and Zenodo

Step-by-step for release `v1.0.0`. Read the privacy section first; the rest is
mechanical.

---

## 0. Privacy gate — do this before `git init`

The working tree contains identifiable patient data. `.gitignore` excludes it,
but **verify rather than trust**, because a file committed once stays in the
history even if deleted later.

```powershell
cd D:\research2026\paper01_glw

# nothing sensitive should be listed
git init
git add -A
git status --short | Select-String "private|/Raw/|real_|patient_key|crosswalk|patient_interval_outcome|calibration_patient_detail|landmark_predictions|ppc_longitudinal_observation|visit_process_patient_counts"
```

That command must return **nothing**. If it returns anything, run
`git rm --cached <file>`, add the pattern to `.gitignore`, and re-check.

Excluded by design:

| Excluded | Why |
|---|---|
| `03_Data/Processed/private/` | `patient_key.csv` maps pseudonyms to **patient names** |
| `03_Data/Raw/` | original clinical extracts, contain names |
| `03_Data/Processed/real_*.csv` | analysed patient records |
| `03_Data/Anonymized/` | pseudonymised, resolvable via the crosswalk |
| 9 CSVs in `08_Model/` | per-patient records (age, duration, covariates) |
| posterior draws, checkpoints, `*.exe` | large and regenerable |

If a commit ever does capture one of these, do not simply delete it — rewrite
history with `git filter-repo` or destroy the repository and start again.

---

## 1. GitHub

```powershell
git commit -m "v1.0.0: analysis code, Stan models and simulation study"
git branch -M main
git remote add origin https://github.com/xulinpan/Threshold-Crossing-Joint-Model.git
git push -u origin main
git tag -a v1.0.0 -m "Release accompanying the SMMR submission"
git push origin v1.0.0
```

Repository settings worth enabling: a description and topics matching the
`keywords` in `CITATION.cff`; **Issues** on, so readers can report
reproduction failures; branch protection on `main` if others will contribute.

GitHub renders `CITATION.cff` as a "Cite this repository" button automatically.

---

## 2. Zenodo

1. Sign in to Zenodo with GitHub, go to **Settings → GitHub**, and flip the
   switch for the repository **on**. Do this *before* creating the release —
   Zenodo only captures releases made after the switch is enabled.
2. On GitHub, **Releases → Draft a new release**, choose tag `v1.0.0`, title
   it `v1.0.0`, and paste the "What changed" section from
   `submission_SMMR/03_README_submission_manifest.md` as the notes.
3. Publish. Zenodo mints a DOI within a few minutes and reads `.zenodo.json`
   for the metadata, so authors, licence, keywords and description are already
   correct.
4. Zenodo issues **two** DOIs: a *concept* DOI that always resolves to the
   newest version, and a *version* DOI fixed to `v1.0.0`. Cite the concept DOI
   in the manuscript so it stays valid across revisions.

If you would rather not connect GitHub, upload
`../cml_threshold_crossing_v1.0.0.zip` manually and paste the fields from
`.zenodo.json`.

---

## 3. After the DOI is assigned

Three places still say the DOI is pending:

- `02_LaTeX/glw_full_manuscript_smmr_rev10.tex` — Data and Code Availability
  says *"URL provided on acceptance; archived at Zenodo with a DOI to be
  assigned"*. Replace with the repository URL and the concept DOI.
- `.zenodo.json` — `related_identifiers[0].identifier` is
  `TO-BE-ASSIGNED-ON-ACCEPTANCE`; set it to the article DOI once published, so
  Zenodo links the software to the paper.
- `CITATION.cff` — `repository-code` is set to the repository URL. ✓

---

## 4. Outstanding before release

- **ORCIDs.** `01_title_page.tex` has `ORCID: [to be completed]`. Zenodo
  disambiguates authors by ORCID; without them, "Wang, Chuanming" and
  "Wang, Xueren" may be conflated. Add them to `CITATION.cff` and
  `.zenodo.json`.
- **Corresponding author affiliation is inconsistent.** The title page gives
  Pennsylvania State University / `xxp5066@psu.edu`; the manuscript author
  block gives Yunnan University / `panxl@ynu.edu.cn`, with Penn State
  commented out. The metadata files use Yunnan. Reconcile before submitting.
- **Anonymisation of the raw extracts** is still an open decision — see
  `04_Code/Python/anonymize_cohort.py`, which is written but whose irreversible
  step has not been taken.
