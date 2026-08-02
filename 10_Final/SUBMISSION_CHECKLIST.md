# Zenodo + GitHub submission checklist

Package: `10_Final/reproducibility_package/`
Archive: `10_Final/cml_threshold_crossing_v1.0.0.zip` (765 KB, 151 files)

## What was included

| Folder | Contents |
|---|---|
| `R/` | 24 analysis scripts (00–17), including the c_F=−4.5 sensitivity, serial-correlation model, multi-state HMC intervals, unified ladder, and LOO diagnostics |
| `stan/` | 7 Stan models, including `glw_unified_joint.stan` (flag-switched likelihood) |
| `threshold_crossing_model/` | standalone threshold-crossing study (R, Stan, Python) |
| `data/` | **synthetic only** (`simulated_*.csv`, IDs `S0001…`), generating parameters, data dictionary, builders |
| `results/` | aggregate fitted results + generated LaTeX tables (no patient-level rows) |
| `figures/` | 9 manuscript figures (vector PDF) |
| root | `README.md`, `LICENSE` (MIT), `CITATION.cff`, `.zenodo.json`, `.gitignore` |

## What was excluded, and why

Verified absent from both the folder and the archive:

- `03_Data/Raw/glw.csv`, `glw_data.xlsx`, `PH+染色体-Table 1.csv` — **contain a patient-name column (姓名)**
- `03_Data/Processed/private/patient_key.csv` — **maps coded IDs → real patient names**
- `real_longitudinal_analysis.csv`, `real_interval_survival_analysis.csv`, `real_patient_level_analysis.csv`, `real_patient_level_complete_covariates.csv` — de-identified but still the real cohort
- `stan_data_real_*.rds` — serialised real-cohort model input
- 9 patient-level result files (per-record predictions, calibration detail)
- compiled Stan binaries (`*.exe`) — platform-specific; excluded from the archive and blocked by `.gitignore`

Pre-publication audit result: **0 real names, 0 real-cohort files, 0 patient-level result files, synthetic IDs only.**

The scripts `R/00_diagnose_cleaning.R` and `R/01_clean_generate_model_data.R` refer to a
`patient_name` *column* because they process the (unshipped) raw file. These are variable
references, not data. They cannot run without the raw file, which is stated in the README.

---

## Before you publish

- [ ] **Delete the two `.exe` files locally** — the sandbox could not remove them from
      `reproducibility_package/threshold_crossing_model/stan/`. They are already excluded
      from the zip and `.gitignore`, but remove them before pushing:
      `del reproducibility_package\threshold_crossing_model\stan\*.exe`
- [ ] Confirm the **ethics approval number** is added to the manuscript (still a placeholder)
- [ ] Decide whether the **institution/city** should be genericised for double-blind review
      (currently "Kunming, Yunnan Province" appears in the Limitations section)
- [ ] Check `data/README_data_dictionary.md` — it documents the *raw* schema including the
      name column; confirm you are content for that schema to be public

## GitHub

```bash
cd 10_Final/reproducibility_package
git init
git add .
git commit -m "Bayesian threshold-crossing joint model: code and synthetic data (v1.0.0)"
git branch -M main
git remote add origin https://github.com/<user>/<repo>.git
git push -u origin main
git tag -a v1.0.0 -m "v1.0.0 — manuscript submission"
git push origin v1.0.0
```

Verify `git status` shows no `.exe`, no `real_*.csv`, no `private/` before the first push.

## Zenodo

**Option A (recommended) — link GitHub:**
1. Zenodo → Account → **GitHub** → toggle the repository **on**
2. On GitHub, create a **Release** from tag `v1.0.0`
3. Zenodo archives it automatically and mints a DOI; `.zenodo.json` supplies the metadata

**Option B — direct upload:** upload `cml_threshold_crossing_v1.0.0.zip`, then set
upload type *Software*, license *MIT*, and paste the description from `.zenodo.json`.

**Reserve the DOI before final manuscript submission** (Zenodo lets you reserve one
pre-publication) so it can be cited in the Data and Code Availability section, which
currently reads "DOI to be assigned".

## After the DOI is minted

- [ ] Replace "URL provided on acceptance; archived at Zenodo with a DOI to be assigned"
      in the manuscript's Data and Code Availability section with the real DOI and URL
- [ ] Add the DOI badge to `README.md`
- [ ] Add `doi:` and `repository-code:` fields to `CITATION.cff`

---

## Independent verification of the archive (28 Jul 2026)

`cml_threshold_crossing_v1.0.0.zip` (846 KB, 147 files) was extracted and audited
byte-by-byte, not merely by filename:

| Test | Result |
|---|---|
| All **97** real patient names from `patient_key.csv` searched in every file (binary-safe) | **PASS — zero found** |
| Forbidden files (`patient_key`, `glw.csv`, `glw_data.xlsx`, `real_*.csv`, `stan_data_real*`, `anon_*`, crosswalk, `private/`, `*.exe`) | **PASS — none present** |
| Patient-ID prefixes in shipped data | **`S0` only — synthetic** |
| All 8 manuscript figures present | **PASS** |
| Archive integrity (`unzip -t`) | **PASS — no errors** |

### Two issues found and fixed during verification

1. **Three manuscript figures were missing.** `figure_08`, `figure_09` and
   `figure_13` live in `08_Model/`, not `05_Figures/`, so the first build omitted
   them. All eight manuscript figures are now included (plus supporting
   figures 10–12, 14).
2. **`figure_05_raw_missingness.pdf` was removed.** It plots the *raw* schema and
   renders a `patient_name` axis label. No patient data, but it needlessly
   discloses that names were collected; it is not used in the manuscript.

### Non-issues (investigated, benign)

- CJK bytes appear inside the figure PDFs — these are **compression artefacts** in
  the binary streams, not text. All PDF text layers contain **no CJK** and no
  identifiers.
- `R/00_diagnose_cleaning.R` and `R/01_clean_generate_model_data.R` reference a
  `patient_name` **column**; they process the (unshipped) raw file. Variable
  names, not data.
- `R/02_eda_figures_tables_latex.R` contains the literal string `P0001` as an
  *example of the ID format* in generated LaTeX prose.
- `染色体` ("chromosome") appears in `README.md` and the data dictionary as part
  of a filename.
