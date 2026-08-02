# GLW CML model-ready datasets

Generated: 2026-07-08 22:36:34 EDT

The files in this folder are regenerated from `03_Data/Raw/glw.csv` and
`03_Data/Raw/PH+染色体-Table 1.csv`. Public modeling files use coded
patient identifiers only. The reversible name map is isolated under
`03_Data/Processed/private/patient_key.csv` and should not be shared.

## Core derived variables

- `patient_id`: privacy-preserving coded identifier.
- `patient_num`: integer patient index for modeling software.
- `t_months`: months from imatinib start to molecular monitoring date.
- `gap_months`: months since previous visit for that patient.
- `log_mrd`: observed LOG-MRD.
- `dmr` and `dmr_from_log_mrd`: 1 when `log_mrd <= -4.5`.
- `cmr`: 1 when `log_mrd <= -5.0`.
- `event_interval`: first at-risk interval ending in DMR.

## Output files

- `real_longitudinal_analysis.csv`: complete longitudinal molecular visits.
- `real_interval_survival_analysis.csv`: at-risk intervals through first DMR.
- `real_patient_level_analysis.csv`: patient outcomes plus baseline covariates.
- `real_patient_level_complete_covariates.csv`: subset with complete core covariates.
- `data_cleaning_audit.csv`: row counts and exclusion counts.

## Current regenerated counts

- Longitudinal observations: 495
- Patients: 87
- At-risk survival intervals: 275
- DMR event patients: 68
- Censored patients: 19
- Complete core baseline covariates: 62

## Cleaning audit

raw_longitudinal_rows: 504
raw_longitudinal_patients: 89
dropped_missing_or_invalid_core_fields: 9
model_ready_longitudinal_rows: 495
model_ready_longitudinal_patients: 87
raw_patient_table_rows: 87
usable_patient_covariate_rows: 71
patient_level_analysis_rows: 87
complete_core_covariate_rows: 62
interval_survival_rows: 275
dmr_event_patients: 68
censored_patients: 19
