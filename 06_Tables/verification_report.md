# Verification checks

- Longitudinal rows: 495
- Patients: 87
- Interval rows: 275
- Patient-level rows: 87
- DMR event patients: 68

## Check results

- longitudinal_required_columns_present: PASS
- longitudinal_complete_required_fields: PASS
- longitudinal_nonnegative_time: PASS
- interval_required_columns_present: PASS
- interval_complete_required_fields: PASS
- interval_time_order: PASS
- one_event_interval_per_patient_max: PASS
- patient_required_columns_present: PASS
- patient_complete_outcome_fields: PASS
- patient_ids_are_coded: PASS
- no_name_columns_in_public_model_files: PASS
- no_cjk_text_in_public_csv_table_latex_outputs: PASS
- figures_exist_and_nonempty: PASS
- tables_exist_and_nonempty: PASS
- latex_pdf_exists_and_nonempty: PASS
- stan_rds_exists_and_nonempty: PASS
