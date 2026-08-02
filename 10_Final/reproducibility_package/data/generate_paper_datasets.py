import os, zipfile, textwrap
from pathlib import Path
import numpy as np
import pandas as pd

OUT = Path('/mnt/data/cml_paper_datasets')
OUT.mkdir(exist_ok=True)

# --------------------------
# 1. REAL DATASET PROCESSING
# --------------------------
long = pd.read_csv('/mnt/data/longitudinal_model.csv')
surv = pd.read_csv('/mnt/data/survival_model.csv')
covs = pd.read_csv('/mnt/data/patient_covariates_model.csv')

# Stable patient order
patient_ids = sorted(long['patient_id'].unique())
patient_map = {pid: i+1 for i, pid in enumerate(patient_ids)}

long = long.copy()
long['patient_num'] = long['patient_id'].map(patient_map)
long['sample_bm'] = (long['sample_type'].astype(str).str.upper() == 'BM').astype(int)
long['log_bcrabl'] = np.log1p(long['bcr_abl_copy'].astype(float))
long['log_abl'] = np.log1p(long['abl_copy'].astype(float))
long['log_IS'] = np.log1p(long['IS'].astype(float))
long = long.sort_values(['patient_id','t_months'])
long['visit_index'] = long.groupby('patient_id').cumcount()+1
long['prev_t_months'] = long.groupby('patient_id')['t_months'].shift(1).fillna(0)
long['gap_months'] = long['t_months'] - long['prev_t_months']
long['log_gap'] = np.log1p(long['gap_months'].clip(lower=0))
long['dmr_from_log_mrd'] = (long['log_mrd'] <= -4.5).astype(int)

# at-risk intervals: include first interval with DMR, exclude after first DMR
parts=[]
for pid, g in long.groupby('patient_id', sort=True):
    g = g.sort_values('t_months').copy()
    cum_prev = g['dmr_from_log_mrd'].cumsum().shift(1).fillna(0)
    g['at_risk'] = (cum_prev == 0).astype(int)
    g['event_interval'] = ((g['dmr_from_log_mrd'] == 1) & (cum_prev == 0)).astype(int)
    parts.append(g[g['at_risk']==1])
intervals = pd.concat(parts, ignore_index=True)
intervals = intervals[['patient_id','patient_num','visit_index','prev_t_months','t_months','gap_months','log_gap','event_interval']]
intervals = intervals.rename(columns={'prev_t_months':'t_start','t_months':'t_end'})

# survival summary based on processed intervals
summary = intervals.groupby(['patient_id','patient_num']).agg(
    n_intervals=('event_interval','size'),
    dmr_event=('event_interval','max'),
    time_to_dmr_or_censor=('t_end','max'),
    followup_months=('t_end','max')
).reset_index()
# replace time to dmr with first event time if event occurred
first_events = intervals[intervals.event_interval==1].groupby('patient_id')['t_end'].first()
summary['time_to_dmr_or_censor'] = summary.apply(
    lambda r: first_events.get(r['patient_id'], r['time_to_dmr_or_censor']), axis=1
)

# merge covariates for analysis; keep all patients, covariates can be missing
covs2 = covs.copy()
covs2['patient_num'] = covs2['patient_id'].map(patient_map)
real_patient = pd.DataFrame({'patient_id': patient_ids})
real_patient['patient_num'] = real_patient['patient_id'].map(patient_map)
real_patient = real_patient.merge(summary[['patient_id','dmr_event','time_to_dmr_or_censor','followup_months','n_intervals']], on='patient_id', how='left')
real_patient = real_patient.merge(covs2.drop(columns=['patient_num']), on='patient_id', how='left')
real_patient['has_complete_covariates'] = real_patient[['age','duration_years','ph_baseline_pct','sex_male']].notna().all(axis=1).astype(int)

# order real longitudinal columns
real_long_cols = ['patient_id','patient_num','visit_index','t_months','gap_months','sample_type','sample_bm',
                  'bcr_abl_copy','abl_copy','ratio','IS','log_bcrabl','log_abl','log_IS','log_mrd','cmr','dmr','dmr_from_log_mrd']
long[real_long_cols].to_csv(OUT/'real_longitudinal_analysis.csv', index=False)
intervals.to_csv(OUT/'real_interval_survival_analysis.csv', index=False)
real_patient.to_csv(OUT/'real_patient_level_analysis.csv', index=False)

# --------------------------
# 2. SIMULATED DATASET
# --------------------------
rng = np.random.default_rng(20260605)
N = 200
max_months = 60

# True parameters for paper simulation
beta0 = -1.05
beta_time = -0.055
beta_logtime = -0.62
beta_age = 0.006
beta_sex = 0.04
sigma_y = 0.32
sigma_b0 = 0.65
sigma_b1 = 0.025
rho_b = -0.35
sigma_u = 0.55
alpha0 = -2.15
alpha_m = -0.85
alpha_slope = -4.2
alpha_gap = 0.34
gamma_age = 0.006
gamma_sex = 0.12
gamma_ph = -0.004

cov_mat = np.array([[sigma_b0**2, rho_b*sigma_b0*sigma_b1],
                    [rho_b*sigma_b0*sigma_b1, sigma_b1**2]])
random_effects = rng.multivariate_normal([0,0], cov_mat, size=N)
frailty = rng.normal(0, sigma_u, size=N)

sim_patients=[]
sim_long=[]
sim_intervals=[]
for i in range(N):
    pid=f'S{i+1:04d}'
    age = float(np.clip(rng.normal(45, 12), 18, 80))
    sex_male = int(rng.binomial(1,0.55))
    duration_years = float(np.clip(rng.gamma(2.0, 1.6), 0.1, 12))
    ph_baseline = float(np.clip(rng.normal(92, 14), 20, 100))
    ph3 = float(np.clip(ph_baseline * rng.beta(1.3, 3.2), 0, 100))
    ph6 = float(np.clip(ph3 * rng.beta(1.1, 3.8), 0, 100))
    ph12 = float(np.clip(ph6 * rng.beta(0.9, 5.0), 0, 100))
    b0,b1=random_effects[i]
    u=frailty[i]
    t=0.0
    visit=0
    event_seen=False
    event_time=None
    rows_long=[]
    rows_int=[]
    while t < max_months and visit < 14:
        # irregular clinic monitoring, more frequent early, then less frequent
        if t < 12:
            gap = float(rng.gamma(2.0, 1.8))
        else:
            gap = float(rng.gamma(2.0, 3.0))
        gap = float(np.clip(gap, 0.8, 12.0))
        t_start=t
        t=min(t+gap,max_months)
        visit+=1
        # latent process at interval start and end
        def mfun(tt):
            return beta0 + b0 + (beta_time+b1)*tt + beta_logtime*np.log1p(tt) + beta_age*(age-45)/10 + beta_sex*sex_male
        def slopefun(tt):
            return (beta_time+b1) + beta_logtime/(1+tt)
        m_start=mfun(t_start)
        slope_start=slopefun(t_start)
        logit_h = alpha0 + alpha_m*m_start + alpha_slope*slope_start + alpha_gap*np.log1p(gap) + gamma_age*(age-45)/10 + gamma_sex*sex_male + gamma_ph*(ph_baseline-90)/10 + u
        hazard=1/(1+np.exp(-logit_h))
        event_interval=int((not event_seen) and (rng.uniform() < hazard))
        if event_interval:
            event_seen=True
            event_time=t
        latent_log_mrd=mfun(t)
        observed_log_mrd=float(rng.normal(latent_log_mrd,sigma_y))
        # generate copy-number-like variables consistent with MRD
        IS=float(np.clip(10**(observed_log_mrd+4.0), 0, 100)) # rough link; not clinical exactness
        abl_copy=float(np.exp(rng.normal(18.0,0.7)))
        ratio=float(IS/100)
        bcr_abl_copy=int(max(0, np.round(ratio*abl_copy)))
        sample_bm=int(rng.binomial(1,0.75 if t<24 else 0.35))
        rows_long.append({
            'patient_id':pid,'patient_num':i+1,'visit_index':visit,'t_months':t,'gap_months':gap,
            'sample_type':'BM' if sample_bm else 'PB','sample_bm':sample_bm,
            'latent_log_mrd':latent_log_mrd,'log_mrd':observed_log_mrd,
            'bcr_abl_copy':bcr_abl_copy,'abl_copy':abl_copy,'ratio':ratio,'IS':IS,
            'cmr':int(observed_log_mrd <= -5.0),'dmr':int(observed_log_mrd <= -4.5)
        })
        rows_int.append({
            'patient_id':pid,'patient_num':i+1,'visit_index':visit,
            't_start':t_start,'t_end':t,'gap_months':gap,'log_gap':np.log1p(gap),
            'event_interval':event_interval,'true_hazard':hazard,
            'latent_m_start':m_start,'latent_slope_start':slope_start
        })
        if event_seen:
            break
    sim_long.extend(rows_long)
    sim_intervals.extend(rows_int)
    sim_patients.append({
        'patient_id':pid,'patient_num':i+1,'age':age,'sex_male':sex_male,'duration_years':duration_years,
        'ph_baseline_pct':ph_baseline,'ph_3m_pct':ph3,'ph_6m_pct':ph6,'ph_12m_pct':ph12,
        'dmr_event':int(event_seen),'time_to_dmr_or_censor':event_time if event_seen else t,
        'followup_months':t,'n_intervals':visit,'b0_true':b0,'b1_true':b1,'u_true':u
    })

sim_long=pd.DataFrame(sim_long)
sim_intervals=pd.DataFrame(sim_intervals)
sim_patients=pd.DataFrame(sim_patients)

sim_long.to_csv(OUT/'simulated_longitudinal_analysis.csv', index=False)
sim_intervals.to_csv(OUT/'simulated_interval_survival_analysis.csv', index=False)
sim_patients.to_csv(OUT/'simulated_patient_level_analysis.csv', index=False)

params = pd.DataFrame({
    'parameter':['beta0','beta_time','beta_logtime','beta_age','beta_sex','sigma_y','sigma_b0','sigma_b1','rho_b','sigma_u','alpha0','alpha_m','alpha_slope','alpha_gap','gamma_age','gamma_sex','gamma_ph'],
    'true_value':[beta0,beta_time,beta_logtime,beta_age,beta_sex,sigma_y,sigma_b0,sigma_b1,rho_b,sigma_u,alpha0,alpha_m,alpha_slope,alpha_gap,gamma_age,gamma_sex,gamma_ph],
    'description':['longitudinal intercept','linear time effect','log time effect','age effect in longitudinal model','male sex effect in longitudinal model','measurement error SD','random intercept SD','random slope SD','random intercept-slope correlation','survival frailty SD','hazard intercept','association with latent LOG-MRD level','association with latent LOG-MRD slope','gap effect in interval hazard','age effect in hazard','male sex effect in hazard','Philadelphia chromosome effect in hazard']
})
params.to_csv(OUT/'simulation_true_parameters.csv', index=False)

# --------------------------
# 3. DATA DICTIONARY + R SCRIPT
# --------------------------
readme = f"""
# CML joint longitudinal–interval survival datasets for paper

Generated: 2026-06-05

This folder contains analysis-ready datasets for a paper on a Bayesian joint longitudinal–interval survival model for time to deep molecular response (DMR) in imatinib-treated chronic myeloid leukemia.

## Real-data files

### `real_longitudinal_analysis.csv`
Long-format molecular monitoring data.

Main columns:
- `patient_id`: anonymized patient identifier.
- `patient_num`: integer patient index for Stan.
- `visit_index`: within-patient visit number.
- `t_months`: months since imatinib start.
- `gap_months`: interval from previous monitoring time to current visit.
- `sample_type`, `sample_bm`: specimen type; `sample_bm = 1` for bone marrow.
- `bcr_abl_copy`, `abl_copy`, `ratio`, `IS`: molecular measurement variables.
- `log_bcrabl`, `log_abl`, `log_IS`: log-transformed molecular variables.
- `log_mrd`: observed LOG-MRD.
- `dmr_from_log_mrd`: DMR indicator defined by `log_mrd <= -4.5`.

### `real_interval_survival_analysis.csv`
Interval-level survival risk set data. Each row is one at-risk interval ending at a molecular monitoring visit.

Main columns:
- `t_start`, `t_end`: start and end of the interval.
- `gap_months`, `log_gap`: interval length and transformed interval length.
- `event_interval`: equals 1 only for the first interval in which DMR is observed.

### `real_patient_level_analysis.csv`
Patient-level outcome and baseline covariate file.

Main columns:
- `dmr_event`: patient-level indicator of first observed DMR.
- `time_to_dmr_or_censor`: event or censoring time in months.
- `age`, `duration_years`, `ph_baseline_pct`, `ph_3m_pct`, `ph_6m_pct`, `ph_12m_pct`, `ph_18m_pct`, `ph_2y_pct`, `sex_male`: baseline and cytogenetic covariates, where available.
- `has_complete_covariates`: indicator for complete core covariates.

## Simulated-data files

### `simulated_longitudinal_analysis.csv`
Synthetic longitudinal molecular monitoring dataset generated from the proposed joint model.

### `simulated_interval_survival_analysis.csv`
Synthetic interval-level survival dataset generated from the proposed interval hazard model.

### `simulated_patient_level_analysis.csv`
Synthetic patient-level covariates, outcomes, and true random effects.

### `simulation_true_parameters.csv`
True parameter values used for the simulated dataset.

## Real-data summary

- Patients: {len(patient_ids)}
- Longitudinal observations: {len(long)}
- At-risk survival intervals: {len(intervals)}
- DMR event patients: {int(summary['dmr_event'].sum())}
- Censored patients: {int(len(summary)-summary['dmr_event'].sum())}

## Simulated-data summary

- Patients: {N}
- Longitudinal observations: {len(sim_long)}
- At-risk survival intervals: {len(sim_intervals)}
- DMR event patients: {int(sim_patients['dmr_event'].sum())}
- Censored patients: {int(N-sim_patients['dmr_event'].sum())}

## Suggested primary analysis

Use the real molecular-only model on all patients:

`logit(lambda_ij) = alpha0 + alpha_m m_i(t_start) + alpha_slope m_i'(t_start) + alpha_gap log(1 + Delta t_ij) + u_i`.

Use the covariate-adjusted model as a secondary analysis because baseline covariates are available for a subset of patients.
"""
(OUT/'README_data_dictionary.md').write_text(readme, encoding='utf-8')

r_script = r'''
# Build Stan data lists for Bayesian joint longitudinal--interval survival model
library(readr)
library(dplyr)

build_stan_data <- function(folder = "cml_paper_datasets", prefix = "real") {
  long <- read_csv(file.path(folder, paste0(prefix, "_longitudinal_analysis.csv")), show_col_types = FALSE)
  interval <- read_csv(file.path(folder, paste0(prefix, "_interval_survival_analysis.csv")), show_col_types = FALSE)
  patient <- read_csv(file.path(folder, paste0(prefix, "_patient_level_analysis.csv")), show_col_types = FALSE)

  # Primary molecular-only analysis: K = 0
  stan_data <- list(
    N_obs = nrow(long),
    N_pat = length(unique(long$patient_num)),
    N_int = nrow(interval),
    id_obs = as.integer(long$patient_num),
    id_int = as.integer(interval$patient_num),
    y = as.numeric(scale(long$log_mrd)),
    t_obs = as.numeric(scale(long$t_months)),
    log_bcrabl = as.numeric(scale(long$log_bcrabl)),
    log_abl = as.numeric(scale(long$log_abl)),
    sample_bm = as.numeric(long$sample_bm),
    t_start = as.numeric(scale(interval$t_start)),
    t_end = as.numeric(scale(interval$t_end)),
    log_gap = as.numeric(scale(interval$log_gap)),
    event_interval = as.integer(interval$event_interval),
    K = 0,
    Z = matrix(0, nrow = length(unique(long$patient_num)), ncol = 0),
    is_event_patient = as.integer(patient$dmr_event[match(sort(unique(long$patient_num)), patient$patient_num)])
  )
  stan_data
}

# Example:
# stan_data_real <- build_stan_data("/mnt/data/cml_paper_datasets", "real")
# stan_data_sim  <- build_stan_data("/mnt/data/cml_paper_datasets", "simulated")
'''
(OUT/'build_stan_data.R').write_text(r_script, encoding='utf-8')

# Python reproducibility script for the simulation only
py_sim_src = Path('/mnt/data/generate_paper_datasets.py').read_text(encoding='utf-8')
(OUT/'generate_paper_datasets.py').write_text(py_sim_src, encoding='utf-8')

# --------------------------
# 4. ZIP ARCHIVE
# --------------------------
zip_path = Path('/mnt/data/cml_paper_datasets.zip')
with zipfile.ZipFile(zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as z:
    for p in sorted(OUT.iterdir()):
        z.write(p, arcname=f'cml_paper_datasets/{p.name}')

print('Created', zip_path)
for p in sorted(OUT.iterdir()):
    print(p.name, p.stat().st_size)
