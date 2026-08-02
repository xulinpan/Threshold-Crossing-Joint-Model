# How to Run the GLW Models

Generated: 2026-07-09

Run commands from the project root:

```powershell
cd D:\research2026\paper01_glw
```

## RStudio Setup

If you run the scripts from RStudio, first set the project root in the R console:

```r
setwd("D:/research2026/paper01_glw")
```

or set an explicit project-root environment variable:

```r
Sys.setenv(GLW_PROJECT_ROOT = "D:/research2026/paper01_glw")
```

Then run scripts with `source()`:

```r
source("04_Code/R/04_prepare_joint_model_data.R")
source("04_Code/R/05_fit_benchmark_models.R")
```

The scripts now also try to detect the project path from the active RStudio document, but setting `setwd()` or `GLW_PROJECT_ROOT` is the most reliable route.

## 1. Rebuild Data, EDA, Figures, Tables, and Benchmark Models

This is the main reproducible project pipeline:

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\run_all.R
```

This regenerates:

- processed model-ready datasets in `03_Data/Processed`
- EDA figures in `05_Figures`
- tables in `06_Tables`
- LaTeX EDA report in `02_LaTeX`
- joint Stan data in `03_Data/Processed/stan_data_real_joint_interval_dmr.rds`
- benchmark model outputs in `08_Model`

## 2. Prepare Only the Stan Data

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\04_prepare_joint_model_data.R
```

Output:

```text
03_Data/Processed/stan_data_real_joint_interval_dmr.rds
08_Model/joint_interval_model_data_summary_real.csv
```

## 3. Run Only the Benchmark Models

These run now because they use installed R packages `survival` and `nlme`:

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\05_fit_benchmark_models.R
```

Outputs:

```text
08_Model/benchmark_interval_survival_coefficients.csv
08_Model/benchmark_longitudinal_coefficients.csv
08_Model/benchmark_model_fit_summary.csv
08_Model/benchmark_model_summary.txt
```

## 4. Run the Renewed Primary Stan Model

The renewed primary Stan model removes the weakly estimated random-effect
correlation and uses independent patient-specific random intercept and slope
terms. The Stan model file is:

```text
04_Code/Stan/glw_joint_interval_dmr_independent.stan
```

The runner is:

```text
04_Code/R/06_fit_stan_joint_model_independent.R
```

Before running it, install `cmdstanr` and CmdStan. The workspace helper is:

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\install_cmdstan_workspace.R
```

Then run the final renewed model from PowerShell:

```powershell
$env:GLW_MODEL_PREFIX="stan_joint_interval_dmr_independent_renewed"
$env:GLW_STAN_WARMUP="2000"
$env:GLW_STAN_SAMPLING="2000"
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\06_fit_stan_joint_model_independent.R
```

Expected Stan outputs:

```text
08_Model/stan_joint_interval_dmr_independent_renewed_fit.rds
08_Model/stan_joint_interval_dmr_independent_renewed_summary.csv
08_Model/stan_joint_interval_dmr_independent_renewed_draws/
```

After fitting, regenerate posterior predictive checks, calibration summaries,
figures, tables, and diagnostic reports:

```powershell
$env:GLW_MODEL_PREFIX="stan_joint_interval_dmr_independent_renewed"
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\07_sensitivity_calibration_ppc.R
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\08_generate_renewed_model_reporting.R
```

The previous correlated-model outputs are archived in:

```text
08_Model/archive_correlated_model_20260709_renewal/
```

## 5. Run Model Comparison and Internal Validation

After the renewed primary model and exact-floor sensitivity model have been
fit, run:

```powershell
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\09_model_comparison_validation.R
```

This generates:

```text
08_Model/model_comparison_performance.csv
08_Model/table_07_model_comparison_validation.tex
08_Model/figure_10_model_comparison_interval_calibration.pdf
08_Model/figure_11_model_comparison_patient_calibration.pdf
08_Model/model_comparison_internal_validation.tex
08_Model/model_comparison_recommendations.md
08_Model/model_comparison_reviewer_explanation.md
```

## 6. Optional Short Test Run

For a faster smoke test after installing CmdStan, reduce the MCMC iterations:

```powershell
$env:GLW_MODEL_PREFIX="stan_joint_interval_dmr_independent_test"
$env:GLW_STAN_CHAINS="2"
$env:GLW_STAN_PARALLEL_CHAINS="2"
$env:GLW_STAN_WARMUP="250"
$env:GLW_STAN_SAMPLING="250"
& "C:\Program Files\R\R-4.5.0\bin\Rscript.exe" 04_Code\R\06_fit_stan_joint_model_independent.R
```

For the final analysis, use 4 chains with 2000 warmup and 2000 sampling
iterations, then regenerate the checks and tables with the same
`GLW_MODEL_PREFIX`.
