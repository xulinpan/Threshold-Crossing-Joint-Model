"""Emit LaTeX (booktabs) tables for the simulation section."""
import pandas as pd, numpy as np
A=pd.read_csv("sim_table_bias.csv"); B=pd.read_csv("sim_table_coverage.csv"); C=pd.read_csv("sim_table_calibration.csv")

pretty={'b1':r'$\beta_{\mathrm{time}}$','tau1':r'$\tau_{b1}$','sigma':r'$\sigma_y$',
        'b0':r'$\beta_0$','b2':r'$\beta_{\mathrm{time}^2}$','tau0':r'$\tau_{b0}$',
        'g0':r'$\gamma_0$','g1':r'$\gamma_{\mathrm{time}}$','g2':r'$\gamma_{\mathrm{gap}}$',
        'alpha':r'$\alpha_{\mathrm{MRD}}$'}

# ---- Table 10: bias/RMSE floor vs exact ----
lines=[r"\begin{table}[!htbp]",r"\centering",
 r"\caption{Simulation bias and root mean squared error (RMSE) for key longitudinal parameters under non-informative monitoring, comparing the assay-floor left-censored estimator with the estimator that treats floor values as exact $-5.0$. Truth equals the fitted-model posterior means. 40 replicates per cell.}",
 r"\label{tab:sim-bias}",
 r"\small",
 r"\begin{tabular}{llrrrrr}",r"\toprule",
 r"$n$ & Parameter & True & \multicolumn{2}{c}{Floor left-censored} & \multicolumn{2}{c}{Floor as exact $-5.0$}\\",
 r"\cmidrule(lr){4-5}\cmidrule(lr){6-7}",
 r" & & & Bias & RMSE & Bias & RMSE\\",r"\midrule"]
for n in [87,150,300]:
    for par in ['b1','tau1','sigma']:
        fl=A[(A.n==n)&(A.estimator=='floor')&(A.parameter==par)].iloc[0]
        ex=A[(A.n==n)&(A.estimator=='exact')&(A.parameter==par)].iloc[0]
        lines.append(f"{n} & {pretty[par]} & {fl['true']:.2f} & {fl['bias']:+.2f} & {fl['rmse']:.2f} & {ex['bias']:+.2f} & {ex['rmse']:.2f}\\\\")
    if n!=300: lines.append(r"\addlinespace")
lines += [r"\bottomrule",r"\end{tabular}",r"\end{table}"]
open("table_10_sim_bias.tex","w").write("\n".join(lines))

# ---- Table 11: coverage (floor, n=150) ----
lines=[r"\begin{table}[!htbp]",r"\centering",
 r"\caption{Simulation Wald 95\% interval coverage and mean asymptotic standard error (SE) versus empirical SD for the assay-floor left-censored joint model, non-informative monitoring, $n=150$, 40 replicates. Parameters $\sigma_y,\tau_{b0},\tau_{b1}$ are summarized on the log scale used for estimation.}",
 r"\label{tab:sim-coverage}",r"\small",
 r"\begin{tabular}{lrrrr}",r"\toprule",
 r"Parameter & True & Mean SE & Emp.\ SD & Coverage$_{95}$\\",r"\midrule"]
for _,r in B.iterrows():
    lines.append(f"{pretty.get(r['parameter'],r['parameter'])} & {r['true']:.2f} & {r['mean_SE']:.2f} & {r['empSD']:.2f} & {r['coverage95']:.2f}\\\\")
lines += [r"\bottomrule",r"\end{tabular}",r"\end{table}"]
open("table_11_sim_coverage.tex","w").write("\n".join(lines))

# ---- Table 12: calibration/Brier ----
lines=[r"\begin{table}[!htbp]",r"\centering",
 r"\caption{Simulation patient-level calibration and Brier score for predicted DMR probability. `floor' is the assay-floor left-censored interval model; `threshold' is the probabilistic latent threshold-crossing model. Mean predicted probability, observed documented-DMR rate, and their difference are averaged over 40 replicates.}",
 r"\label{tab:sim-calibration}",r"\small",
 r"\begin{tabular}{llrrrrr}",r"\toprule",
 r"Monitoring & Model ($n$) & Brier & Cal.\ int. & Cal.\ slope & Mean pred. & Obs.\ rate\\",r"\midrule"]
disp=[('ni_n87','floor','Non-informative','floor ($n{=}87$)'),
      ('ni_n150','floor','Non-informative','floor ($n{=}150$)'),
      ('ni_n300','floor','Non-informative','floor ($n{=}300$)'),
      ('ni_n150','threshold','Non-informative','threshold ($n{=}150$)'),
      ('inf_n150','floor','Informative','floor ($n{=}150$)')]
Craw=pd.read_csv("sim_results.csv")
for scen,est,mon,lab in disp:
    d=Craw[(Craw.scen==scen)&(Craw.est==est)]
    lines.append(f"{mon} & {lab} & {d.brier.mean():.3f} & {d.cal_int.mean():+.2f} & {d.cal_slope.mean():.2f} & {d.mean_pred.mean():.2f} & {d.obs_rate.mean():.2f}\\\\")
lines += [r"\bottomrule",r"\end{tabular}",r"\end{table}"]
open("table_12_sim_calibration.tex","w").write("\n".join(lines))
print("wrote table_10_sim_bias.tex table_11_sim_coverage.tex table_12_sim_calibration.tex")
