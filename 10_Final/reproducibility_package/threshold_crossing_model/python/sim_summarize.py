"""Summarize simulation results -> tidy CSVs, LaTeX tables, and figures."""
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sim_gen import TRUTH

df = pd.read_csv("sim_results.csv")

TRUE_NAT = dict(b1=TRUTH['b1'], tau1=TRUTH['tau1'], sigma=TRUTH['sigma'],
                b0=TRUTH['b0'], alpha=TRUTH['alpha'], tau0=TRUTH['tau0'])

def nat(df, name):
    if name in ('b0','b1','b2','g0','g1','g2','alpha'): return df[f'est_{name}']
    if name=='sigma': return np.exp(df['est_log_sigma'])
    if name=='tau0': return np.exp(df['est_log_tau0'])
    if name=='tau1': return np.exp(df['est_log_tau1'])

# ---------- Table A: bias in key longitudinal params, floor vs exact, by n ----------
rowsA=[]
for n,scen in [(87,'ni_n87'),(150,'ni_n150'),(300,'ni_n300')]:
    for est in ['floor','exact']:
        d=df[(df.scen==scen)&(df.est==est)]
        for par,true in [('b1',TRUTH['b1']),('tau1',TRUTH['tau1']),('sigma',TRUTH['sigma'])]:
            v=nat(d,par).values
            rowsA.append(dict(n=n,estimator=est,parameter=par,true=true,
                              mean_est=v.mean(),bias=v.mean()-true,
                              rel_bias_pct=100*(v.mean()-true)/abs(true),
                              empSD=v.std(ddof=1),rmse=np.sqrt(np.mean((v-true)**2))))
A=pd.DataFrame(rowsA); A.to_csv("sim_table_bias.csv",index=False)

# ---------- Table B: coverage (floor, ni_n150) ----------
rowsB=[]
d=df[(df.scen=='ni_n150')&(df.est=='floor')]
for par in ['b0','b1','b2','sigma','tau0','tau1','g0','g1','g2','alpha']:
    if par in ('sigma','tau0','tau1'):
        key={'sigma':'log_sigma','tau0':'log_tau0','tau1':'log_tau1'}[par]
        est=d[f'est_{key}'].values; se=d[f'se_{key}'].values
        true=np.log({'sigma':TRUTH['sigma'],'tau0':TRUTH['tau0'],'tau1':TRUTH['tau1']}[par])
    else:
        est=d[f'est_{par}'].values; se=d[f'se_{par}'].values; true=TRUTH[par if par!='alpha' else 'alpha']
        true={'b0':TRUTH['b0'],'b1':TRUTH['b1'],'b2':TRUTH['b2'],'g0':TRUTH['g0'],
              'g1':TRUTH['g1'],'g2':TRUTH['g2'],'alpha':TRUTH['alpha']}[par]
    ok=np.isfinite(se)&(se>0)
    lo=est-1.96*se; hi=est+1.96*se
    cov=np.mean((true>=lo)&(true<=hi))[()] if ok.any() else np.nan
    cov=np.mean(((true>=lo)&(true<=hi))[ok])
    rowsB.append(dict(parameter=par,true=true,mean_est=est[ok].mean(),
                      mean_SE=se[ok].mean(),empSD=est[ok].std(ddof=1),
                      coverage95=cov,nrep=int(ok.sum())))
B=pd.DataFrame(rowsB); B.to_csv("sim_table_coverage.csv",index=False)

# ---------- Table C: calibration & Brier, floor vs threshold, by n + informative ----------
rowsC=[]
for scen,n,lab in [('ni_n87',87,'non-inf'),('ni_n150',150,'non-inf'),
                    ('ni_n300',300,'non-inf'),('inf_n150',150,'informative')]:
    ests = ['floor','threshold'] if scen!='inf_n150' else ['floor']
    for est in ests:
        d=df[(df.scen==scen)&(df.est==est)]
        if len(d)==0: continue
        rowsC.append(dict(scenario=lab,n=n,monitoring=lab,estimator=est,
                          brier=d.brier.mean(),
                          cal_intercept=d.cal_int.mean(),cal_slope=d.cal_slope.mean(),
                          mean_pred=d.mean_pred.mean(),obs_rate=d.obs_rate.mean(),
                          pred_minus_obs=(d.mean_pred-d.obs_rate).mean()))
C=pd.DataFrame(rowsC); C.to_csv("sim_table_calibration.csv",index=False)

# ---------- Figure 1: random-slope SD bias, floor vs exact, by n ----------
fig,axes=plt.subplots(1,2,figsize=(9,3.8))
for ax,par,true,ttl in [(axes[0],'tau1',TRUTH['tau1'],r'Random-slope SD $\tau_{b1}$'),
                        (axes[1],'b1',TRUTH['b1'],r'Time slope $\beta_{time}$')]:
    xs=[87,150,300];
    for est,mk in [('floor','o'),('exact','s')]:
        means=[]; sds=[]
        for n,scen in zip(xs,['ni_n87','ni_n150','ni_n300']):
            v=nat(df[(df.scen==scen)&(df.est==est)],par).values
            means.append(v.mean()); sds.append(v.std(ddof=1)/np.sqrt(len(v)))
        ax.errorbar(xs,means,yerr=1.96*np.array(sds),marker=mk,capsize=3,
                    label=('floor left-censored' if est=='floor' else 'floor as exact -5.0'))
    ax.axhline(true,ls='--',color='k',lw=1,label='truth')
    ax.set_title(ttl); ax.set_xlabel('sample size n'); ax.set_xticks(xs)
axes[0].legend(fontsize=8,loc='best')
fig.tight_layout(); fig.savefig("figure_15_sim_floor_bias.pdf"); fig.savefig("figure_15_sim_floor_bias.png",dpi=150)

# ---------- Figure 2: calibration (pred vs obs) floor vs threshold + informative ----------
fig,ax=plt.subplots(figsize=(5,4))
order=[('ni_n87','floor','floor n=87'),('ni_n150','floor','floor n=150'),
       ('ni_n300','floor','floor n=300'),('ni_n150','threshold','threshold n=150'),
       ('inf_n150','floor','floor informative n=150')]
for scen,est,lab in order:
    d=df[(df.scen==scen)&(df.est==est)]
    ax.scatter(d.mean_pred.mean(),d.obs_rate.mean(),s=60)
    ax.annotate(lab,(d.mean_pred.mean(),d.obs_rate.mean()),fontsize=7,
                xytext=(4,4),textcoords='offset points')
ax.plot([0,1],[0,1],'k--',lw=1); ax.set_xlim(0.3,0.95); ax.set_ylim(0.3,0.95)
ax.set_xlabel('mean predicted DMR probability'); ax.set_ylabel('observed DMR rate')
ax.set_title('Patient-level DMR calibration-in-the-large')
fig.tight_layout(); fig.savefig("figure_16_sim_calibration.pdf"); fig.savefig("figure_16_sim_calibration.png",dpi=150)

pd.set_option('display.width',200,'display.max_columns',20)
print("=== BIAS (floor vs exact) ===\n", A.round(3).to_string(index=False))
print("\n=== COVERAGE (floor, n=150) ===\n", B.round(3).to_string(index=False))
print("\n=== CALIBRATION / BRIER ===\n", C.round(3).to_string(index=False))
