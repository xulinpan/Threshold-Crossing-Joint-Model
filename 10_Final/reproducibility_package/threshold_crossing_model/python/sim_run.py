"""
Checkpointed simulation runner. Each invocation processes pending jobs for
~TIME_BUDGET seconds, appending one CSV row per (scenario, estimator, rep),
then exits. Re-run until 'ALL DONE'.
"""
import os, time, csv, numpy as np
from sim_gen import simulate_cohort, TRUTH
import sim_fit as F

R = 40
TIME_BUDGET = 38.0
OUT = "sim_results.csv"

# scenarios: (name, n, informative)
SCEN = [
    ("ni_n87",  87,  False),
    ("ni_n150", 150, False),
    ("ni_n300", 300, False),
    ("inf_n150",150, True),
]
# which estimators per scenario
EST = {
    "ni_n87":  ["floor","exact","threshold"],
    "ni_n150": ["floor","exact","threshold"],
    "ni_n300": ["floor","exact","threshold"],
    "inf_n150":["floor"],
}
# compute Wald coverage (Hessian) only here:
HESS = {("ni_n150","floor")}

PARAMS_JL = ["b0","b1","b2","log_sigma","log_tau0","log_tau1","g0","g1","g2","alpha"]
PARAMS_TH = ["b0","b1","b2","log_sigma","log_tau0","log_tau1","log_sthr"]
TRUE_NAT = dict(b0=TRUTH['b0'],b1=TRUTH['b1'],b2=TRUTH['b2'],
                log_sigma=np.log(TRUTH['sigma']),log_tau0=np.log(TRUTH['tau0']),
                log_tau1=np.log(TRUTH['tau1']),g0=TRUTH['g0'],g1=TRUTH['g1'],
                g2=TRUTH['g2'],alpha=TRUTH['alpha'])

def logit(p): p=np.clip(p,1e-6,1-1e-6); return np.log(p/(1-p))

def cal_fit(pred, obs):
    """Logistic calibration: obs ~ intercept + slope*logit(pred). Newton."""
    x=logit(pred); X=np.column_stack([np.ones_like(x), x]); b=np.zeros(2)
    for _ in range(50):
        eta=X@b; mu=1/(1+np.exp(-eta)); W=mu*(1-mu)+1e-9
        g=X.T@(obs-mu); H=(X*W[:,None]).T@X
        try: step=np.linalg.solve(H,g)
        except np.linalg.LinAlgError: break
        b=b+step
        if np.max(np.abs(step))<1e-8: break
    return b[0], b[1]

def jobs():
    for name,n,inf in SCEN:
        for est in EST[name]:
            for rep in range(R):
                yield (name,n,inf,est,rep)

def done_set():
    s=set()
    if os.path.exists(OUT):
        with open(OUT) as f:
            for row in csv.DictReader(f):
                s.add((row["scen"],row["est"],int(row["rep"])))
    return s

HEADER=["scen","n","informative","est","rep","conv","brier","cal_int","cal_slope",
        "mean_pred","obs_rate"] + \
        [f"est_{p}" for p in PARAMS_JL] + [f"se_{p}" for p in PARAMS_JL] + ["est_log_sthr","se_log_sthr"]

def run():
    t0=time.time(); did=0
    have=done_set(); first=not os.path.exists(OUT)
    fh=open(OUT,"a",newline=""); w=csv.writer(fh)
    if first: w.writerow(HEADER)
    total=sum(len(EST[s[0]]) for s in SCEN)*R
    for (name,n,inf,est,rep) in jobs():
        if (name,est,rep) in have: continue
        if time.time()-t0>TIME_BUDGET: break
        seed=hash((name,est,rep))%(2**31)
        D=F.pack(simulate_cohort(n, informative=inf, seed=seed))
        res=F.fit(D,est)
        pr,ob=F.predict_dmr(res.x,D,est)
        brier=float(np.mean((pr-ob)**2))
        ci,cs=cal_fit(pr,ob)
        row=dict(scen=name,n=n,informative=int(inf),est=est,rep=rep,
                 conv=int(res.success),brier=brier,cal_int=ci,cal_slope=cs,
                 mean_pred=float(pr.mean()),obs_rate=float(ob.mean()))
        # params
        do_h=(name,est) in HESS
        se=F.wald_se(res.x,D,est) if do_h else np.full(res.x.size,np.nan)
        if est in ("floor","exact"):
            for i,p in enumerate(PARAMS_JL):
                row[f"est_{p}"]=res.x[i]; row[f"se_{p}"]=se[i]
            row["est_log_sthr"]=np.nan; row["se_log_sthr"]=np.nan
        else:
            for i,p in enumerate(PARAMS_TH[:6]):
                row[f"est_{p}"]=res.x[i]; row[f"se_{p}"]=se[i] if do_h else np.nan
            for p in ["g0","g1","g2","alpha"]:
                row[f"est_{p}"]=np.nan; row[f"se_{p}"]=np.nan
            row["est_log_sthr"]=res.x[6]; row["se_log_sthr"]=se[6] if do_h else np.nan
        w.writerow([row.get(h,"") for h in HEADER]); fh.flush()
        did+=1
    fh.close()
    remaining=total-len(done_set())
    print(f"did={did} thiscall  done={len(done_set())}/{total}  remaining={remaining}")
    if remaining<=0: print("ALL DONE")

if __name__=="__main__":
    run()
