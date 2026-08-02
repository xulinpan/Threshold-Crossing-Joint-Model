"""
Estimation engine (fully vectorized) for the CML joint-model simulation.
Marginal MAP under weakly-informative priors via 2D Gauss-Hermite quadrature.
Estimators: 'floor' (assay-floor left-censored), 'exact' (floor as exact -5.0),
'threshold' (probabilistic latent threshold-crossing).
"""
import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import logsumexp, log_ndtr, ndtr
from sim_gen import TRUTH, CF, CDMR

LOG2PI = np.log(2*np.pi)

def gh_grid(nq=8):
    x, w = hermegauss(nq)
    w = w/np.sqrt(2*np.pi)
    X0, X1 = np.meshgrid(x, x)
    return X0.ravel(), X1.ravel(), np.log(np.outer(w, w).ravel())
GX0, GX1, LOGGW = gh_grid(8)         # Q=64 nodes
Q = GX0.size

def _logphi(z): return -0.5*(z*z + LOG2PI)

def pack(data):
    """Flatten to contiguous arrays with per-patient segment boundaries."""
    L, I, P = data['long'], data['intr'], data['pat']; n = data['n']
    order = np.argsort(L[:,0], kind='stable'); L = L[order]
    order = np.argsort(I[:,0], kind='stable'); I = I[order]
    lt = np.log1p(L[:,1])
    long = dict(lt=lt, lt2=lt**2, y=L[:,2], cens=L[:,3].astype(bool), pid=L[:,0].astype(int))
    midlog = np.log1p((I[:,1]+I[:,2])/2.0); gaplog = np.log1p(I[:,2]-I[:,1])
    intr = dict(iL=I[:,1], iR=I[:,2], Llog=np.log1p(I[:,1]), Rlog=np.log1p(I[:,2]),
                midlog=midlog, gaplog=gaplog, delta=(I[:,2]-I[:,1]),
                ev=I[:,3].astype(bool), pid=I[:,0].astype(int))
    # reduceat start indices (patients contiguous & complete 0..n-1)
    def starts(pid):
        s = np.searchsorted(pid, np.arange(n)); return s
    long['starts'] = starts(long['pid']); intr['starts'] = starts(intr['pid'])
    pat = dict(Cend=P[:,3], doc=P[:,2].astype(int), onset=P[:,1], n=n)
    return dict(long=long, intr=intr, pat=pat, n=n)

def _seg_sum(mat, starts):
    """Sum (Q,M) over patient segments -> (Q,n) using reduceat on axis 1."""
    return np.add.reduceat(mat, starts, axis=1)

def log_prior(th, kind):
    b0,b1,b2,ls,lt0,lt1 = th[:6]
    sig,t0,t1 = np.exp(ls),np.exp(lt0),np.exp(lt1)
    lp = _logphi((b0+2.5)/2.0)+_logphi((b1+1.0))+_logphi(b2/0.5)
    lp += (-sig)+ls + (-t0)+lt0 + (-t1)+lt1
    if kind in ('floor','exact'):
        g0,g1,g2,al = th[6:10]
        lp += _logphi((g0+2.0)/2.0)+_logphi(g1)+_logphi(g2)+_logphi((al+0.5)/0.75)
    else:
        sthr=np.exp(th[6]); lp += (-sthr)+th[6]
    return lp

def neg_log_post(th, D, kind):
    b0,b1,b2 = th[:3]; sig=np.exp(th[3]); t0=np.exp(th[4]); t1=np.exp(th[5])
    Lg=D['long']; In=D['intr']; n=D['n']
    u0=t0*GX0[:,None]; u1=t1*GX1[:,None]           # (Q,1)
    # longitudinal (Q,M)
    mu = b0 + b1*Lg['lt'][None,:] + b2*Lg['lt2'][None,:] + u0 + u1*Lg['lt'][None,:]
    z=(Lg['y'][None,:]-mu)/sig
    ll_obs=_logphi(z)-np.log(sig)
    if kind=='floor':
        ll_cell=np.where(Lg['cens'][None,:], log_ndtr((CF-mu)/sig), ll_obs)
    else:
        ll_cell=ll_obs
    ll_long=_seg_sum(ll_cell, Lg['starts'])         # (Q,n)
    # event (Q,K)
    if kind in ('floor','exact'):
        g0,g1,g2,al=th[6:10]
        mmid=b0+b1*In['midlog'][None,:]+b2*In['midlog'][None,:]**2+u0+u1*In['midlog'][None,:]
        logh=g0+g1*In['midlog'][None,:]+g2*In['gaplog'][None,:]+al*mmid
        h=np.exp(np.clip(logh,-30,30)); p=np.clip(1-np.exp(-h*In['delta'][None,:]),1e-12,1-1e-12)
        ev_cell=np.where(In['ev'][None,:], np.log(p), np.log1p(-p))
    else:
        sthr=np.exp(th[6])
        mL=b0+b1*In['Llog'][None,:]+b2*In['Llog'][None,:]**2+u0+u1*In['Llog'][None,:]
        mR=b0+b1*In['Rlog'][None,:]+b2*In['Rlog'][None,:]**2+u0+u1*In['Rlog'][None,:]
        PncL=np.clip(ndtr((mL-CDMR)/sthr),1e-12,1-1e-12)
        PncR=np.clip(ndtr((mR-CDMR)/sthr),1e-12,1-1e-12)
        pcr=np.clip(PncL-PncR,1e-12,1-1e-12); surv=np.clip(PncR/PncL,1e-12,1-1e-12)
        ev_cell=np.where(In['ev'][None,:], np.log(pcr), np.log(surv))
    ll_ev=_seg_sum(ev_cell, In['starts'])           # (Q,n)
    tot = ll_long + ll_ev + LOGGW[:,None]           # (Q,n)
    pat_ll = logsumexp(tot, axis=0)                 # (n,)
    return -(pat_ll.sum() + log_prior(th, kind))

def start_theta(kind):
    base=[TRUTH['b0'],TRUTH['b1'],TRUTH['b2'],np.log(TRUTH['sigma']),
          np.log(TRUTH['tau0']),np.log(TRUTH['tau1'])]
    base += [TRUTH['g0'],TRUTH['g1'],TRUTH['g2'],TRUTH['alpha']] if kind in ('floor','exact') else [np.log(0.6)]
    return np.array(base)

def fit(D, kind, jitter=0.0, seed=0):
    x0=start_theta(kind)
    if jitter>0: x0=x0+np.random.default_rng(seed).normal(0,jitter,size=x0.size)
    return minimize(neg_log_post, x0, args=(D,kind), method='L-BFGS-B',
                    options=dict(maxiter=300, maxfun=3000, ftol=1e-7, gtol=1e-5))

def wald_se(th, D, kind):
    g0=neg_log_post(th,D,kind); k=th.size; H=np.zeros((k,k)); e=1e-3
    fp=np.array([neg_log_post(_bump(th,i,e),D,kind) for i in range(k)])
    fm=np.array([neg_log_post(_bump(th,i,-e),D,kind) for i in range(k)])
    for a in range(k):
        for b in range(a,k):
            if a==b:
                H[a,a]=(fp[a]-2*g0+fm[a])/(e*e)
            else:
                fpp=neg_log_post(_bump(_bump(th,a,e),b,e),D,kind)
                fmm=neg_log_post(_bump(_bump(th,a,-e),b,-e),D,kind)
                H[a,b]=H[b,a]=(fpp-fp[a]-fp[b]+2*g0-(-fmm+fm[a]+fm[b]))/(2*e*e)
    try:
        se=np.sqrt(np.clip(np.diag(np.linalg.inv(H)),0,None))
    except np.linalg.LinAlgError:
        se=np.full(k,np.nan)
    return se

def _bump(th,i,e):
    t=th.copy(); t[i]+=e; return t

def predict_dmr(th, D, kind):
    b0,b1,b2=th[:3]; sig=np.exp(th[3]); t0=np.exp(th[4]); t1=np.exp(th[5])
    Lg=D['long']; In=D['intr']; n=D['n']
    u0=t0*GX0[:,None]; u1=t1*GX1[:,None]
    mu=b0+b1*Lg['lt'][None,:]+b2*Lg['lt2'][None,:]+u0+u1*Lg['lt'][None,:]
    z=(Lg['y'][None,:]-mu)/sig; ll_obs=_logphi(z)-np.log(sig)
    if kind=='floor':
        ll_cell=np.where(Lg['cens'][None,:], log_ndtr((CF-mu)/sig), ll_obs)
    else:
        ll_cell=ll_obs
    ll_long=_seg_sum(ll_cell, Lg['starts'])+LOGGW[:,None]   # (Q,n) longitudinal posterior (unnorm)
    W=np.exp(ll_long-ll_long.max(0)); W/=W.sum(0)           # (Q,n)
    if kind in ('floor','exact'):
        g0,g1,g2,al=th[6:10]
        mmid=b0+b1*In['midlog'][None,:]+b2*In['midlog'][None,:]**2+u0+u1*In['midlog'][None,:]
        logh=g0+g1*In['midlog'][None,:]+g2*In['gaplog'][None,:]+al*mmid
        h=np.exp(np.clip(logh,-30,30)); p=np.clip(1-np.exp(-h*In['delta'][None,:]),1e-12,1-1e-12)
        logsurv=_seg_sum(np.log1p(-p), In['starts'])        # (Q,n)
        pdmr=1-np.exp(logsurv)                              # (Q,n)
    else:
        sthr=np.exp(th[6]); lC=np.log1p(D['pat']['Cend'])[None,:]
        mC=b0+b1*lC+b2*lC**2+u0+u1*lC
        pdmr=ndtr((CDMR-mC)/sthr)                           # (Q,n)
    pred=(W*pdmr).sum(0)
    return pred, D['pat']['doc']

if __name__=='__main__':
    import time
    from sim_gen import simulate_cohort
    d=pack(simulate_cohort(120, informative=False, seed=3))
    for kind in ('floor','exact','threshold'):
        t=time.time(); res=fit(d,kind); dt=time.time()-t
        pr,ob=predict_dmr(res.x,d,kind); brier=np.mean((pr-ob)**2)
        print(f"{kind:9s} conv={res.success} nll={res.fun:.1f} t={dt:.1f}s "
              f"b1={res.x[1]:.2f}/{TRUTH['b1']} tau1={np.exp(res.x[5]):.2f}/{TRUTH['tau1']} "
              f"alpha={res.x[9] if kind!='threshold' else float('nan'):.2f} "
              f"brier={brier:.3f} pred={pr.mean():.2f} obs={ob.mean():.2f}")
