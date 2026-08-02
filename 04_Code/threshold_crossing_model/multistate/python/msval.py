import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr, logsumexp
rng_global=np.random.default_rng(7)

# ---- truth: latent trajectory m(l)=a+g*l+b2*l^2, l=log(1+t) ----
TRUTH=dict(b0=-2.0,b1=-3.6,b2=1.30,tau0=0.55,tau1=0.75,sigma=0.8,sigmthr=0.15)
cD=-4.5; cF=-5.0; W=1.0   # durability window (years)
def ell(t): return np.log1p(t)
def mval(a,g,b2,l): return a+g*l+b2*l*l

def roots(a,g,b2,c):
    disc=g*g-4*b2*(a-c)
    if disc<=0: return None
    s=np.sqrt(disc); return ((-g-s)/(2*b2),(-g+s)/(2*b2))  # (l_down, l_up)

def simulate(n,seed=1,p=TRUTH):
    rng=np.random.default_rng(seed)
    long=[]; pat=[]
    for i in range(n):
        u0=rng.normal(0,p['tau0']); u1=rng.normal(0,p['tau1'])
        a=p['b0']+u0; g=p['b1']+u1; b2=p['b2']
        Ci=rng.uniform(3,10)
        Cthr=rng.normal(cD,p['sigmthr'])       # per-patient threshold
        # visits
        t=[]; cur=0.0
        while True:
            cur+=rng.lognormal(np.log(0.55),0.5)
            if cur>Ci: break
            t.append(cur)
        if len(t)==0 or t[0]>0.3: t=[rng.uniform(.03,.12)]+list(t)
        t=np.sort(np.array(t)); lv=ell(t)
        y=mval(a,g,b2,lv)+rng.normal(0,p['sigma'],len(t))
        fl=(y<=cF).astype(int); yobs=np.where(fl==1,cF,y)
        # true onset/relapse from roots vs Cthr
        rr=roots(a,g,b2,Cthr)
        lvtx=-g/(2*b2)
        if rr is None or rr[0]<0:   # never crosses down after 0 -> non-responder (or already below at 0)
            t_on=np.inf; t_off=np.inf
            if mval(a,g,b2,0.0)<=Cthr: t_on=0.0
        else:
            t_on=np.expm1(rr[0]); t_off=np.expm1(rr[1])
        onset=int(t_on<=Ci)
        relapse=int(onset and t_off<=Ci)
        sustained=int(onset and (t_off-t_on>=W))
        # documentation intervals (bracket by visits)
        def bracket(tt):
            if tt==0.0: return (0.0, t[0])
            k=np.searchsorted(t,tt)
            L=t[k-1] if k>0 else 0.0
            R=t[k] if k<len(t) else Ci
            return (L,R)
        onL,onR=bracket(t_on) if onset else (np.nan,np.nan)
        rlL,rlR=bracket(t_off) if relapse else (np.nan,np.nan)
        for tk,yk,fk in zip(t,yobs,fl): long.append((i,tk,ell(tk),yk,fk))
        pat.append((i,Ci,onset,onL,onR,relapse,rlL,rlR,sustained,t_on,t_off))
    L=np.array(long,float)
    P=np.array(pat,float)
    return dict(long=L,pat=P,n=n)

# ---- Gauss-Hermite grid ----
def gh(nq=9):
    x,w=hermegauss(nq); w=w/np.sqrt(2*np.pi)
    X0,X1=np.meshgrid(x,x); return X0.ravel(),X1.ravel(),np.log(np.outer(w,w).ravel())
GX0,GX1,LGW=gh(9)

def pack(d):
    n=d['n']; L=d['long']; P=d['pat']; pats=[]
    for i in range(n):
        lm=L[L[:,0]==i]; pr=P[P[:,0]==i][0]
        pats.append(dict(lt=lm[:,2],y=lm[:,3],cens=lm[:,4].astype(bool),
            Ci=pr[1],onset=int(pr[2]),onL=pr[3],onR=pr[4],
            relapse=int(pr[5]),rlL=pr[6],rlR=pr[7]))
    return pats

def runmin(a,g,b2,lam):   # running min of convex quad over [0,lam] (vectorized over nodes)
    lv=-g/(2*b2); lc=np.clip(lv,0,lam)
    return mval(a,g,b2,lc)

def neg_log_post(th,pats):
    b0,b1,b2=th[0],th[1],th[2]; sig=np.exp(th[3]); t0=np.exp(th[4]); t1=np.exp(th[5]); sthr=np.exp(th[6])
    a=b0+t0*GX0; g=b1+t1*GX1   # (Q,)
    ll=np.zeros(len(GX0))+LGW
    tot=0.0
    for pt in pats:
        lt=pt['lt'][:,None]  # (J,1) over nodes broadcast
        mu=b0+ (b1)*0  # placeholder
        # longitudinal
        muobs=(a[None,:]+g[None,:]*pt['lt'][:,None]+b2*pt['lt'][:,None]**2)  # (J,Q)
        z=(pt['y'][:,None]-muobs)/sig
        llo=-0.5*(z*z+np.log(2*np.pi))-np.log(sig)
        llc=log_ndtr((cF-muobs)/sig)
        lllong=np.where(pt['cens'][:,None],llc,llo).sum(0)  # (Q,)
        # onset term
        lamC=ell(pt['Ci'])
        if pt['onset']==1:
            MA=runmin(a,g,b2,ell(pt['onL'])); MB=runmin(a,g,b2,ell(pt['onR']))
            aa=(MA-cD)/sthr; bb=(MB-cD)/sthr
        
            llon=logsumexp(np.stack([log_ndtr(aa),log_ndtr(bb)]),axis=0,b=np.array([1.0,-1.0])[:,None])
        else:
            MC=runmin(a,g,b2,lamC); llon=log_ndtr((MC-cD)/sthr)
        # relapse term (only if onset)
        if pt['onset']==1 and pt['relapse']==1:
            mR2=mval(a,g,b2,ell(pt['rlR'])); mL2=mval(a,g,b2,ell(pt['rlL']))
            aa=(mR2-cD)/sthr; bb=(mL2-cD)/sthr
            llrl=logsumexp(np.stack([log_ndtr(aa),log_ndtr(bb)]),axis=0,b=np.array([1.0,-1.0])[:,None])
        elif pt['onset']==1 and pt['relapse']==0:
            mC=mval(a,g,b2,lamC); llrl=log_ndtr((cD-mC)/sthr)   # no relapse by C
        else:
            llrl=0.0
        tot_node=lllong+llon+(llrl if np.ndim(llrl) else llrl)+LGW
        tot+=logsumexp(tot_node)
    # weak priors
    lp=-0.5*((b0+2.5)/2)**2-0.5*((b1+1)/1)**2-0.5*((b2)/1.5)**2
    lp+=-sig+th[3]-t0+th[4]-t1+th[5]-sthr+th[6]
    return -(tot+lp)

if __name__=="__main__":
    import time
    d=simulate(450,seed=3)
    P=d['pat']
    print("n=%d  onset_rate=%.3f  relapse|onset=%.3f  sustained|onset=%.3f  floor=%.3f  visits/pt=%.1f"%(
        d['n'],P[:,2].mean(),P[P[:,2]==1][:,5].mean(),P[P[:,2]==1][:,8].mean(),
        d['long'][:,4].mean(),len(d['long'])/d['n']))
    pats=pack(d)
    x0=np.array([TRUTH['b0'],TRUTH['b1'],TRUTH['b2'],np.log(TRUTH['sigma']),
                 np.log(TRUTH['tau0']),np.log(TRUTH['tau1']),np.log(TRUTH['sigmthr'])])+rng_global.normal(0,0.1,7)
    t=time.time()
    res=minimize(neg_log_post,x0,args=(pats,),method='L-BFGS-B',options=dict(maxiter=300,maxfun=3000))
    print("fit %.1fs conv=%s"%(time.time()-t,res.success))
    names=['b0','b1','b2','sigma','tau0','tau1','sigmthr']
    est=res.x.copy(); est[3:]=np.exp(est[3:])
    tru=[TRUTH['b0'],TRUTH['b1'],TRUTH['b2'],TRUTH['sigma'],TRUTH['tau0'],TRUTH['tau1'],TRUTH['sigmthr']]
    print(f"{'param':8s}{'est':>9s}{'truth':>9s}")
    for nm,e,tv in zip(names,est,tru): print(f"{nm:8s}{e:9.3f}{tv:9.3f}")
