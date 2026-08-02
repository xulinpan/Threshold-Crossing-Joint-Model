import numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr, logsumexp
cD,cF=-4.5,-5.0
TR=dict(b0=-2.0,b1=-3.7,b2=1.2,sig=0.8,t0=0.55,t1=0.75,sthr=0.15,srho=0.35)
srel_true=np.sqrt(TR['sthr']**2+TR['srho']**2)
def ell(t): return np.log1p(t)
def mq(a,g,b2,l): return a+g*l+b2*l*l
def simulate(n,seed):
    rng=np.random.default_rng(seed); pats=[]
    for i in range(n):
        u0=rng.normal(0,TR['t0']); u1=rng.normal(0,TR['t1']); a=TR['b0']+u0; g=TR['b1']+u1; b2=TR['b2']
        Ci=rng.uniform(3,10)
        C0=rng.normal(cD,TR['sthr']); Cr=rng.normal(cD,srel_true)
        t=[]; cur=0
        while True:
            cur+=rng.lognormal(np.log(0.55),0.5)
            if cur>Ci: break
            t.append(cur)
        if len(t)==0 or t[0]>0.3: t=[rng.uniform(.03,.12)]+list(t)
        t=np.sort(np.array(t)); y=mq(a,g,b2,ell(t))+rng.normal(0,TR['sig'],len(t)); fl=(y<=cF)
        yobs=np.where(fl,cF,y); edges=np.concatenate([[0],t])
        d0=g*g-4*b2*(a-C0); t_on=np.inf
        if d0>0:
            lon=(-g-np.sqrt(d0))/(2*b2); t_on=np.expm1(lon) if lon>=0 else (0 if mq(a,g,b2,0)<=C0 else np.inf)
        dr=g*g-4*b2*(a-Cr); t_off=np.inf
        if dr>0: t_off=np.expm1((-g+np.sqrt(dr))/(2*b2))
        onset=int(t_on<=Ci); relapse=int(onset and t_off<=Ci and t_off>t_on)
        def br(tt):
            k=np.searchsorted(t,tt); L=t[k-1] if k>0 else 0.0; R=t[k] if k<len(t) else Ci; return ell(L),ell(R)
        onL,onR=br(t_on) if onset else (ell(Ci),ell(Ci))
        rlL,rlR=br(t_off) if relapse else (ell(Ci),ell(Ci))
        pats.append(dict(ell=ell(t),y=yobs,cens=fl,bm=np.zeros(len(t)),onset=onset,onL=onL,onR=onR,
            ellC=ell(Ci),rel=relapse,rlL=rlL,rlR=rlR))
    return pats
pats=simulate(500,seed=5)
print("onset=%.2f relapse|onset=%.2f (srel_true=%.3f)"%(np.mean([p['onset'] for p in pats]),
    np.mean([p['rel'] for p in pats if p['onset']]),srel_true))
x,w=hermegauss(10); w=w/np.sqrt(2*np.pi)
GX0,GX1=np.meshgrid(x,x); GX0=GX0.ravel(); GX1=GX1.ravel(); LGW=np.log(np.outer(w,w).ravel())
def runmin(a,g,b2,lam): lv=-g/(2*b2); lc=np.clip(lv,0,lam); return a+g*lc+b2*lc*lc
def parts(th,fr):
    b0,b1=th[0],th[1]; b2=np.exp(th[2]); sig=np.exp(th[3]); t0=np.exp(th[4]); t1=np.exp(th[5]); sthr=np.exp(th[6])
    srel=np.sqrt(sthr**2+np.exp(th[7])**2) if fr else sthr
    return b0,b1,b2,sig,t0,t1,sthr,srel
def negpost(th,fr):
    b0,b1,b2,sig,t0,t1,sthr,srel=parts(th,fr); a=b0+t0*GX0; g=b1+t1*GX1; tot=0
    for pt in pats:
        muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2
        z=(pt['y'][:,None]-muo)/sig
        ll=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)
        if pt['onset']:
            MA=runmin(a,g,b2,pt['onL']); MB=runmin(a,g,b2,pt['onR'])
            ll=ll+np.log(np.clip(ndtr((MA-cD)/sthr)-ndtr((MB-cD)/sthr),1e-12,1))
            if pt['rel']:
                ll=ll+np.log(np.clip(ndtr((mq(a,g,b2,pt['rlR'])-cD)/srel)-ndtr((mq(a,g,b2,pt['rlL'])-cD)/srel),1e-12,1))
            else: ll=ll+log_ndtr((cD-mq(a,g,b2,pt['ellC']))/srel)
        else: ll=ll+log_ndtr((runmin(a,g,b2,pt['ellC'])-cD)/sthr)
        tot+=logsumexp(ll+LGW)
    lp=-0.5*((b0+2.5)/2)**2-0.5*((b1+1))**2-0.5*b2**2+th[2]-sig+th[3]-t0+th[4]-t1+th[5]-sthr+th[6]
    if fr: lp+=-np.exp(th[7])+th[7]
    return -(tot+lp)
def relcal(th,fr):
    b0,b1,b2,sig,t0,t1,sthr,srel=parts(th,fr); a=b0+t0*GX0; g=b1+t1*GX1; pr=[]; ob=[]
    for pt in pats:
        muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2; z=(pt['y'][:,None]-muo)/sig
        W=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)+LGW
        W=np.exp(W-W.max()); W/=W.sum()
        pon=1-ndtr((runmin(a,g,b2,pt['ellC'])-cD)/sthr)
        pr.append((W*ndtr((mq(a,g,b2,pt['ellC'])-cD)/srel)*pon).sum()); ob.append(pt['rel'] if pt['onset'] else 0)
    pr=np.array(pr); ob=np.array(ob); return pr.mean(),ob.mean(),np.mean((pr-ob)**2)
x0=np.array([TR['b0'],TR['b1'],np.log(TR['b2']),np.log(TR['sig']),np.log(TR['t0']),np.log(TR['t1']),np.log(0.2),np.log(0.3)])
for fr in [False,True]:
    r=minimize(negpost,x0,args=(fr,),method='L-BFGS-B',options=dict(maxiter=500,maxfun=5000))
    p=parts(r.x,fr); rp,ro,rb=relcal(r.x,fr)
    tag='frailty' if fr else 'base'
    print("\n== %s == nll=%.1f srel=%.3f sthr=%.3f srho=%s | relapse pred=%.3f obs=%.3f Brier=%.3f"%(
        tag,r.fun,p[7],p[6],("%.3f"%np.exp(r.x[7]) if fr else "--"),rp,ro,rb))
