import pandas as pd, numpy as np, itertools, time
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr, logsumexp
cD,cF=-4.5,-5.0; KN=np.log1p(2.0)   # knot at 2 years (ell scale)

# ---- build real-data per-patient ----
d=pd.read_csv('rl.csv'); d['t']=d['t_months']/12.0; d=d.sort_values(['patient_num','t'])
pats=[]
for pid,g in d.groupby('patient_num'):
    t=g['t'].values; y=g['log_mrd'].values; bm=g['sample_bm'].values
    edges=np.concatenate([[0.0],t]); Ci=t.max(); dv=np.where(y<=cD)[0]; onset=int(len(dv)>0)
    onL=onR=Ci; rel=0; rlL=rlR=Ci
    if onset:
        k=dv[0]; onL=edges[k]; onR=t[k]
        for r in range(k+1,len(t)-1):
            if y[r]>cD and y[r+1]>cD: rel=1; rlL=t[r-1]; rlR=t[r]; break
    pats.append(dict(ell=np.log1p(t),y=np.where(y<=cF,cF,y),cens=(y<=cF),bm=bm,onset=onset,
        onL=np.log1p(onL),onR=np.log1p(onR),ellC=np.log1p(Ci),rel=rel,rlL=np.log1p(rlL),rlR=np.log1p(rlR)))
print("onset=%d confirmed_relapse=%d"%(sum(p['onset'] for p in pats),sum(p['rel'] for p in pats)))

# ---- 3D Gauss-Hermite ----
nq=6; x,w=hermegauss(nq); w=w/np.sqrt(2*np.pi)
G=np.array(list(itertools.product(x,x,x)))            # (nq^3,3)
LW=np.log(np.array([wi*wj*wk for wi,wj,wk in itertools.product(w,w,w)]))
Z0,Z1,Z2=G[:,0],G[:,1],G[:,2]

def qmin(A,B,C,lo,hi):          # min of convex quad over [lo,hi], vectorized
    lv=-B/(2*A); lc=np.clip(lv,lo,hi); return A*lc*lc+B*lc+C
def mev(A,B,C,l): return A*l*l+B*l+C

def negpost(th,spline):
    b0,b1=th[0],th[1]; b2=np.exp(th[2]); bbm=th[3]; sig=np.exp(th[4])
    t0=np.exp(th[5]); t1=np.exp(th[6]); sz=np.exp(th[7]) if spline else 0.0; sthr=np.exp(th[8])
    a=b0+t0*Z0; g=b1+t1*Z1; ze=(sz*Z2) if spline else np.zeros_like(a)   # (Q,)
    A=b2; B1=g; C1=a; B2=g+ze; C2=a-ze*KN
    tot=0.0
    for pt in pats:
        # longitudinal: piecewise mean per obs
        below=pt['ell']<=KN
        muo=np.empty((len(pt['ell']),a.shape[0]))
        L=pt['ell'][:,None]
        muo=np.where(below[:,None], A*L*L+B1[None,:]*L+C1[None,:], A*L*L+B2[None,:]*L+C2[None,:])+bbm*pt['bm'][:,None]
        z=(pt['y'][:,None]-muo)/sig
        ll=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)
        # onset via piecewise running-min up to onL/onR
        def rmin(lam):
            m1=qmin(A,B1,C1,0.0,np.minimum(lam,KN))
            m2=np.where(lam>KN, qmin(A,B2,C2,KN,np.maximum(lam,KN)), np.inf)
            return np.minimum(m1,m2)
        if pt['onset']:
            MA=rmin(pt['onL']); MB=rmin(pt['onR'])
            ll=ll+np.log(np.clip(ndtr((MA-cD)/sthr)-ndtr((MB-cD)/sthr),1e-12,1))
            if pt['rel']:
                mR=mev(A,B2,C2,pt['rlR']); mL=mev(A,B2,C2,pt['rlL'])
                ll=ll+np.log(np.clip(ndtr((mR-cD)/sthr)-ndtr((mL-cD)/sthr),1e-12,1))
            else:
                mC=mev(A,B2,C2,pt['ellC']); ll=ll+log_ndtr((cD-mC)/sthr)
        else:
            ll=ll+log_ndtr((rmin(pt['ellC'])-cD)/sthr)
        tot+=logsumexp(ll+LW)
    lp=-0.5*((b0+2.5)/2)**2-0.5*(b1+1)**2-0.5*b2**2+th[2]-sig+th[4]-t0+th[5]-t1+th[6]-sthr+th[8]
    if spline: lp+=-sz+th[7]
    return -(tot+lp)

def relcal(th,spline):
    b0,b1=th[0],th[1]; b2=np.exp(th[2]); bbm=th[3]; sig=np.exp(th[4]); t0=np.exp(th[5]); t1=np.exp(th[6])
    sz=np.exp(th[7]) if spline else 0.0; sthr=np.exp(th[8])
    a=b0+t0*Z0; g=b1+t1*Z1; ze=(sz*Z2) if spline else np.zeros_like(a); A=b2; B1=g; C1=a; B2=g+ze; C2=a-ze*KN
    pr=[]; ob=[]
    for pt in pats:
        below=pt['ell']<=KN; L=pt['ell'][:,None]
        muo=np.where(below[:,None],A*L*L+B1[None,:]*L+C1[None,:],A*L*L+B2[None,:]*L+C2[None,:])+bbm*pt['bm'][:,None]
        z=(pt['y'][:,None]-muo)/sig
        W=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)+LW
        W=np.exp(W-W.max()); W/=W.sum()
        def rmin(lam):
            m1=qmin(A,B1,C1,0.0,np.minimum(lam,KN)); m2=np.where(lam>KN,qmin(A,B2,C2,KN,np.maximum(lam,KN)),np.inf); return np.minimum(m1,m2)
        pon=1-ndtr((rmin(pt['ellC'])-cD)/sthr)
        prel=ndtr((mev(A,B2,C2,pt['ellC'])-cD)/sthr)*pon
        pr.append((W*prel).sum()); ob.append(pt['rel'] if pt['onset'] else 0)
    pr=np.array(pr); ob=np.array(ob); return pr.mean(),ob.mean(),np.mean((pr-ob)**2)

x0=np.array([-1.77,-3.6,np.log(0.46),0.73,np.log(1.84),np.log(1.05),np.log(2.0),np.log(0.5),np.log(0.18)])
for spline in [False,True]:
    t=time.time(); r=minimize(negpost,x0,args=(spline,),method='L-BFGS-B',options=dict(maxiter=400,maxfun=4000))
    rp,ro,rb=relcal(r.x,spline); tag='spline' if spline else 'quadratic'
    print("%-9s nll=%.1f %.1fs conv=%s | sigma_zeta=%s | relapse pred=%.3f obs=%.3f Brier=%.3f"%(
        tag,r.fun,time.time()-t,r.success,("%.3f"%np.exp(r.x[7]) if spline else "0"),rp,ro,rb))
