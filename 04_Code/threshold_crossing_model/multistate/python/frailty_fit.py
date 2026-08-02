import pandas as pd, numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import log_ndtr, ndtr, logsumexp
cD,cF=-4.5,-5.0
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
    pats.append(dict(ell=np.log1p(t),y=np.where(y<=cF,cF,y),cens=(y<=cF),bm=bm,
        onset=onset,onL=np.log1p(onL),onR=np.log1p(onR),ellC=np.log1p(Ci),
        rel=rel,rlL=np.log1p(rlL),rlR=np.log1p(rlR)))
print("onset=%d confirmed_relapse=%d"%(sum(p['onset'] for p in pats),sum(p['rel'] for p in pats)))
x,w=hermegauss(12); w=w/np.sqrt(2*np.pi)
GX0,GX1=np.meshgrid(x,x); GX0=GX0.ravel(); GX1=GX1.ravel(); LGW=np.log(np.outer(w,w).ravel())
def runmin(a,g,b2,lam): lv=-g/(2*b2); lc=np.clip(lv,0,lam); return a+g*lc+b2*lc*lc
def mq(a,g,b2,l): return a+g*l+b2*l*l
names=['beta0','beta_time','beta_time2','beta_bm','sigma_y','tau0','tau1','sigma_thr','sigma_rho']

def parts(th,frailty):
    b0,b1=th[0],th[1]; b2=np.exp(th[2]); bbm=th[3]; sig=np.exp(th[4]); t0=np.exp(th[5]); t1=np.exp(th[6]); sthr=np.exp(th[7])
    srel=np.sqrt(sthr**2+np.exp(th[8])**2) if frailty else sthr
    return b0,b1,b2,bbm,sig,t0,t1,sthr,srel

def negpost(th,frailty):
    b0,b1,b2,bbm,sig,t0,t1,sthr,srel=parts(th,frailty)
    a=b0+t0*GX0; g=b1+t1*GX1; tot=0.0
    for pt in pats:
        muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2+bbm*pt['bm'][:,None]
        z=(pt['y'][:,None]-muo)/sig
        ll=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)
        if pt['onset']:
            MA=runmin(a,g,b2,pt['onL']); MB=runmin(a,g,b2,pt['onR'])
            ll=ll+np.log(np.clip(ndtr((MA-cD)/sthr)-ndtr((MB-cD)/sthr),1e-12,1))
            if pt['rel']:
                mR=mq(a,g,b2,pt['rlR']); mL=mq(a,g,b2,pt['rlL'])
                ll=ll+np.log(np.clip(ndtr((mR-cD)/srel)-ndtr((mL-cD)/srel),1e-12,1))
            else:
                ll=ll+log_ndtr((cD-mq(a,g,b2,pt['ellC']))/srel)
        else:
            ll=ll+log_ndtr((runmin(a,g,b2,pt['ellC'])-cD)/sthr)
        tot+=logsumexp(ll+LGW)
    lp=-0.5*((b0+2.5)/2)**2-0.5*((b1+1)/1)**2-0.5*(bbm)**2-0.5*(b2)**2+th[2]
    lp+=-sig+th[4]-t0+th[5]-t1+th[6]-sthr+th[7]
    if frailty: lp+=-np.exp(th[8])+th[8]  # exp(1) on sigma_rho
    return -(tot+lp)

def calib(th,frailty):
    b0,b1,b2,bbm,sig,t0,t1,sthr,srel=parts(th,frailty)
    a=b0+t0*GX0; g=b1+t1*GX1; pon=[]; prel=[]; oon=[]; orel=[]
    for pt in pats:
        muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2+bbm*pt['bm'][:,None]
        z=(pt['y'][:,None]-muo)/sig
        W=np.where(pt['cens'][:,None],log_ndtr((cF-muo)/sig),-0.5*(z*z+np.log(2*np.pi))-np.log(sig)).sum(0)+LGW
        W=np.exp(W-W.max()); W/=W.sum()
        ponq=1-ndtr((runmin(a,g,b2,pt['ellC'])-cD)/sthr)
        prelq=ndtr((mq(a,g,b2,pt['ellC'])-cD)/srel)*ponq
        pon.append((W*ponq).sum()); prel.append((W*prelq).sum())
        oon.append(pt['onset']); orel.append(pt['rel'] if pt['onset'] else 0)
    pon=np.array(pon); prel=np.array(prel); oon=np.array(oon); orel=np.array(orel)
    return dict(onset_pred=pon.mean(),onset_obs=oon.mean(),onset_brier=np.mean((pon-oon)**2),
                rel_pred=prel.mean(),rel_obs=orel.mean(),rel_brier=np.mean((prel-orel)**2))

x0=np.array([-1.766,-3.623,np.log(0.456),0.730,np.log(1.842),np.log(1.052),np.log(2.055),np.log(0.178),np.log(0.40)])
import time,json
res={}
for frailty in [False,True]:
    t=time.time(); r=minimize(negpost,x0,args=(frailty,),method='L-BFGS-B',options=dict(maxiter=500,maxfun=5000))
    est=list(parts(r.x,frailty)); c=calib(r.x,frailty)
    tag='frailty' if frailty else 'base'
    print("\n== %s == nll=%.1f %.1fs conv=%s"%(tag,r.fun,time.time()-t,r.success))
    for nm,e in zip(names,est): print("  %-11s %8.3f"%(nm,e))
    print("  onset:   pred=%.3f obs=%.3f brier=%.3f"%(c['onset_pred'],c['onset_obs'],c['onset_brier']))
    print("  relapse: pred=%.3f obs=%.3f brier=%.3f"%(c['rel_pred'],c['rel_obs'],c['rel_brier']))
    res[tag]=dict(est={n:float(e) for n,e in zip(names,est)},calib={k:float(v) for k,v in c.items()},nll=float(r.fun),x=list(map(float,r.x)))
json.dump(res,open('frailty_results.json','w'),indent=2)
print("\nsaved frailty_results.json")
