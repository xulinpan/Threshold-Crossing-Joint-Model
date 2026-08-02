import pandas as pd, numpy as np
from numpy.polynomial.hermite_e import hermegauss
from scipy.optimize import minimize
from scipy.special import log_ndtr, logsumexp, ndtr
cD,cF=-4.5,-5.0

# rebuild per-patient with BOTH raw and confirmed relapse
d=pd.read_csv('rl.csv'); d['t']=d['t_months']/12.0; d=d.sort_values(['patient_num','t'])
pat=[]; longrows=[]
for pid,g in d.groupby('patient_num'):
    g=g.sort_values('t'); t=g['t'].values; y=g['log_mrd'].values; bm=g['sample_bm'].values
    edges=np.concatenate([[0.0],t]); Ci=t.max()
    dv=np.where(y<=cD)[0]; onset=int(len(dv)>0)
    onL=onR=Ci; rlL_r=rlR_r=Ci; rel_raw=0; rlL_c=rlR_c=Ci; rel_conf=0
    if onset:
        k=dv[0]; onL=edges[k]; onR=t[k]
        post=np.where(y[k+1:]>cD)[0]
        if len(post)>0:
            r=k+1+post[0]; rel_raw=1; rlL_r=t[r-1]; rlR_r=t[r]
        # confirmed: two consecutive >cD after onset
        for r in range(k+1,len(t)-0):
            if y[r]>cD and (r==len(t)-1 and False): pass
        for r in range(k+1,len(t)-1):
            if y[r]>cD and y[r+1]>cD:
                rel_conf=1; rlL_c=t[r-1]; rlR_c=t[r]; break
    pat.append((pid,Ci,onset,onL,onR,rel_raw,rlL_r,rlR_r,rel_conf,rlL_c,rlR_c))
    for tt,yy,bb in zip(t,y,bm):
        fl=int(yy<=cF); longrows.append((pid,np.log1p(tt),yy if fl==0 else cF,fl,int(bb)))
P=pd.DataFrame(pat,columns=['pid','Cend','onset','onL','onR','rel_raw','rlLr','rlRr','rel_conf','rlLc','rlRc'])
L=pd.DataFrame(longrows,columns=['pid','ell','y','floor','bm'])
print("onset=%d  relapse_raw=%d  relapse_confirmed=%d"%(P['onset'].sum(),P['rel_raw'].sum(),P['rel_conf'].sum()))

# ---- marginal-ML fit ----
x,w=hermegauss(10); w=w/np.sqrt(2*np.pi)
GX0,GX1=np.meshgrid(x,x); GX0=GX0.ravel(); GX1=GX1.ravel(); LGW=np.log(np.outer(w,w).ravel())
def ellf(t): return np.log1p(t)
def runmin(a,g,b2,lam):
    lv=-g/(2*b2); lc=np.clip(lv,0,lam); return a+g*lc+b2*lc*lc
def mq(a,g,b2,l): return a+g*l+b2*l*l
pats=[]
for pid,pr in P.iterrows():
    lm=L[L['pid']==pr['pid']]
    pats.append(dict(ell=lm['ell'].values,y=lm['y'].values,cens=lm['floor'].values.astype(bool),
        bm=lm['bm'].values,Ci=pr['Cend'],onset=int(pr['onset']),onL=ellf(pr['onL']),onR=ellf(pr['onR']),
        rawrel=int(pr['rel_raw']),rlLr=ellf(pr['rlLr']),rlRr=ellf(pr['rlRr']),
        confrel=int(pr['rel_conf']),rlLc=ellf(pr['rlLc']),rlRc=ellf(pr['rlRc']),ellC=ellf(pr['Cend'])))

def negpost(th,relmode='conf'):
    b0,b1=th[0],th[1]; b2=np.exp(th[2]); bbm=th[3]; sig=np.exp(th[4]); t0=np.exp(th[5]); t1=np.exp(th[6]); sthr=np.exp(th[7])
    a=b0+t0*GX0; g=b1+t1*GX1
    tot=0.0
    for pt in pats:
        muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2+bbm*pt['bm'][:,None]
        z=(pt['y'][:,None]-muo)/sig
        llo=-0.5*(z*z+np.log(2*np.pi))-np.log(sig)
        llc=log_ndtr((cF-muo)/sig)
        ll=np.where(pt['cens'][:,None],llc,llo).sum(0)
        # onset
        if pt['onset']==1:
            MA=runmin(a,g,b2,pt['onL']); MB=runmin(a,g,b2,pt['onR'])
            diff=np.clip(ndtr((MA-cD)/sthr)-ndtr((MB-cD)/sthr),1e-12,1.0); ll=ll+np.log(diff)
        else:
            MC=runmin(a,g,b2,pt['ellC']); ll=ll+log_ndtr((MC-cD)/sthr)
        # relapse
        if pt['onset']==1:
            rel=pt['confrel'] if relmode=='conf' else pt['rawrel']
            rlL=pt['rlLc'] if relmode=='conf' else pt['rlLr']; rlR=pt['rlRc'] if relmode=='conf' else pt['rlRr']
            if rel==1:
                mR=mq(a,g,b2,rlR); mL=mq(a,g,b2,rlL)
                diff=np.clip(ndtr((mR-cD)/sthr)-ndtr((mL-cD)/sthr),1e-12,1.0); ll=ll+np.log(diff)
            else:
                mC=mq(a,g,b2,pt['ellC']); ll=ll+log_ndtr((cD-mC)/sthr)
        tot+=logsumexp(ll+LGW)
    lp=-0.5*((b0+2.5)/2)**2-0.5*((b1+1)/1)**2-0.5*((b2)/1)**2-0.5*(bbm/1)**2
    lp+=-sig+th[4]-t0+th[5]-t1+th[6]-sthr+th[7]  # exp Jacobians; b2 half-normal via -0.5 b2^2 in lp above
    return -(tot+lp)

x0=np.array([-2.0,-3.5,np.log(0.5),0.3,np.log(1.7),np.log(1.0),np.log(2.0),np.log(0.4)])
names=['beta0','beta_time','beta_time2','beta_bm','sigma_y','tau0','tau1','sigma_thr']
def fit(mode):
    r=minimize(negpost,x0,args=(mode,),method='L-BFGS-B',options=dict(maxiter=400,maxfun=4000))
    est=r.x.copy(); est[2]=np.exp(est[2]); est[4:]=np.exp(est[4:])
    # Hessian SE (finite diff) on unconstrained, then delta approx: report on natural where simple
    return r,est
import time
import time,json
mode='conf'; t=time.time(); r,est=fit(mode)
print("== %s == nll=%.1f %.1fs conv=%s"%(mode,r.fun,time.time()-t,r.success))
# Wald SE via numeric Hessian on unconstrained params
k=len(r.x); H=np.zeros((k,k)); e=1e-3; f0=negpost(r.x,mode)
fp=np.array([negpost(r.x+np.eye(k)[i]*e,mode) for i in range(k)])
fm=np.array([negpost(r.x-np.eye(k)[i]*e,mode) for i in range(k)])
for i in range(k):
  for j in range(i,k):
    if i==j: H[i,i]=(fp[i]-2*f0+fm[i])/e**2
    else:
      fpp=negpost(r.x+(np.eye(k)[i]+np.eye(k)[j])*e,mode); fmm=negpost(r.x-(np.eye(k)[i]+np.eye(k)[j])*e,mode)
      H[i,j]=H[j,i]=(fpp-fp[i]-fp[j]+2*f0-(-fmm+fm[i]+fm[j]))/(2*e**2)
try: se_u=np.sqrt(np.clip(np.diag(np.linalg.inv(H)),0,None))
except: se_u=np.full(k,np.nan)
# transform SE to natural scale: params 2,4,5,6,7 are on log scale
se_nat=se_u.copy()
for idx in [2,4,5,6,7]: se_nat[idx]=est[idx]*se_u[idx]
out=[]
for nm,ev,sv in zip(names,est,se_nat):
    lo,hi=ev-1.96*sv,ev+1.96*sv
    print("  %-11s %8.3f  (%.3f, %.3f)"%(nm,ev,lo,hi)); out.append(dict(param=nm,est=float(ev),se=float(sv),lo=float(lo),hi=float(hi)))
json.dump(out,open('real_fit_conf.json','w'),indent=2)
np.save('real_fit_conf_x.npy',r.x)
print('saved real_fit_conf.json')


# ---------- sensitivity: raw relapse ----------
import json
r_raw=minimize(negpost,x0,args=('raw',),method='L-BFGS-B',options=dict(maxiter=400,maxfun=4000))
est_raw=r_raw.x.copy(); est_raw[2]=np.exp(est_raw[2]); est_raw[4:]=np.exp(est_raw[4:])
print("\n== raw == nll=%.1f conv=%s"%(r_raw.fun,r_raw.success))
for nm,e in zip(names,est_raw): print("  %-11s %8.3f"%(nm,e))

# ---------- transition calibration on real cohort (confirmed) using fitted params ----------
th=np.load('real_fit_conf_x.npy'); b0,b1=th[0],th[1]; b2=np.exp(th[2]); bbm=th[3]
sig=np.exp(th[4]); t0=np.exp(th[5]); t1=np.exp(th[6]); sthr=np.exp(th[7])
from scipy.special import ndtr
a=b0+t0*GX0; g=b1+t1*GX1
pon=[]; prel=[]; obs_on=[]; obs_rel=[]
for pt in pats:
    muo=a[None,:]+g[None,:]*pt['ell'][:,None]+b2*pt['ell'][:,None]**2+bbm*pt['bm'][:,None]
    z=(pt['y'][:,None]-muo)/sig
    llo=-0.5*(z*z+np.log(2*np.pi))-np.log(sig); llc=log_ndtr((cF-muo)/sig)
    W=np.where(pt['cens'][:,None],llc,llo).sum(0)+LGW; W=np.exp(W-W.max()); W/=W.sum()
    MC=runmin(a,g,b2,pt['ellC']); ponq=1-ndtr((MC-cD)/sthr)
    mC=mq(a,g,b2,pt['ellC']); prelq=ndtr((mC-cD)/sthr)*ponq
    pon.append(float((W*ponq).sum())); prel.append(float((W*prelq).sum()))
    obs_on.append(pt['onset']); obs_rel.append(pt['confrel'] if pt['onset']==1 else 0)
pon=np.array(pon); prel=np.array(prel); obs_on=np.array(obs_on); obs_rel=np.array(obs_rel)
print("\ncalibration (confirmed): onset pred=%.3f obs=%.3f | relapse pred=%.3f obs=%.3f"%(
    pon.mean(),obs_on.mean(),prel.mean(),obs_rel.mean()))
print("onset Brier=%.3f"%np.mean((pon-obs_on)**2))

# ---------- results table CSV ----------
res=json.load(open('real_fit_conf.json'))
import csv
with open('real_multistate_fit.csv','w',newline='') as fh:
    w=csv.writer(fh); w.writerow(['param','estimate','se','lo95','hi95','estimate_raw_relapse'])
    for row,er in zip(res,est_raw):
        w.writerow([row['param'],round(row['est'],4),round(row['se'],4),round(row['lo'],4),round(row['hi'],4),round(er,4)])
print("wrote real_multistate_fit.csv")

# ---------- Stan-ready data (JSON) for full HMC (confirmed relapse) ----------
# build Stan long arrays in patient order
pid=[]; ell=[]; y=[]; fl=[]; bm=[]
for i,pt in enumerate(pats,1):
    for e_,yy,cc,bb in zip(pt['ell'],pt['y'],pt['cens'],pt['bm']):
        pid.append(i); ell.append(float(e_)); y.append(float(yy)); fl.append(int(cc)); bm.append(int(bb))
def col(key): return [float(pt[key]) for pt in pats]
stan=dict(N=len(pats),Nobs=len(pid),pid=pid,ell=ell,y=y,floor_ind=fl,bm=bm,
    cF=cF,cD=cD,W=1.0,
    onset=[int(pt['onset']) for pt in pats],onL=col('onL'),onR=col('onR'),ellC=col('ellC'),
    relapse=[int(pt['confrel']) for pt in pats],rlL=col('rlLc'),rlR=col('rlRc'))
json.dump(stan,open('real_multistate_standata.json','w'))
print("wrote real_multistate_standata.json  (N=%d Nobs=%d)"%(stan['N'],stan['Nobs']))
