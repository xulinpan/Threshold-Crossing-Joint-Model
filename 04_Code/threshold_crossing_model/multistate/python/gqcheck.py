import numpy as np
from scipy.stats import norm
from msval import simulate, TRUTH, cD, W, ell, mval

d=simulate(4000,seed=11,p=TRUTH); P=d['pat']; sthr=TRUTH['sigmthr']; b2=TRUTH['b2']
# recover per-patient true a,g from stored? not stored; recompute from pat via roots is circular.
# Instead regenerate a,g deterministically using same seed logic is hard; so compute predicted
# probabilities from the LATENT quantities we can reconstruct: use onset/relapse observed vs
# model formulas evaluated at each patient's own trajectory using stored t_on/t_off is trivial.
# Simplest valid check: the marginal predicted rates from the formula using the population
# distribution of (a,g) should match observed onset/sustained/relapse rates.
rng=np.random.default_rng(0); Nsim=20000
u0=rng.normal(0,TRUTH['tau0'],Nsim); u1=rng.normal(0,TRUTH['tau1'],Nsim)
a=TRUTH['b0']+u0; g=TRUTH['b1']+u1
C=rng.uniform(3,10,Nsim); lamC=ell(C)
def runmin(a,g,b2,lam):
    lv=-g/(2*b2); lc=np.clip(lv,0,lam); return a+g*lc+b2*lc**2
MC=runmin(a,g,b2,lamC)
p_onset=1-norm.cdf((MC-cD)/sthr)
mC=mval(a,g,b2,lamC)
p_relapse=norm.cdf((mC-cD)/sthr)*p_onset
disc=g**2-4*b2*(a-cD); p_sus=np.zeros(Nsim)
ok=disc>0
lon=np.where(ok,(-g-np.sqrt(np.where(ok,disc,0)))/(2*b2),0)
ton=np.expm1(np.clip(lon,0,None)); lW=np.log1p(ton+W)
p_sus=np.where(ok, norm.cdf((cD-mval(a,g,b2,lW))/sthr)*(1-norm.cdf((runmin(a,g,b2,lamC)-cD)/sthr)),0)
# empirical from a matched simulation using per-patient single Cthr draw
Cthr=rng.normal(cD,sthr,Nsim)
disc2=g**2-4*b2*(a-Cthr); ok2=disc2>0
lon2=np.where(ok2,(-g-np.sqrt(np.where(ok2,disc2,0)))/(2*b2),np.nan)
loff2=np.where(ok2,(-g+np.sqrt(np.where(ok2,disc2,0)))/(2*b2),np.nan)
ton2=np.expm1(np.clip(lon2,0,None)); toff2=np.expm1(loff2)
onset=(ok2)&(np.nan_to_num(ton2,nan=1e9)<=C)
relapse=onset&(np.nan_to_num(toff2,nan=1e9)<=C)
sustained=onset&((toff2-ton2)>=W)
print("            model_mean   empirical")
print("onset      %.3f       %.3f"%(p_onset.mean(),onset.mean()))
print("relapse    %.3f       %.3f"%(p_relapse.mean(),relapse.mean()))
print("sustained  %.3f       %.3f"%(p_sus.mean(),sustained.mean()))
