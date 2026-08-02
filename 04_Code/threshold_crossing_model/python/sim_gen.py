"""
Data-generating model for the CML joint longitudinal--interval /
threshold-crossing simulation study (Section 3.8 of the manuscript).

Latent log-MRD trajectory (years, time on log(1+t) scale):
    m_i(t) = b0 + b1*log(1+t) + b2*{log(1+t)}^2 + u0_i + u1_i*log(1+t)
    u0_i ~ N(0, tau0^2), u1_i ~ N(0, tau1^2) independent.
Observation:  y = m_i(t) + eps, eps ~ N(0, sigma^2), left-censored at floor cF=-5.0.
True DMR onset:  first t with m_i(t) <= cDMR=-4.5  (latent threshold crossing).
Documented DMR: first visit with observed y <= -4.5 (measured endpoint).

Visits:
  non-informative: renewal process, gaps ~ LogNormal(median ~0.5y)
  informative:     non-homogeneous Poisson, intensity depends on latent MRD
                   lambda(t) = exp(d0 + d1*(m_i(t)+4.5) + d2*log(1+t) + w_i)
"""
import numpy as np

# ---- ground-truth parameters (posterior means of the fitted renewed model) ----
TRUTH = dict(
    b0=-2.074, b1=-3.561, b2=0.502,      # longitudinal fixed effects
    sigma=1.813, tau0=0.968, tau1=2.203, # residual sd, random-effect sds
    g0=-4.095, g1=-0.228, g2=-1.831, alpha=-1.275,  # interval-hazard params
)
CF = -5.0     # assay floor
CDMR = -4.5   # DMR threshold

def _m(logt, u0, u1, p):
    return p['b0'] + p['b1']*logt + p['b2']*logt**2 + u0 + u1*logt

def simulate_cohort(n, informative=False, dense=False, seed=0, p=TRUTH):
    rng = np.random.default_rng(seed)
    # administrative follow-up horizon per patient (years): heterogeneous
    C = rng.uniform(1.5, 7.0, size=n)
    # visit gap target (years); median ~0.5y (6 months); dense -> 0.33y (~4 months)
    gap_med = 0.33 if dense else 0.5
    gap_sig = 0.5  # lognormal shape
    d0, d1, d2, wsig = -0.7, 0.35, 0.2, 0.5  # informative-visit intensity params

    long_rows = []   # (pid, t, y, cens)  cens=1 if left-censored at floor
    interval_rows = []  # (pid, L, R, event)
    patient_rows = []   # (pid, true_onset, documented, C)

    for i in range(n):
        u0 = rng.normal(0, p['tau0'])
        u1 = rng.normal(0, p['tau1'])
        Ci = C[i]
        # ---- generate visit times in (0, Ci] ----
        if not informative:
            t, cur = [], 0.0
            while True:
                g = rng.lognormal(mean=np.log(gap_med), sigma=gap_sig)
                cur += g
                if cur > Ci:
                    break
                t.append(cur)
        else:
            # thinning of NHPP with intensity lambda(t)
            w = rng.normal(0, wsig)
            lam_max = None
            # coarse upper bound on intensity over [0,Ci]
            grid = np.linspace(1e-3, Ci, 50)
            lg = np.log1p(grid)
            mg = _m(lg, u0, u1, p)
            lam_grid = np.exp(d0 + d1*(mg + 4.5) + d2*lg + w)
            lam_max = lam_grid.max()*1.2 + 1e-6
            t, cur = [], 0.0
            while True:
                cur += rng.exponential(1.0/lam_max)
                if cur > Ci:
                    break
                lg = np.log1p(cur)
                mm = _m(lg, u0, u1, p)
                lam = np.exp(d0 + d1*(mm+4.5) + d2*lg + w)
                if rng.uniform() < lam/lam_max:
                    t.append(cur)
        # ensure a baseline visit
        if len(t) == 0 or t[0] > 0.25:
            t = [rng.uniform(0.02, 0.12)] + list(t)
        t = np.sort(np.array(t))
        logt = np.log1p(t)
        m = _m(logt, u0, u1, p)
        y = m + rng.normal(0, p['sigma'], size=len(t))
        cens = (y <= CF).astype(int)
        y_obs = np.where(cens == 1, CF, y)

        # ---- true latent onset (first crossing of -4.5) on a fine grid ----
        fg = np.linspace(1e-3, Ci, 400)
        mfg = _m(np.log1p(fg), u0, u1, p)
        below = np.where(mfg <= CDMR)[0]
        true_onset = fg[below[0]] if below.size else np.inf

        # ---- documented DMR: first visit with y_obs <= -4.5 ----
        dmr_visit = np.where(y_obs <= CDMR)[0]
        documented = int(dmr_visit.size > 0)

        # ---- at-risk intervals for interval model ----
        # intervals between consecutive visits (start at 0 -> first visit)
        edges = np.concatenate([[0.0], t])
        if documented:
            k_event = dmr_visit[0]  # index into t
            # intervals up to and including the event interval
            for k in range(k_event+1):
                L, R = edges[k], edges[k+1]
                ev = 1 if k == k_event else 0
                interval_rows.append((i, L, R, ev))
        else:
            for k in range(len(t)):
                L, R = edges[k], edges[k+1]
                interval_rows.append((i, L, R, 0))

        for tk, yk, ck in zip(t, y_obs, cens):
            long_rows.append((i, tk, yk, ck))
        patient_rows.append((i, true_onset, documented, Ci))

    import numpy as _np
    L = _np.array(long_rows, dtype=float)
    I = _np.array(interval_rows, dtype=float)
    P = _np.array(patient_rows, dtype=float)
    return dict(long=L, intr=I, pat=P, n=n)

if __name__ == "__main__":
    for inf in (False, True):
        d = simulate_cohort(120, informative=inf, seed=1)
        L, I, P = d['long'], d['intr'], d['pat']
        floor_rate = L[:,3].mean()
        dmr_rate = P[:,2].mean()
        vpp = len(L)/d['n']
        onset = P[:,1]
        obs_onset = onset[np.isfinite(onset)]
        print(f"informative={inf}: n={d['n']} obs/pt={vpp:.1f} floor_rate={floor_rate:.3f} "
              f"docDMR_rate={dmr_rate:.3f} med_true_onset(y)={np.median(obs_onset):.2f} "
              f"n_intervals={len(I)} event_intervals={int(I[:,3].sum())}")
