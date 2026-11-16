import sys, json; sys.path.insert(0, '.')
import numpy as np, fs_core as F
from scipy.optimize import minimize_scalar, minimize
from npa_float import pguess_full_dual
from core import Hmin_chsh
B = F.Beta('3'); out = {}
for es in F.ETAS:
    eta = float(es); d = F.lossy_device(es); c = F.chsh_device(eta)
    def Hb(mu, dv): return -np.log2(B(mu*dv['g']) - mu*dv['omega'])
    rb = minimize_scalar(lambda m: -Hb(m, d), bounds=(0, 30), method='bounded', options=dict(xatol=1e-9))
    rc = minimize_scalar(lambda m: -Hb(m, c), bounds=(0, 30), method='bounded', options=dict(xatol=1e-9))
    val, stt, fstar = pguess_full_dual('3', ('A', 1), d['p'])
    # Lemma 2 with f = -fstar? check sign: P <= beta(f) - f.p ; try both
    Hfull = max(-np.log2(B(s*fstar) - s*fstar @ d['p']) for s in (1, -1))
    valc, _, fc = pguess_full_dual('3', ('A', 1), c['p']); Hfullc = -np.log2(valc)
    out[es] = dict(mu_b=rb.x, H_b=-rb.fun, H_lo_t4=d['H_lo'], H_full_primal=float(-np.log2(val)), H_full_dual=float(Hfull),
                   S_chshdev=c['omega'], S_up=d['S_up'], mu_a=rc.x, H_a=-rc.fun, H_a_closed=Hmin_chsh(c['omega']), H_cap=d['H_cap'],
                   H_full_at_chshdev=float(Hfullc), chsh_x=c['x'], fstar=list(fstar))
    print(es, json.dumps({k: v for k, v in out[es].items() if k not in ('chsh_x', 'fstar')}), flush=True)
json.dump(out, open('asym.json', 'w'), indent=1)
