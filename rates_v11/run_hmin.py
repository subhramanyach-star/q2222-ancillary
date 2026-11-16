"""H_min(A1|X=1,E) = -log2 P_guess, NPA level-2 relaxation (real moment matrix, v7 npa.py), for the theta envelope:
constraints CHSH score, tilted-CHSH score, full statistics. Solver values only (NOT dual-certified); flagged as numerical."""
import numpy as np, json, warnings
warnings.filterwarnings('ignore')
from npa_hmin import pguess
from families import ptheta, gtilt
gch = np.array([0, 0, 0, 0, 1, 1, 1, -1.0]); out = []
THETAS = [('pi/16', np.pi/16), ('0.1pi', 0.1*np.pi), ('pi/8', np.pi/8), ('27.76deg', 27.76/180*np.pi), ('3pi/16', 3*np.pi/16), ('0.22pi', 0.22*np.pi), ('pi/4', np.pi/4)]
VS = [0.999, 0.994, 0.99, 0.98, 0.97, 0.95, 0.92, 0.90, 0.86]
for lab, th in THETAS:
    p, c0, _ = ptheta(th); gt = gtilt(th)
    for v in VS:
        rec = dict(label=lab, theta=th, v=v)
        for nm, kw in [('CHSH', dict(g=gch, omega=v*(gch @ p))), ('tilted', dict(g=gt, omega=v*(gt @ p))), ('full', dict(pobs=v*p))]:
            try:
                pg, st = pguess('2', ('A', 1), **kw); rec[nm] = float(-np.log2(min(max(pg, 0.5), 1))); rec[nm + '_status'] = st
            except Exception as e: rec[nm] = None; rec[nm + '_status'] = str(e)[:60]
        out.append(rec); print(rec, flush=True)
        json.dump(out, open('out/hmin_env.json', 'w'), indent=1)
