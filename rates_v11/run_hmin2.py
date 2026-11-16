"""NPA H_min cross-check for the global runs (same targets / protocols as run_global.py). -> out/global/hmin_global.jsonl"""
import numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
from npa_hmin2 import pguess
from targets2 import target, lift, best_key
from families import gtilt
from run_global_cfg import TH, VS, EQ
OUT = 'out/global/hmin_global.jsonl'; done = set()
if os.path.exists(OUT):
    for l in open(OUT): d = json.loads(l); done.add((d['label'], d['v'], d['proto']))
for v in VS:
    for lab, th in TH:
        h0, xs, phi = best_key(th); z = target(th, phi); g = lift(gtilt(th))
        for proto in ('witness', 'full'):
            if (lab, v, proto) in done: continue
            if proto == 'full': pg, st = pguess(xs, 2, list(range(11)), v * z)
            else: pg, st = pguess(xs, 2, EQ, v * z[EQ], g=g, omega=v * (g @ z))
            rec = dict(label=lab, v=v, proto=proto, xs=xs, Pg=pg, Hmin=float(-np.log2(pg)) if pg == pg else None, status=st)
            open(OUT, 'a').write(json.dumps(rec) + '\n'); print(rec, flush=True)
