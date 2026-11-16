"""Certified BFF lower bounds on GLOBAL randomness H(AB|X=x*,Y=2,E), 2x3 scenario (Bob's extra key setting B2).
Targets: tilted-CHSH optimal strategy at theta + B2 chosen to maximise ideal H(AB) (targets2.best_key); theta=pi/4 is
maximally entangled CHSH (Tsirelson) with B2 = sx, x*=0 (ideal 2 bits).  White noise z -> v z.
Protocols: 'witness' = exposing witness score (tilted CHSH; CHSH at pi/4) >= v*g.p on the 2222 block + the three
B2-row correlators <B2>,<A0B2>,<A1B2> fixed; 'full' = all 11 correlators fixed.  m=8 Gauss-Radau.
Resumable: appends to out/global/bff_global.jsonl."""
import numpy as np, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from bffcert2 import Rate2
from targets2 import target, lift, best_key, pab, H
from families import gtilt
os.makedirs('out/global', exist_ok=True); OUT = 'out/global/bff_global.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        d = json.loads(line); done.add((d['label'], d['v'], d['proto']))
TH = [('pi/4', np.pi/4), ('pi/8', np.pi/8), ('27.76deg', 27.76/180*np.pi), ('3pi/16', 3*np.pi/16), ('0.22pi', 0.22*np.pi), ('0.1pi', 0.1*np.pi), ('pi/16', np.pi/16)]
VS = [0.99, 0.98, 0.999, 0.97, 0.95, 0.994]
R = {}
EQ = [4, 7, 10]
for v in VS:
    for lab, th in TH:
        h0, xs, phi = best_key(th); z = target(th, phi)
        if xs not in R: R[xs] = Rate2(m=8, nB=3, xs=xs, ys=2)
        g = lift(gtilt(th))
        for proto in ('witness', 'full'):
            if (lab, v, proto) in done: continue
            t0 = time.time()
            if proto == 'full': r = R[xs].bound(list(range(11)), v * z)
            else: r = R[xs].bound(EQ, v * z[EQ], g=g, omega=v * (g @ z))
            rec = dict(label=lab, theta=th, v=v, proto=proto, xs=xs, phi=phi, m=8, H_cert=r.get('H_cert'), H_num=r.get('H_num'),
                       status=r['status'], H_shannon_ideal=H(pab(v * z, xs, 2)), K0=r.get('K0'),
                       kappa=None if r.get('kappa') is None else list(map(float, r['kappa'])),
                       max_node_gap=None if 'nodes' not in r else max(n['gap'] for n in r['nodes']), secs=time.time() - t0)
            with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
            print('%-8s v=%.3f %-7s cert %.4f num %.4f shannon %.4f (%.0fs)' % (lab, v, proto, rec['H_cert'] or np.nan, rec['H_num'] or np.nan, rec['H_shannon_ideal'], rec['secs']), flush=True)
open('out/global/global.done', 'w').write('done\n')
