"""Envelope over theta: extra theta values (0.1pi, 0.22pi, pi/4) for certified BFF H(A1|X=1,E), tilted score and full stats.
Appends to out/bff_main.jsonl with the same schema as run_main.py (resumable)."""
import numpy as np, json, os, time, warnings
warnings.filterwarnings('ignore')
from bffcert import Rate
from families import ptheta, gtilt
OUT = 'out/bff_main.jsonl'; done = set()
for line in open(OUT):
    d = json.loads(line); done.add((d['label'], d['v'], d['expr'], d['m']))
THETAS = [('0.1pi', 0.1*np.pi), ('0.22pi', 0.22*np.pi), ('pi/4', np.pi/4)]
VS = [0.999, 0.994, 0.99, 0.98, 0.97, 0.95, 0.92, 0.90, 0.86]
R = {8: Rate(m=8, xstar=1), 12: Rate(m=12, xstar=1)}
for lab, th in THETAS:
    p, c0, _ = ptheta(th); gt = gtilt(th)
    for v in VS:
        for m in (8, 12):
            if m == 12 and v < 0.98: continue
            for expr in (('full',) if lab == 'pi/4' else ('tilted', 'full')):
                if (lab, v, expr, m) in done: continue
                t0 = time.time()
                r = R[m].bound('score', g=gt, omega=v*(gt @ p)) if expr == 'tilted' else R[m].bound('full', z=v*p)
                rec = dict(label=lab, theta=th, v=v, expr=expr, m=m, H_cert=r.get('H_cert'), H_num=r.get('H_num'), status=r['status'],
                           K0=r.get('K0'), kappa=None if r.get('kappa') is None else list(map(float, r['kappa'])),
                           g=list(map(float, gt)) if expr == 'tilted' else None,
                           max_node_gap=None if 'nodes' not in r else max(n['gap'] for n in r['nodes']), secs=time.time() - t0)
                with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
                print('%s v=%.3f m=%d %-6s cert %.5f (%.0fs)' % (lab, v, m, expr, rec['H_cert'] or np.nan, rec['secs']), flush=True)
