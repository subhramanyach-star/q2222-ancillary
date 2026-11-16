"""Main BFF runs: certified H(A1|X=1,E) for the tilted target p_theta under white noise, constraint = tilted-CHSH score
or full statistics (m=8 for all v, m=12 additionally for v>=0.97). Resumable: results appended to out/bff_main.jsonl."""
import numpy as np, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from bffcert import Rate
from families import ptheta, gtilt
OUT = 'out/bff_main.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        d = json.loads(line); done.add((round(d['theta'], 10), d['v'], d['expr'], d['m']))
THETAS = [('pi/8', np.pi/8), ('27.76deg', 27.76/180*np.pi), ('pi/16', np.pi/16), ('3pi/16', 3*np.pi/16)]
VS = [0.999, 0.995, 0.994, 0.99, 0.985, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 0.92, 0.90, 0.88, 0.86]
which = sys.argv[1:] if len(sys.argv) > 1 else [t[0] for t in THETAS]
R = {8: Rate(m=8, xstar=1), 12: Rate(m=12, xstar=1)}
for lab, th in THETAS:
    if lab not in which: continue
    p, c0, _ = ptheta(th); gt = gtilt(th)
    for v in VS:
        for m in (8, 12):
            if m == 12 and v < 0.97: continue
            for expr in ('tilted', 'full'):
                key = (round(th, 10), v, expr, m)
                if key in done: continue
                t0 = time.time()
                if expr == 'tilted': r = R[m].bound('score', g=gt, omega=v*(gt @ p))
                else: r = R[m].bound('full', z=v*p)
                rec = dict(label=lab, theta=th, v=v, expr=expr, m=m, H_cert=r.get('H_cert'), H_num=r.get('H_num'), status=r['status'],
                           K0=r.get('K0'), kappa=None if r.get('kappa') is None else list(map(float, r['kappa'])),
                           g=list(map(float, gt)) if expr == 'tilted' else None,
                           max_node_gap=None if 'nodes' not in r else max(n['gap'] for n in r['nodes']), secs=time.time() - t0)
                with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
                print('%s v=%.3f m=%d %-6s cert %.5f num %.5f (%.0fs)' % (lab, v, m, expr, rec['H_cert'] or np.nan, rec['H_num'] or np.nan, rec['secs']), flush=True)
