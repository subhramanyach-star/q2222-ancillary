"""Search for the best sigma-symmetric expression in the exact normal cone C_theta (hexagon regime) at each (theta, v):
g = l1 (V1+V1')/2 + l2 V2 + l3 (V3+V3')/2 + l4 V4, l in simplex (all such g satisfy g.p_theta = 1 and lie in C_theta).
Search with a reduced-trial certified bound (every evaluated number is still a rigorous lower bound), then re-certify the
winner with the full trial set at m=8 and m=12. Appends to out/bff_cone.jsonl (resumable)."""
import numpy as np, json, os, sys, time, itertools, warnings
warnings.filterwarnings('ignore')
from bffcert import Rate
from families import ptheta, gtilt, hexagon, hL
OUT = 'out/bff_cone.jsonl'
done = set()
if os.path.exists(OUT):
    for line in open(OUT):
        d = json.loads(line); done.add((d['label'], d['v']))
JOBS = [('pi/8', np.pi/8, [0.99, 0.98, 0.999, 0.97]), ('27.76deg', 27.76/180*np.pi, [0.99, 0.98, 0.999, 0.97, 0.95]),
        ('3pi/16', 3*np.pi/16, [0.99])]
fast = Rate(m=8, xstar=1); fast.TRIALS = fast.TRIALS[:2]; fast.HARD_NODES = 1; fast.ROWSCALE_ITERS = 1
full = {8: Rate(m=8, xstar=1), 12: Rate(m=12, xstar=1)}
for lab, th, VS in JOBS:
    p, c0, _ = ptheta(th); P, order, V = hexagon(th)
    G = [(order[0] + order[2])/2, order[1], (order[3] + order[5])/2, order[4]]
    gt = gtilt(th)
    # tilted CHSH in these coordinates (it is sigma-symmetric): least squares, for a starting point
    lt = np.linalg.lstsq(np.array(G).T, gt/(gt @ p), rcond=None)[0]
    for v in VS:
        if (lab, v) in done: continue
        t0 = time.time(); cache = {}
        def F(l):
            l = np.maximum(np.asarray(l, float), 0); l = l/l.sum(); key = tuple(np.round(l, 4))
            if key in cache: return cache[key]
            g = sum(a*b for a, b in zip(l, G))
            if hL(g) >= v - 1e-9: cache[key] = 0.0; return 0.0       # classical value reaches the noisy score: certifies nothing
            r = fast.bound('score', g=g, omega=v); val = r['H_cert'] if r['status'] == 'ok' else -1
            cache[key] = val; return val
        cands = [l for l in itertools.product(range(5), repeat=4) if sum(l) == 4]
        cands = [np.array(l, float)/4 for l in cands] + [np.maximum(lt, 0)/np.maximum(lt, 0).sum()]
        vals = [F(l) for l in cands]
        order_ = np.argsort(vals)[::-1]; best_l = cands[order_[0]]; best = vals[order_[0]]
        # local coordinate search on the simplex
        step = 0.125
        while step > 0.01:
            imp = False
            for i in range(4):
                for j in range(4):
                    if i == j: continue
                    l = best_l.copy(); d = min(step, l[j]); 
                    if d <= 0: continue
                    l[i] += d; l[j] -= d; val = F(l)
                    if val > best + 1e-6: best, best_l, imp = val, l, True
            if not imp: step /= 2
        g = sum(a*b for a, b in zip(best_l, G))
        rec = dict(label=lab, theta=th, v=v, lam=list(map(float, best_l)), g=list(map(float, g)), hL=float(hL(g)), H_search=best, n_evals=len(cache))
        for m in (8, 12):
            r = full[m].bound('score', g=g, omega=v)
            rec['H_cert_m%d' % m] = r.get('H_cert'); rec['K0_m%d' % m] = r.get('K0'); rec['kappa_m%d' % m] = None if r.get('kappa') is None else list(map(float, r['kappa']))
        rec['secs'] = time.time() - t0
        with open(OUT, 'a') as f: f.write(json.dumps(rec) + '\n')
        print(lab, v, 'best lam', np.round(best_l, 3), 'hL %.4f' % rec['hL'], 'H m8 %.5f m12 %.5f' % (rec['H_cert_m8'], rec['H_cert_m12']), 'evals', len(cache), '%.0fs' % rec['secs'], flush=True)
