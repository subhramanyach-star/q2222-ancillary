"""EAT finite-size rates for protocols (a) CHSH score [CHSH device], (b) optimal witness g [lossycone device],
(c) full statistics [lossycone device], (c') full statistics [CHSH device] -- the Bancal-Sheridan-Scarani setting.
All use the same min-tradeoff construction (fs_core docstring) and the same EAT (fs_core.eat_rate = rates_v11/eat.py).
Free parameters optimised per (n, eps): (a),(b): mu >= 0 and tangent point s;  (c),(c'): f in R^8 and s.
eps_h = eps_EAT = eps (soundness/smoothing), completeness eps_c = 1e-6.  Output: out/eat_<eta>_<prot>.jsonl (resumable)."""
import sys, os, json, time; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, fs_core as F
from scipy.optimize import minimize
es, prot = sys.argv[1], sys.argv[2]; eta = float(es)
NS = [1e8, 1e9, 1e10, 1e11, 1e12]; EPS = [1e-6, 1e-10]
os.makedirs('out', exist_ok=True); OUT = 'out/eat_%s_%s.jsonl' % (es, prot)
done = {(r['n'], r['eps']): r for r in map(json.loads, open(OUT))} if os.path.exists(OUT) else {}
A = json.load(open('asym.json'))[es]
B = F.Beta('3')
if prot in ('a', 'cp'):
    from core import z_honest
    z = z_honest(A['chsh_x'], eta)
else:
    z = F.lossy_device(es)['p']
g = F.lossy_device(es)['g'] if prot == 'b' else F.CHSH
def kk(f, ls):
    be = B(f); Bh = be - f @ z; s = Bh*np.exp(ls)
    return F.tangent(be, f, s)
def rate(f, ls, n, eps):
    K0, ka = kk(f, ls); return F.eat_rate(K0, ka, z, n, eps, eps)
# asymptotic starting points
if prot in ('a', 'b'):
    mu0 = A['mu_a'] if prot == 'a' else A['mu_b']
    def obj(x, n, eps): return -rate(abs(x[0])*g, x[1], n, eps)
    x0s = [np.array([mu0, 0.0])]
else:
    fstar = np.array(A['fstar'])
    if prot == 'cp':
        from npa_float import pguess_full_dual
        _, _, fstar = pguess_full_dual('3', ('A', 1), z)
    sg = max((1, -1), key=lambda s: -np.log2(B(s*fstar) - s*fstar @ z)); fstar = sg*fstar
    def obj(x, n, eps): return -rate(x[:8], x[8], n, eps)
    x0s = [np.concatenate([fstar, [0.0]])]
    if prot == 'c':   # also start from the best (b) solution at each (n,eps) (so (c) >= (b) by construction)
        Bres = {(r['n'], r['eps']): r for r in map(json.loads, open('out/eat_%s_b.jsonl' % es))} if os.path.exists('out/eat_%s_b.jsonl' % es) else {}
prev = {}
for eps in EPS:
    for n in NS:
        if (n, eps) in done: prev[eps] = np.array(done[(n, eps)]['x']); continue
        t = time.time(); starts = list(x0s)
        if eps in prev: starts.append(prev[eps])
        if prot == 'c' and (n, eps) in Bres:
            rb = Bres[(n, eps)]; starts.append(np.concatenate([abs(rb['x'][0])*F.lossy_device(es)['g'], [rb['x'][1]]]))
        best = None
        for x0 in starts:
            if prot in ('a', 'b'):
                for ls0 in (x0[1], 0.0, 0.5, 1.5):
                    for m in (1.0, 0.5, 0.25):
                        y0 = np.array([x0[0]*m, ls0]); v = obj(y0, n, eps)
                        if best is None or v < best[0]: best = (v, y0)
            else:
                v = obj(x0, n, eps)
                if best is None or v < best[0]: best = (v, x0)
        x = best[1]
        for rnd in range(3 if prot in ('a', 'b') else 2):
            r = minimize(obj, x, args=(n, eps), method='Nelder-Mead',
                         options=dict(maxfev=300 if prot in ('a', 'b') else 1500, xatol=1e-7, fatol=1e-10, adaptive=True))
            if r.fun <= best[0]: best = (r.fun, r.x); x = r.x
        x = best[1]; f = abs(x[0])*g if prot in ('a', 'b') else x[:8]; ls = x[-1]
        K0, ka = kk(f, ls); R, info = F.eat_rate(K0, ka, z, n, eps, eps, full=True)
        rec = dict(eta=es, prot=prot, n=n, eps=eps, rate=R, x=list(map(float, x)), f=list(map(float, f)), beta_npa3=B(f),
                   K0=float(K0), kappa=list(map(float, ka)), z=list(map(float, z)), secs=time.time() - t, **info)
        open(OUT, 'a').write(json.dumps(rec) + '\n'); prev[eps] = x
        print(es, prot, '%.0e %.0e rate %.6f h %.6f gamma %.2e (%.0fs)' % (n, eps, R, info.get('h', 0), info.get('gamma', 0), rec['secs']), flush=True)
