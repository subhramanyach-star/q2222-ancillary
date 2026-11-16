"""NPA level-2 (real) guessing probability for the GLOBAL outcome pair (A_{x*},B_{y*}) in the 2x3 scenario (numerical
cross-check; H_min = -log2 P_g).  +-1 observables; z = (A0,A1,B0,B1,B2, A0B0,A0B1,A0B2, A1B0,A1B1,A1B2)."""
import itertools, numpy as np, cvxpy as cp
OA = ['A0', 'A1']; OB = ['B0', 'B1', 'B2']
def red(s):
    out = []
    for o in s:
        if out and out[-1] == o: out.pop()
        else: out.append(o)
    return tuple(out)
def canon(w):
    a = red([o for o in w if o[0] == 'A']); b = red([o for o in w if o[0] == 'B'])
    return min((a, b), (a[::-1], b[::-1]))
W = set()
for n in range(3):
    for w in itertools.product(OA + OB, repeat=n):
        W.add(canon(w)[0] + canon(w)[1]) if False else None
        a = red([o for o in w if o[0] == 'A']); b = red([o for o in w if o[0] == 'B']); W.add(a + b)
W = sorted(W, key=lambda w: (len(w), w)); n = len(W)
keys = {}; idx = np.zeros((n, n), int)
for i, u in enumerate(W):
    for j, v in enumerate(W):
        k = canon(tuple(reversed(u)) + v)
        if k not in keys: keys[k] = len(keys)
        idx[i, j] = keys[k]
m = len(keys)
F = [np.zeros((n, n)) for _ in range(m)]
for i in range(n):
    for j in range(n): F[idx[i, j]][i, j] = 1
k1 = keys[((), ())]
kz = [keys[(('A0',), ())], keys[(('A1',), ())]] + [keys[((), (b,))] for b in OB] + [keys[((a,), (b,))] for a in OA for b in OB]
def pguess(xs, ys, Eq, zeq, g=None, omega=None):
    ztot = 0; norm = 0; obj = 0; cons = []
    for sa in (1, -1):
        for sb in (1, -1):
            y = cp.Variable(m); G = sum(y[k] * F[k] for k in range(m)); cons.append(G >> 0)
            z = cp.hstack([y[k] for k in kz]); ztot = ztot + z; norm = norm + y[k1]
            obj = obj + (y[k1] + sa * z[xs] + sb * z[2 + ys] + sa * sb * z[5 + 3 * xs + ys]) / 4
    cons.append(norm == 1)
    if Eq: cons.append(ztot[Eq] == np.asarray(zeq))
    if g is not None: cons.append(np.asarray(g) @ ztot >= omega)
    prob = cp.Problem(cp.Maximize(obj), cons)
    for opts in [dict(solver='CLARABEL', tol_gap_abs=1e-9, tol_gap_rel=1e-9, tol_feas=1e-9, max_iter=500), dict(solver='SDPA')]:
        try:
            prob.solve(**opts)
            if prob.status in ('optimal', 'optimal_inaccurate'): return prob.value, prob.status + '/' + opts['solver']
        except Exception: pass
    return np.nan, 'failed'
