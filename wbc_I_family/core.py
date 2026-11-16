"""
Fixed-overlap NDWU slices K_c of Q_2222 (Ramanathan, arXiv:2610.08749, Eq. (15)).

Coordinates: z = (alpha0, alpha1, beta0, beta1, E00, E01, E10, E11)
   alpha_x = <A_x>, beta_y = <B_y>, E_xy = <A_x B_y>.
Functional g = (a0, a1, b0, b1, c00, c01, c10, c11), g.z = sum (no constant term).
Overlap c is BOB's (common-overlap condition on Bob's side, set A_B).

For each (a,x), a = +-1, x = 0,1:
   t_{a|x} = 1 + a alpha_x,  p_{a|x} = beta0 + a E_x0,  q_{a|x} = beta1 + a E_x1.
K_c (|c|<1):   || (q - c p, s p) ||_2 <= s t   for all (a,x),  s = sqrt(1-c^2)
K_{+-1}:       q = c p ,  |p| <= t            for all (a,x)   (LP; positivity explicit)

Dual (derived in dual_derivation.tex):
   h_{K_c}(g) = min_{delta in R^2} Phi_c(delta),
   Phi_c(delta) = sum_x max_{a=+-1} [ a a_x + || B_x(delta) + a c_x ||_{G_c} ],
   B_0 = b/2 + delta, B_1 = b/2 - delta, c_x = (c_x0, c_x1),
   ||v||_{G_c} = sqrt(v0^2 + v1^2 + 2 c v0 v1)   (valid for all c in [-1,1]).
"""
import numpy as np
import cvxpy as cp

IDX = dict(a0=0, a1=1, b0=2, b1=3, c00=4, c01=5, c10=6, c11=7)
PT = np.array([0, 0, 0, 0, 1, 1, 1, -1]) / np.sqrt(2) * np.array([0, 0, 0, 0, 1, 1, 1, 1]) + 0.0
PT = np.array([0, 0, 0, 0, 1/np.sqrt(2), 1/np.sqrt(2), 1/np.sqrt(2), -1/np.sqrt(2)])

def _tpq(z, a, x):
    t = 1 + a * z[x]
    p = z[2] + a * z[4 + 2 * x]
    q = z[3] + a * z[5 + 2 * x]
    return t, p, q

def primal_hK(g, c, solver="CLARABEL", return_z=False, form="paper"):
    """max g.z over K_c. form='paper' uses Eq.(15) literally; form='scaled' divides by s."""
    g = np.asarray(g, float)
    z = cp.Variable(8)
    cons = []
    if abs(c) >= 1:
        cs = np.sign(c)
        for a in (1, -1):
            for x in (0, 1):
                t, p, q = _tpq(z, a, x)
                cons += [q == cs * p, p <= t, -p <= t]
    else:
        s = np.sqrt(1 - c * c)
        for a in (1, -1):
            for x in (0, 1):
                t, p, q = _tpq(z, a, x)
                if form == "paper":
                    cons += [cp.SOC(s * t, cp.hstack([q - c * p, s * p]))]
                else:
                    cons += [cp.SOC(t, cp.hstack([(q - c * p) / s, p]))]
    prob = cp.Problem(cp.Maximize(g @ z), cons)
    prob.solve(solver=solver)
    if return_z:
        return prob.value, z.value
    return prob.value

def gram_norm(v, c):
    v = np.asarray(v, float)
    return np.sqrt(max(v[0] ** 2 + v[1] ** 2 + 2 * c * v[0] * v[1], 0.0))

def Phi(g, c, delta):
    """Dual objective (closed form after eliminating lambda). Upper bound on h_{K_c}(g) for ANY delta."""
    g = np.asarray(g, float)
    a = g[0:2]; b = g[2:4]; C = g[4:8].reshape(2, 2)
    delta = np.asarray(delta, float)
    tot = 0.0
    for x in (0, 1):
        Bx = b / 2 + (1 if x == 0 else -1) * delta
        tot += max(sgn * a[x] + gram_norm(Bx + sgn * C[x], c) for sgn in (1, -1))
    return tot

def dual_hK(g, c, solver="CLARABEL", return_y=False):
    """Conic dual: min sum lambda_{a|x} s.t. lambda >= ||eta||_{G_c}, a_x = -sum_a a lambda_{a|x},
       b_y = sum_{a,x} eta_y^{a|x}, c_xy = sum_a a eta_y^{a|x}.  Valid for all c in [-1,1]."""
    g = np.asarray(g, float)
    a = g[0:2]; b = g[2:4]; C = g[4:8].reshape(2, 2)
    s = np.sqrt(max(1 - c * c, 0.0))
    R = np.array([[1.0, c], [0.0, s]])  # R^T R = G_c
    lam = {(aa, x): cp.Variable() for aa in (1, -1) for x in (0, 1)}
    eta = {(aa, x): cp.Variable(2) for aa in (1, -1) for x in (0, 1)}
    cons = []
    for k in lam:
        cons += [cp.SOC(lam[k], R @ eta[k])]
    for x in (0, 1):
        cons += [a[x] == -(lam[(1, x)] - lam[(-1, x)])]
        cons += [C[x] == eta[(1, x)] - eta[(-1, x)]]
    cons += [b == sum(eta.values())]
    prob = cp.Problem(cp.Minimize(sum(lam.values())), cons)
    prob.solve(solver=solver)
    if return_y:
        e0p, e0m = eta[(1, 0)].value, eta[(-1, 0)].value
        delta = (e0p + e0m) - b / 2   # B_0 = eta^{+|0}+eta^{-|0}
        return prob.value, dict(delta=delta, lam={k: v.value for k, v in lam.items()},
                                eta={k: v.value for k, v in eta.items()})
    return prob.value

def dual_LP_endpoint(g, c, solver="HIGHS"):
    """Endpoint c=+-1: LP dual in probability form.
       min sum_{a,b,x} mu_{ab|x}, mu>=0, s.t.
         a_x      = -sum_{a,b} a mu_{ab|x}
         bt0      = -sum_{a,b,x} b mu_{ab|x}       (bt0 = b0 + c b1)
         ct_x0    = -sum_{a,b} a b mu_{ab|x}       (ct_x0 = c_x0 + c c_x1)
       multipliers of the 8 constraints  t_{a|x} + b p_{a|x} >= 0 (i.e. 4 P(a,b|x,0) >= 0)."""
    g = np.asarray(g, float)
    a = g[0:2]; b = g[2:4]; C = g[4:8].reshape(2, 2)
    bt0 = b[0] + c * b[1]; ct = C[:, 0] + c * C[:, 1]
    mu = {(aa, bb, x): cp.Variable(nonneg=True) for aa in (1, -1) for bb in (1, -1) for x in (0, 1)}
    cons = []
    for x in (0, 1):
        cons += [a[x] == -sum(aa * mu[(aa, bb, x)] for aa in (1, -1) for bb in (1, -1))]
        cons += [ct[x] == -sum(aa * bb * mu[(aa, bb, x)] for aa in (1, -1) for bb in (1, -1))]
    cons += [bt0 == -sum(bb * v for (aa, bb, x), v in mu.items())]
    prob = cp.Problem(cp.Minimize(sum(mu.values())), cons)
    prob.solve(solver=solver)
    return prob.value

def hK_fast(g, c, x0=None, restarts=4, seed=0):
    """h_{K_c}(g) via the 2-variable convex (nonsmooth) minimisation of Phi; returns (value, delta).
    Any delta gives a rigorous-in-exact-arithmetic UPPER bound Phi_c(delta) >= h_{K_c}(g)."""
    from scipy.optimize import minimize
    f = lambda d: Phi(g, c, d)
    rng = np.random.default_rng(seed)
    starts = [np.zeros(2)] + ([np.asarray(x0, float)] if x0 is not None else [])
    starts += [rng.normal(size=2) for _ in range(restarts)]
    best = None
    for s0 in starts:
        r = minimize(f, s0, method="Nelder-Mead", options=dict(xatol=1e-12, fatol=1e-14, maxiter=6000))
        for _ in range(3):   # restart from optimum (Nelder-Mead stalls on kinks)
            r2 = minimize(f, r.x + 1e-3 * rng.normal(size=2), method="Nelder-Mead",
                          options=dict(xatol=1e-13, fatol=1e-15, maxiter=6000))
            if r2.fun < r.fun: r = r2
        if best is None or r.fun < best.fun:
            best = r
    return best.fun, best.x

def hQ(g, ngrid=401, refine=True, method="socp"):
    """h_Q(g) = max_c h_{K_c}(g). Grid over c plus bounded local refinement."""
    from scipy.optimize import minimize_scalar
    cs = np.linspace(-1, 1, ngrid)
    f = (lambda c: dual_hK(g, c)) if method == "socp" else (lambda c: hK_fast(g, c)[0])
    vals = np.array([f(c) for c in cs])
    i = int(np.argmax(vals)); cbest, vbest = cs[i], vals[i]
    if refine:
        lo, hi = cs[max(i - 1, 0)], cs[min(i + 1, ngrid - 1)]
        if hi > lo:
            r = minimize_scalar(lambda c: -f(c), bounds=(lo, hi), method="bounded", options=dict(xatol=1e-10))
            if -r.fun > vbest:
                cbest, vbest = r.x, -r.fun
    return vbest, cbest, (cs, vals)


def hK_path(g, cgrid, restarts=1, seed=0):
    """h_{K_c}(g) along an increasing c-grid, warm-starting delta from the previous c (fast)."""
    from scipy.optimize import minimize
    rng = np.random.default_rng(seed)
    vals = np.empty(len(cgrid)); deltas = np.empty((len(cgrid), 2))
    prev = np.zeros(2)
    for i, c in enumerate(cgrid):
        f = lambda d: Phi(g, c, d)
        cands = [prev, np.zeros(2)]
        best = None
        for s0 in cands:
            r = minimize(f, s0, method="Nelder-Mead", options=dict(xatol=1e-11, fatol=1e-14, maxiter=3000))
            for _ in range(restarts):
                r2 = minimize(f, r.x + 1e-3 * rng.normal(size=2), method="Nelder-Mead",
                              options=dict(xatol=1e-12, fatol=1e-15, maxiter=3000))
                if r2.fun < r.fun: r = r2
            if best is None or r.fun < best.fun: best = r
        vals[i] = best.fun; deltas[i] = best.x; prev = best.x
    return vals, deltas
