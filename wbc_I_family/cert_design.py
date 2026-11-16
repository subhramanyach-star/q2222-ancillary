"""Numerical design of a polynomial global certificate delta(c) = d0 + d1 (c-c0) + (c-c0)^2 sum_j R_j c^j
for a fixed functional g (SOCP: maximise eps s.t. Phi_c(delta(c)) <= 1 - eps (c-c0)^2 (1-c^2) on a c-grid,
and endpoint mixed terms M_{aa'eta}(delta(+-1)) <= 1 - marg). Coefficients rounded to rationals."""
import numpy as np, cvxpy as cp
from fractions import Fraction
from common import Phi
def Rmat(c):
    s = np.sqrt(max(1 - c*c, 0.0)); return np.array([[1.0, c], [0.0, s]])
def mixed_terms(g, dl, eps):
    """Endpoint lemma mixed terms M_{a a' eta} at c = eps (list, for numeric dl) -- expressions affine in dl."""
    f = np.array([1.0, eps]); b = g[2:4]; out = []
    for a in (1, -1):
        for a2 in (1, -1):
            for eta in (1, -1):
                u = b/2 + dl + a*g[4:6]; v = b/2 - dl + a2*g[6:8]
                out.append(a*g[0] + a2*g[1] + eta*((u - v) @ f))
    return out
def design(g, d0, d1, c0, k=3, nc=401, marg=0.02, weight=True):
    R = cp.Variable((k+1, 2)); eps = cp.Variable(); cons = []
    a = g[0:2]; b = g[2:4]; Cm = g[4:8].reshape(2, 2)
    cs = np.unique(np.concatenate([np.linspace(-1, 1, nc), c0 + np.linspace(-0.05, 0.05, 41)]))
    cs = cs[(cs >= -1) & (cs <= 1)]
    def dl_of(c): return d0 + d1*(c - c0) + (c - c0)**2 * (np.array([c**j for j in range(k+1)]) @ R)
    for c in cs:
        dl = dl_of(c); Rc = Rmat(c); tx = []
        for x in (0, 1):
            Bx = b/2 + (dl if x == 0 else -dl); t = cp.Variable(); tx.append(t)
            for sg in (1, -1): cons.append(sg*a[x] + cp.norm(Rc @ (Bx + sg*Cm[x])) <= t)
        w = (c - c0)**2 * (1 - c*c) if weight else (c - c0)**2
        cons.append(tx[0] + tx[1] <= 1 - eps*w)
    for e in (1, -1):
        for M in mixed_terms(g, dl_of(float(e)), e): cons.append(M <= 1 - marg)
    cons.append(cp.abs(R) <= 5)
    cp.Problem(cp.Maximize(eps), cons).solve(solver='CLARABEL')
    return eps.value, R.value
def rationalise(R, den=20000):
    return [[Fraction(float(x)).limit_denominator(den) for x in row] for row in R]
def check_float(g, d0, d1, c0, Rq, n=20001):
    Rf = np.array([[float(x) for x in row] for row in Rq]); k = len(Rf) - 1
    cs = np.concatenate([np.linspace(-1, 1, n), 1 - np.logspace(-10, -2, 30), -1 + np.logspace(-10, -2, 30)])
    worst = -1; wc = None
    for c in cs:
        dl = d0 + d1*(c - c0) + (c - c0)**2 * (np.array([c**j for j in range(k+1)]) @ Rf)
        v = Phi(g, c, dl) - 1
        if v > worst: worst, wc = v, c
    M = max(max(mixed_terms(g, d0 + d1*(e - c0) + (e - c0)**2*(np.array([e**j for j in range(k+1)]) @ Rf), e)) for e in (1., -1.))
    return worst, wc, M
