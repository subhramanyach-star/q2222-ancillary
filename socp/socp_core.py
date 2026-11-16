"""Exact quantum maximum h_Q(g) of any Bell functional in the CHSH scenario (marginals included) via
   h_Q(g) = max_{c in [-1,1]} min_{delta in R^2} Phi_c(delta),
   Phi_c(delta) = sum_x max_{a=+-1} [ a*a_x + || b/2 + (-1)^x delta + a e_x ||_{G_c} ],  G_c = [[1,c],[c,1]].
Coordinates: g = (a0,a1,b0,b1,e00,e01,e10,e11) acting on z = (<A0>,<A1>,<B0>,<B1>,<A0B0>,<A0B1>,<A1B0>,<A1B1>).
Inner problem: 2-variable SOCP (compiled once with cvxpy/Clarabel, then polished by direct minimisation of Phi;
any delta gives an upper bound on h_{K_c}). Outer problem: 1-D maximisation over c (grid + Brent)."""
import numpy as np, cvxpy as cp, time
from scipy.optimize import minimize, minimize_scalar

def Lc(c):
    return np.array([[1.0, 0.0], [c, np.sqrt(max(0.0, 1 - c * c))]])   # G_c = L L^T

def gnorm(u, c):
    return np.sqrt(max(0.0, u[0] ** 2 + u[1] ** 2 + 2 * c * u[0] * u[1]))

def Phi(g, c, d):
    g = np.asarray(g, float); a = g[0:2]; b = g[2:4]; e = [g[4:6], g[6:8]]
    tot = 0.0
    for x in (0, 1):
        s = 1 if x == 0 else -1
        tot += max(aa * a[x] + gnorm(b / 2 + s * d + aa * e[x], c) for aa in (1, -1))
    return tot

class SliceSolver:
    def __init__(self, g):
        self.g = np.asarray(g, float); g = self.g
        a = g[0:2]; b = g[2:4]; e = [g[4:6], g[6:8]]
        self.LT = cp.Parameter((2, 2))
        self.d = cp.Variable(2); t = cp.Variable(2)
        cons = []
        for x in (0, 1):
            s = 1 if x == 0 else -1
            for aa in (1, -1):
                u = b / 2 + s * self.d + aa * e[x]
                cons.append(t[x] >= aa * a[x] + cp.norm(self.LT @ u))
        self.prob = cp.Problem(cp.Minimize(cp.sum(t)), cons)
    def value(self, c, polish=True):
        self.LT.value = Lc(c).T
        try:
            self.prob.solve(solver='CLARABEL', tol_gap_abs=1e-12, tol_gap_rel=1e-12, tol_feas=1e-12)
            d0 = np.array(self.d.value, float)
        except Exception:
            d0 = np.zeros(2)
        v = Phi(self.g, c, d0)
        if polish:
            d1 = polish_slsqp(self.g, c, d0)
            v1 = Phi(self.g, c, d1)
            if v1 < v: v, d0 = v1, d1
        return v, d0

def polish_slsqp(g, c, d0):
    """Refine the 2-variable dual by SLSQP on its epigraph form; any delta is a valid upper bound."""
    a = g[0:2]; b = g[2:4]; e = [g[4:6], g[6:8]]; G = np.array([[1, c], [c, 1.]])
    def terms(d):
        out = []
        for x in (0, 1):
            s = 1 if x == 0 else -1
            for aa in (1, -1):
                u = b / 2 + s * d + aa * e[x]; out.append(aa * a[x] + np.sqrt(max(u @ G @ u, 1e-300)))
        return np.array(out)
    t0 = terms(d0); x0 = np.r_[d0, max(t0[0], t0[1]), max(t0[2], t0[3])]
    con = {'type': 'ineq', 'fun': lambda v: np.r_[v[2] - terms(v[:2])[:2], v[3] - terms(v[:2])[2:]]}
    try:
        r = minimize(lambda v: v[2] + v[3], x0, constraints=[con], method='SLSQP', options=dict(ftol=1e-16, maxiter=500))
        return r.x[:2]
    except Exception:
        return d0

def hQ_socp(g, ngrid=48):
    t0 = time.time()
    S = SliceSolver(g)
    cs = np.cos(np.pi * np.arange(ngrid + 1) / ngrid)[::-1]          # Chebyshev grid, dense near +-1
    vals = [S.value(c, polish=False)[0] for c in cs]
    k = int(np.argmax(vals))
    lo, hi = cs[max(k - 1, 0)], cs[min(k + 1, ngrid)]
    r = minimize_scalar(lambda c: -S.value(c)[0], bounds=(lo, hi), method='bounded',
                        options=dict(xatol=1e-13))
    best_c, best = r.x, -r.fun
    for c in (lo, hi, cs[k]):
        v = S.value(c)[0]
        if v > best: best, best_c = v, c
    v, d = S.value(best_c)
    return dict(value=max(best, v), c=best_c, delta=d, time=time.time() - t0)

# --- independent lower bound: explicit two-qubit strategy (A0=Z, A1 in xz-plane, B0=Z, B1 in xz-plane)
X = np.array([[0, 1], [1, 0.]]); Z = np.diag([1, -1.]); I2 = np.eye(2)
def W(g, a, b):
    A = [Z, np.cos(a) * Z + np.sin(a) * X]; B = [Z, np.cos(b) * Z + np.sin(b) * X]
    M = g[0] * np.kron(A[0], I2) + g[1] * np.kron(A[1], I2) + g[2] * np.kron(I2, B[0]) + g[3] * np.kron(I2, B[1])
    k = 4
    for x in (0, 1):
        for y in (0, 1):
            M = M + g[k] * np.kron(A[x], B[y]); k += 1
    return M
def hQ_qubit(g, starts=60, seed=1):
    g = np.asarray(g, float); rng = np.random.default_rng(seed)
    f = lambda p: -np.linalg.eigvalsh(W(g, p[0], p[1]))[-1]
    best = -np.inf
    for _ in range(starts):
        r = minimize(f, rng.uniform(0, 2 * np.pi, 2), method='Nelder-Mead', options=dict(xatol=1e-13, fatol=1e-15, maxiter=6000))
        best = max(best, -r.fun)
    return best

def hL(g):
    g = np.asarray(g, float); best = -np.inf
    for A0 in (1, -1):
        for A1 in (1, -1):
            for B0 in (1, -1):
                for B1 in (1, -1):
                    z = [A0, A1, B0, B1, A0 * B0, A0 * B1, A1 * B0, A1 * B1]
                    best = max(best, float(np.dot(g, z)))
    return best
