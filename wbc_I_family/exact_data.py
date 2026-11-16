"""Exact (sympy) data for zero-marginal maximally entangled xz-plane points:
p, c0, slice generators, overlap vector m, Lambda, g_corr, marginal directions M_x, polygon vertices,
forced certificate data d0 = delta(c0) and d1 = delta'(c0)."""
import sympy as sp, itertools
KEYS = [(1, 0), (1, 1), (-1, 0), (-1, 1)]
DET = [sp.Matrix([A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1]) for A0, A1, B0, B1 in itertools.product((1, -1), repeat=4)]
def lab(L): return ''.join('+' if s > 0 else '-' for s in L[:4])
LS = {lab(L): L for L in DET}
def simp(e): return sp.nsimplify(sp.radsimp(sp.simplify(e)))
def point(E, c0):
    p = sp.Matrix([0, 0, 0, 0] + list(E)); return p, c0
def slice_data(p, c0):
    G = sp.Matrix([[1, c0], [c0, 1]]); Gi = G.inv()
    gens, vs = [], []
    for a, x in KEYS:
        t = 1 + a*p[x]; w = sp.Matrix([p[2] + a*p[4+2*x], p[3] + a*p[5+2*x]])
        assert sp.simplify((w.T*Gi*w)[0] - t**2) == 0, 'cone not active'
        v = (Gi*w/t).applyfunc(sp.simplify)
        g = sp.zeros(8, 1); g[x] = -a; g[2] = v[0]; g[3] = v[1]; g[4+2*x] = a*v[0]; g[5+2*x] = a*v[1]
        gens.append(g); vs.append(v)
    M = sp.Matrix.hstack(*gens).T
    zd = sp.zeros(8, 1); zd[3] = (p[2] - c0*p[3])/(1 - c0**2)
    for x in (0, 1): zd[5+2*x] = (p[4+2*x] - c0*p[5+2*x])/(1 - c0**2)
    m = (M*zd).applyfunc(sp.simplify)
    return M, vs, m
def vertex(M, m, L1, L2):
    mu = sp.symbols('m0:4')
    eqs = [sum(mu) - 1, sum(mu[i]*m[i] for i in range(4))] + [sum(mu[i]*(M[i, :]*L)[0] for i in range(4)) - 1 for L in (L1, L2)]
    sol = sp.solve(eqs, mu, dict=True)[0]
    muv = sp.Matrix([sp.radsimp(sp.simplify(sol[k])) for k in mu])
    return muv, (M.T*muv).applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
def forced_d0(muv, vs, g):
    """delta(c0): eta^{a|x} = mu_{a|x} v_{a|x}; B_0 = eta^{+|0} + eta^{-|0}; delta = B_0 - b/2."""
    B0 = muv[0]*vs[0] + muv[2]*vs[2]
    return (B0 - sp.Matrix([g[2], g[3]])/2).applyfunc(lambda e: sp.radsimp(sp.simplify(e)))
c, e0, e1, f0, f1 = sp.symbols('c e0 e1 f0 f1')
def branch_AB(g, dl):
    out = {}
    for s0 in (1, -1):
        u = [g[2+i]/2 + dl[i] + s0*g[4+i] for i in range(2)]
        out[('A', s0)] = u[0]**2 + u[1]**2 + 2*c*u[0]*u[1]
    for s1 in (1, -1):
        w = [g[2+i]/2 - dl[i] + s1*g[6+i] for i in range(2)]
        out[('B', s1)] = w[0]**2 + w[1]**2 + 2*c*w[0]*w[1]
    return out
def forced_d1(g, d0, c0):
    """delta'(c0) from F_c(c0)=0 for all four sign pairs: A'_s/(2 sqrt A_s) = -B'_t/(2 sqrt B_t) for all s,t.
    Unknowns (f0,f1) = delta'(c0); equations a_+ = a_-, a_+ = -b_+, a_+ = -b_- (linear in f)."""
    dl = [e0 + f0*(c - c0), e1 + f1*(c - c0)]
    AB = branch_AB(g, dl)
    sub = {c: c0, e0: d0[0], e1: d0[1]}
    rate = {}
    for k, Q in AB.items():
        val = sp.radsimp(sp.simplify(Q.subs(sub)))
        der = sp.expand(sp.diff(Q, c).subs(sub))
        rate[k] = der/(2*sp.sqrt(val))
    eqs = [rate[('A', 1)] - rate[('A', -1)], rate[('A', 1)] + rate[('B', 1)], rate[('A', 1)] + rate[('B', -1)]]
    eqs = [sp.radsimp(sp.expand(e)) for e in eqs]
    sol = sp.solve(eqs[:2], [f0, f1], dict=True)[0]
    res = sp.radsimp(sp.simplify(eqs[2].subs(sol)))
    return sp.Matrix([sp.radsimp(sp.simplify(sol[f0])), sp.radsimp(sp.simplify(sol[f1]))]), res
