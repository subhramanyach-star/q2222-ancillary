"""Closed forms (rational in a parameter) for the WBC families, in the d-coordinates
   g = g_corr + d_0 M_0 + d_1 M_1,  M_x = (G^{-1}E_x).(B_0,B_1) - A_x,
and the polygon vertices, forced certificate data d0(c0), d1(c0).
I-family: t = tan(delta/2):  sin d = 2t/(1+t^2), cos d = (1-t^2)/(1+t^2), E = (0, cos d, cos d, -sin 2d), c0 = -sin d.
J-family: mu = cos(gamma + pi/6): E = (3mu - 4mu^3, mu, mu, -mu), c0 = 1 - 2mu^2.
          randomness branch gamma in [0, pi/12] <-> mu in [1/sqrt2, sqrt3/2]; key branch gamma in [pi/12, pi/6] <-> mu in [1/2, 1/sqrt2]."""
import sympy as sp, itertools, pickle, sys
from exact_data import KEYS, DET, lab, LS
t, mu = sp.symbols('t mu', positive=True)
def fam(name):
    if name == 'I':
        S = 2*t/(1+t**2); C = (1-t**2)/(1+t**2); return t, [sp.Integer(0), C, C, -2*S*C], -S
    if name == 'J':
        return mu, [3*mu - 4*mu**3, mu, mu, -mu], 1 - 2*mu**2
def F(e): return sp.factor(sp.cancel(sp.together(e)))
def build(name):
    x, E, c0 = fam(name)
    p = sp.Matrix([0, 0, 0, 0] + E)
    G = sp.Matrix([[1, c0], [c0, 1]]); Gi = G.inv()
    gens = []; vs = []
    for a, xx in KEYS:
        w = a*sp.Matrix([E[2*xx], E[2*xx+1]]); tt = 1
        assert F((w.T*Gi*w)[0] - 1) == 0
        v = (Gi*w).applyfunc(F); vs.append(v)
        g = sp.zeros(8, 1); g[xx] = -a; g[2] = v[0]; g[3] = v[1]; g[4+2*xx] = a*v[0]; g[5+2*xx] = a*v[1]; gens.append(g)
    M = sp.Matrix.hstack(*gens).T
    zd = sp.zeros(8, 1)
    for xx in (0, 1): zd[5+2*xx] = (p[4+2*xx] - c0*p[5+2*xx])/(1 - c0**2)
    m = (M*zd).applyfunc(F)
    Lam0 = F(-m[1]/(m[0] - m[1])); Lam1 = F(1 - Lam0)
    gcorr = F(Lam0*(M[0, :] + M[2, :])/2 + Lam1*(M[1, :] + M[3, :])/2).T  # mu_{+|x}=mu_{-|x}=Lam_x/2
    gcorr = gcorr.applyfunc(F)
    # marginal directions: d_x = mu_{+|x} - mu_{-|x}; dg/dd_x = (gamma_{+|x} - gamma_{-|x})/2
    Mx = [((M[KEYS.index((1, xx)), :] - M[KEYS.index((-1, xx)), :])/2).T.applyfunc(F) for xx in (0, 1)]
    return dict(x=x, E=E, c0=c0, p=p, M=M, vs=vs, m=m, Lam=(Lam0, Lam1), gcorr=gcorr, Mx=Mx)
def local_lines(B):
    """g.L <= 1  <=>  d.ML <= 1 - gcorr.L  for each local L (lines in the d-plane)."""
    out = {}
    for L in DET:
        out[lab(L)] = (F((B['Mx'][0].T*L)[0]), F((B['Mx'][1].T*L)[0]), F(1 - (B['gcorr'].T*L)[0]))
    return out
def vertex_d(B, l1, l2):
    LL = local_lines(B); a1, b1, r1 = LL[l1]; a2, b2, r2 = LL[l2]
    det = a1*b2 - a2*b1
    return F((r1*b2 - r2*b1)/det), F((a1*r2 - a2*r1)/det)
def vertex_data(B, l1, l2):
    d = vertex_d(B, l1, l2)
    g = (B['gcorr'] + d[0]*B['Mx'][0] + d[1]*B['Mx'][1]).applyfunc(F)
    Lam = B['Lam']; muv = [F((Lam[0] + d[0])/2), F((Lam[1] + d[1])/2), F((Lam[0] - d[0])/2), F((Lam[1] - d[1])/2)]
    vs = B['vs']
    Bv0 = muv[0]*vs[0] + muv[2]*vs[2]
    d0 = (Bv0 - sp.Matrix([g[2], g[3]])/2).applyfunc(F)
    # forced derivative: sqrt(A_s(c0)) = 2 mu_{s|0}, sqrt(B_s(c0)) = 2 mu_{s|1}
    c, f0, f1 = sp.symbols('c f0 f1'); c0 = B['c0']
    dl = [d0[0] + f0*(c - c0), d0[1] + f1*(c - c0)]; rate = {}
    for s in (1, -1):
        u = [g[2+i]/2 + dl[i] + s*g[4+i] for i in range(2)]; A = u[0]**2 + u[1]**2 + 2*c*u[0]*u[1]
        w = [g[2+i]/2 - dl[i] + s*g[6+i] for i in range(2)]; Bq = w[0]**2 + w[1]**2 + 2*c*w[0]*w[1]
        rate[('A', s)] = sp.diff(A, c).subs(c, c0)/(4*muv[KEYS.index((s, 0))])
        rate[('B', s)] = sp.diff(Bq, c).subs(c, c0)/(4*muv[KEYS.index((s, 1))])
        assert F(A.subs(c, c0) - 4*muv[KEYS.index((s, 0))]**2) == 0 and F(Bq.subs(c, c0) - 4*muv[KEYS.index((s, 1))]**2) == 0
    eqs = [rate[('A', 1)] - rate[('A', -1)], rate[('A', 1)] + rate[('B', 1)], rate[('A', 1)] + rate[('B', -1)]]
    sol = sp.solve(eqs[:2], [f0, f1], dict=True)[0]
    res = F(eqs[2].subs(sol))
    d1 = sp.Matrix([F(sol[f0]), F(sol[f1])])
    return dict(d=d, g=g, mu=muv, d0=d0, d1=d1, consistency=res, tight=(l1, l2))
if __name__ == '__main__':
    out = {}
    for name in sys.argv[1:] or ['I', 'J']:
        B = build(name); print(name, 'Lambda0 =', B['Lam'][0], '\n gcorr =', list(B['gcorr']), '\n M0 =', list(B['Mx'][0]), '\n M1 =', list(B['Mx'][1]))
        out[name] = B
        out[name]['lines'] = local_lines(B)
    __import__('polygon_io').dump(out, 'out/sym_family.json')
