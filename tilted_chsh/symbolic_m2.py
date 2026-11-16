"""Milestone 2 (uniform in theta): exact closed forms in (s, C, r), s = sin 2theta, C = cos 2theta = sqrt(1-s^2),
r = sqrt(1+s^2) for: slice generators, Danskin vector m, the six hexagon vertices (mu, g), the forced c0-data
delta(c0), delta'(c0), and the forced endpoint value t* for the axis vertices.  Saved to m2_symbolic.pkl."""
import sympy as sp, pickle, time, itertools
s, C, r = sp.symbols('s C r', positive=True)
REL = {C**2: 1 - s**2, r**2: 1 + s**2}
def red(e):
    e = sp.together(sp.expand(e))
    n, d = sp.fraction(e)
    def rp(x):
        x = sp.expand(x)
        for _ in range(12):
            x2 = sp.Poly(x, C, r)
            terms = 0
            for (i, j), cf in x2.terms():
                terms += cf * (1 - s**2)**(i // 2) * C**(i % 2) * (1 + s**2)**(j // 2) * r**(j % 2)
            x2 = sp.expand(terms)
            if x2 == x: break
            x = x2
        return x
    n, d = rp(n), rp(d)
    # rationalise denominator: d = d0 + d1 C + d2 r + d3 C r
    for sym in (C, r):
        dc = sp.Poly(d, sym)
        if dc.degree() == 1:
            a1, a0 = dc.all_coeffs()
            conj = a0 - a1 * sym
            n, d = rp(n * conj), rp(d * conj)
    g = sp.gcd(sp.Poly(n, s, C, r), sp.Poly(d, s, C, r))
    n = sp.cancel(n / g.as_expr()); d = sp.cancel(d / g.as_expr())
    return sp.factor(n) / sp.factor(d)

KEYS = [(1, 0), (1, 1), (-1, 0), (-1, 1)]
PAIRS = {'V1': ((1, 1, 1, 1), (1, -1, -1, 1)), 'V2': ((1, 1, 1, -1), (1, -1, -1, 1)), 'V1s': ((1, 1, 1, -1), (1, -1, 1, 1)),
         'V3': ((1, -1, 1, 1), (-1, 1, 1, -1)), 'V4': ((-1, 1, 1, -1), (-1, -1, -1, 1)), 'V3s': ((1, 1, 1, 1), (-1, -1, -1, 1))}
def Lv(A0, A1, B0, B1): return sp.Matrix([A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1])
DET = list(itertools.product((1, -1), repeat=4))

def build():
    t0 = time.time()
    p = [C, 0, C / r, C / r, 1 / r, 1 / r, s**2 / r, -s**2 / r]
    c0 = (1 - s**2) / (1 + s**2)
    Gi = sp.Matrix([[1, -c0], [-c0, 1]]) / (1 - c0**2)
    gens, vs = {}, {}
    for a, x in KEYS:
        t = 1 + a * p[x]; w = sp.Matrix([p[2] + a*p[4 + 2*x], p[3] + a*p[5 + 2*x]])
        v = (Gi * w / t).applyfunc(red)
        assert red((w.T * Gi * w)[0] - t**2) == 0
        g = [0]*8; g[x] = -a; g[2], g[3] = v; g[4 + 2*x] = a*v[0]; g[5 + 2*x] = a*v[1]
        gens[(a, x)] = sp.Matrix(g); vs[(a, x)] = v
    zd = sp.zeros(8, 1)
    zd[3] = (p[2] - c0*p[3]) / (1 - c0**2); zd[5] = (p[4] - c0*p[5]) / (1 - c0**2); zd[7] = (p[6] - c0*p[7]) / (1 - c0**2)
    m = [red((gens[kk].T * zd)[0]) for kk in KEYS]
    print("v:", {k: list(v) for k, v in vs.items()}); print("m:", m, "(%.1fs)" % (time.time() - t0), flush=True)
    out = dict(p=p, c0=c0, gens={k: list(v) for k, v in gens.items()}, v={k: list(v) for k, v in vs.items()}, m=m, verts={})
    for name, (l1, l2) in PAIRS.items():
        mu = sp.symbols('mu0:4')
        g = sum((mu[i]*gens[KEYS[i]] for i in range(4)), sp.zeros(8, 1))
        eqs = [sum(mu) - 1, sum(mu[i]*m[i] for i in range(4)), (g.T*Lv(*l1))[0] - 1, (g.T*Lv(*l2))[0] - 1]
        sol = sp.solve(eqs, mu, dict=True)[0]
        muv = [red(sol[x]) for x in mu]
        gv = [red(sum(muv[i]*gens[KEYS[i]][j] for i in range(4))) for j in range(8)]
        b = sp.Matrix(gv[2:4])
        eta = {kk: muv[i]*vs[kk] for i, kk in enumerate(KEYS)}
        d0 = [red(e) for e in (eta[(1, 0)] + eta[(-1, 0)] - b/2)]
        # delta'(c0): for each x, sum_s s/nr_s [ sx * gram(u_s, d') + u_s0 u_s1 ] = 0, u_s = 2 eta^{s|x}, nr_s = 2 mu_{s|x}
        dp = sp.symbols('dp0 dp1'); eqs2 = []
        for x in (0, 1):
            sx = 1 if x == 0 else -1; e = 0
            for sg in (1, -1):
                u = 2*eta[(sg, x)]; nr = 2*muv[KEYS.index((sg, x))]
                gram = u[0]*dp[0] + u[1]*dp[1] + c0*(u[0]*dp[1] + u[1]*dp[0])
                e += sg / nr * (sx*gram + u[0]*u[1])
            eqs2.append(e)
        sol2 = sp.solve(eqs2, dp, dict=True)[0]
        d1 = [red(sol2[dp[0]]), red(sol2[dp[1]])]
        # non-tight constraints
        cons = {('det',) + L: red(1 - (sp.Matrix(gv).T*Lv(*L))[0]) for L in DET if L not in (l1, l2)}
        out['verts'][name] = dict(mu=muv, g=gv, d0=d0, d1=d1, cons=cons)
        print("== %s (%.1fs)\n   g = %s\n   mu = %s\n   d0 = %s\n   d1 = %s" % (name, time.time() - t0, gv, muv, d0, d1), flush=True)
    pickle.dump(out, open('m2_symbolic.pkl', 'wb'))
    return out

if __name__ == "__main__":
    build()
