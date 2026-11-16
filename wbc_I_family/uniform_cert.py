"""Uniform-in-parameter certificates for a polygon vertex of a one-parameter family (I: x=t=tan(delta/2), J: x=mu),
delta(c,x) = d0(x) + d1(x)(c - c0(x)) + (c - c0(x))^2 sum_{j<=k,l<=L} R_jl c^j sig^l, sig = (x - XM)/XW,
design by SOCP on an x-grid, rational rounding, then a rigorous 2D Arb branch-and-bound on [-1,1] x [xa,xb]
(same tests as prl anc thm3 verifier: T0 enclosure, T1 mean value, T2 convexity at the forced touching curve c=c0(x),
T3 monotone approach to c=+-1 with the endpoint lemma checked on x-boxes)."""
import sympy as sp, numpy as np, cvxpy as cp, flint, json, sys, time, pickle
from flint import arb, arb_series, fmpq
from fractions import Fraction
from multiprocessing import Pool
from arbeval import compile_expr
from cert_design import Rmat, mixed_terms
from common import Phi
flint.ctx.prec = 110
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

def ssqrt(A):
    """sqrt of a Gram-norm square (A >= 0 mathematically, |c| <= 1); if the enclosure meets 0, return an
    enclosure of [0, sqrt(upper(A))] (rigorous superset)."""
    if isinstance(A, arb) and not (A > 0):
        ub = A.upper()
        if not (ub > 0): return arb(0)
        return arb(0).union(arb(ub).sqrt())
    return A.sqrt()
def q(x): x = Fraction(x); return arb(fmpq(x.numerator, x.denominator))
class Family:
    def __init__(self, xsym, g, d0, d1, c0, XM, XW, k, L, R=None):
        self.xs = xsym; self.k, self.L = k, L; self.XM, self.XW = Fraction(XM), Fraction(XW)
        self.ex = dict(g=list(g), d0=list(d0), d1=list(d1), c0=c0)
        self.fg = [compile_expr(e, [xsym]) for e in g]; self.fd0 = [compile_expr(e, [xsym]) for e in d0]
        self.fd1 = [compile_expr(e, [xsym]) for e in d1]; self.fc0 = compile_expr(c0, [xsym])
        self.nf = sp.lambdify(xsym, sp.Matrix(list(g) + list(d0) + list(d1) + [c0]), 'numpy')
        self.R = R
    def num(self, x):
        v = np.array(self.nf(x), float).ravel(); return v[:8], v[8:10], v[10:12], v[12]
    def basis(self, c, x):
        sg = (x - float(self.XM))/float(self.XW); return np.array([c**j * sg**l for j in range(self.k+1) for l in range(self.L+1)])
    def design(self, xgrid, nc=241, marg=0.02, bound=5):
        R = cp.Variable(((self.k+1)*(self.L+1), 2)); eps = cp.Variable(); cons = []
        for x in xgrid:
            g, d0, d1, c0 = self.num(x); a = g[0:2]; b = g[2:4]; Cm = g[4:8].reshape(2, 2)
            cs = np.unique(np.concatenate([np.linspace(-1, 1, nc), c0 + np.linspace(-0.05, 0.05, 21)])); cs = cs[(cs >= -1) & (cs <= 1)]
            dl_of = lambda c: d0 + d1*(c - c0) + (c - c0)**2 * (self.basis(c, x) @ R)
            for c in cs:
                dl = dl_of(c); Rc = Rmat(c); tx = []
                for xx in (0, 1):
                    Bx = b/2 + (dl if xx == 0 else -dl); tv = cp.Variable(); tx.append(tv)
                    for sg in (1, -1): cons.append(sg*a[xx] + cp.norm(Rc @ (Bx + sg*Cm[xx])) <= tv)
                cons.append(tx[0] + tx[1] <= 1 - eps*(c - c0)**2*(1 - c*c))
            for e in (1., -1.):
                for M in mixed_terms(g, dl_of(e), e): cons.append(M <= 1 - marg)
        cons.append(cp.abs(R) <= bound)
        cp.Problem(cp.Maximize(eps), cons).solve(solver='CLARABEL')
        self.R = [[Fraction(float(v)).limit_denominator(20000) for v in row] for row in R.value]
        return eps.value
    def float_check(self, xs, n=4001):
        Rf = np.array([[float(v) for v in row] for row in self.R]); worst = -1; arg = None; Mmax = -1
        for x in xs:
            g, d0, d1, c0 = self.num(x)
            for c in np.concatenate([np.linspace(-1, 1, n), 1 - np.logspace(-9, -2, 15), -1 + np.logspace(-9, -2, 15)]):
                dl = d0 + d1*(c - c0) + (c - c0)**2*(self.basis(c, x) @ Rf); v = Phi(g, c, dl) - 1
                if v > worst: worst, arg = v, (x, c)
            for e in (1., -1.):
                Mmax = max(Mmax, max(mixed_terms(g, d0 + d1*(e - c0) + (e - c0)**2*(self.basis(e, x) @ Rf), e)))
        return worst, arg, Mmax
    # ---------- rigorous ----------
    def arbdata(self, xb, series=False):
        X = arb_series([xb, 1], prec=2) if series else xb
        return dict(g=[f([X]) for f in self.fg], d0=[f([X]) for f in self.fd0], d1=[f([X]) for f in self.fd1], c0=self.fc0([X]),
                    sig=(X - q(self.XM))/q(self.XW))
    def delta(self, dat, c):
        t = c - dat['c0']; poly = [0, 0]
        for j in range(self.k+1):
            for l in range(self.L+1):
                w = (c**j if j else 1)*(dat['sig']**l if l else 1); row = self.R[j*(self.L+1) + l]
                poly[0] = poly[0] + q(row[0])*w; poly[1] = poly[1] + q(row[1])*w
        return [dat['d0'][i] + dat['d1'][i]*t + t*t*poly[i] for i in range(2)]
    def F(self, dat, c):
        g = dat['g']; dl = self.delta(dat, c); b = [g[2]/2, g[3]/2]; out = []
        for s0, s1 in COMBOS:
            u = [b[i] + dl[i] + s0*g[4+i] for i in range(2)]; w = [b[i] - dl[i] + s1*g[6+i] for i in range(2)]
            A = u[0]*u[0] + u[1]*u[1] + 2*c*u[0]*u[1]; B = w[0]*w[0] + w[1]*w[1] + 2*c*w[0]*w[1]
            out.append(1 - s0*g[0] - s1*g[1] - ssqrt(A) - ssqrt(B))
        return out
    def endpoint_ok(self, xb):
        dat = self.arbdata(xb); g = dat['g']; b = [g[2]/2, g[3]/2]
        for e in (1, -1):
            c = arb(e); dl = self.delta(dat, c)
            for a in (1, -1):
                for a2 in (1, -1):
                    for eta in (1, -1):
                        u = [b[i] + dl[i] + a*g[4+i] for i in range(2)]; v = [b[i] - dl[i] + a2*g[6+i] for i in range(2)]
                        M = a*g[0] + a2*g[1] + eta*((u[0]-v[0]) + e*(u[1]-v[1]))
                        if not (M < 1): return False
        return True
def ball(a, b): return arb((a + b)/2, (b - a)/2)
def check_box(Fm, ca, cb, xa, xb, cache):
    key = (xa, xb)
    if key not in cache:
        X = ball(xa, xb); cache[key] = (Fm.arbdata(X), Fm.arbdata(X, series=True), Fm.arbdata(arb((xa + xb)/2)))
    dB, dS, dM = cache[key]; cB = ball(ca, cb); rc = arb((cb - ca)/2); rx = arb((xb - xa)/2)
    Fd = Fm.F(dB, cB)
    Fmid = Fm.F(dM, arb((ca + cb)/2))
    Fc = Fm.F(dB, arb_series([cB, 1], prec=3)); Fx = Fm.F(dS, cB)
    hb = cB.union(dB['c0']); Fh = Fm.F(dB, arb_series([hb, 1], prec=3)) if hb.rad() < 0.3 else None
    for i in range(4):
        if Fd[i] > 0: continue
        fc = Fc[i].coeffs(); fx = Fx[i].coeffs()
        dc = fc[1] if len(fc) > 1 else arb(0); dx = fx[1] if len(fx) > 1 else arb(0)
        if Fmid[i] - abs(dc)*rc - abs(dx)*rx > 0: continue
        if cb == 1 and dc < 0: continue
        if ca == -1 and dc > 0: continue
        if Fh is not None:
            ch = Fh[i].coeffs()
            if len(ch) > 2 and ch[2] > 0: continue
        return False
    return True
FAM = None
def run_strip(args):
    xa, xb, nc0 = args; Fm = FAM; cache = {}
    stack = [(fmpq(-1) + fmpq(2*j, nc0), fmpq(-1) + fmpq(2*(j+1), nc0), xa, xb) for j in range(nc0)]; nacc = 0
    while stack:
        ca, cb, a, b = stack.pop()
        if check_box(Fm, ca, cb, a, b, cache): nacc += 1; continue
        wc, wx = cb - ca, b - a
        if wc < fmpq(1, 2**24) and wx < fmpq(1, 2**30): return dict(ok=False, box=(float(ca), float(cb), float(a), float(b)))
        if (wc*4 > wx*64 and wc > fmpq(1, 2**24)) or wx < fmpq(1, 2**30):
            m = (ca + cb)/2; stack += [(ca, m, a, b), (m, cb, a, b)]
        else:
            m = (a + b)/2; stack += [(ca, cb, a, m), (ca, cb, m, b)]
        if len(cache) > 300: cache.clear()
    return dict(ok=True, nboxes=nacc)
def endpoint_strips(Fm, xa, xb, n=128):
    stack = [(xa + (xb - xa)*fmpq(j, n), xa + (xb - xa)*fmpq(j+1, n)) for j in range(n)]; cnt = 0
    while stack:
        a, b = stack.pop()
        if Fm.endpoint_ok(ball(a, b)): cnt += 1; continue
        if b - a < fmpq(1, 10**9): return False, (float(a), float(b))
        m = (a + b)/2; stack += [(a, m), (m, b)]
    return True, cnt
def prove(Fm, xa, xb, NS=64, nc0=32, procs=6):
    global FAM; FAM = Fm
    ok_e, ne = endpoint_strips(Fm, xa, xb)
    strips = [(xa + (xb - xa)*fmpq(j, NS), xa + (xb - xa)*fmpq(j+1, NS), nc0) for j in range(NS)]
    t0 = time.time()
    with Pool(procs) as P: res = P.map(run_strip, strips)
    return dict(endpoint_ok=ok_e, endpoint_boxes=ne, bnb_ok=all(r['ok'] for r in res), nboxes=sum(r.get('nboxes', 0) for r in res),
                fails=[r for r in res if not r['ok']], time_s=round(time.time() - t0, 1))
