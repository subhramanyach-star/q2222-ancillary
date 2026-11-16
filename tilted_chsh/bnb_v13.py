"""Goal 2, vertices V1, V3: rigorous 2D interval branch-and-bound proof that
      F_{s0 s1}(c, theta) := R_{s0 s1} - sqrt(A_{s0}) - sqrt(B_{s1}) >= 0     on [-1,1] x Theta, all 4 sign pairs,
i.e. Phi_c(delta(c,theta)) <= 1 = g_v(theta).p_theta, for the certificate of certs_v13.json.
All arithmetic is arb ball arithmetic (python-flint); theta enters via exact closed forms (thm3_data.json).
A box [ca,cb] x [ta,tb] (theta = pi * x, x rational) is ACCEPTED by one of:
  (T1) mean-value form: F(cm, xm) - |F_c(box)| rc - |F_x(box)| rx > 0;
  (T2) Taylor at the forced touching curve c = c0(theta): F''_cc > 0 on hull([ca,cb] u c0([ta,tb])) x [ta,tb].
       Valid because F(c0(theta),theta) = 0 and F_c(c0(theta),theta) = 0 identically (Lemmas: complementary
       slackness + forced delta'(c0); checked as enclosures of 0 below);
  (T3) cb = 1 and F_c < 0 on the box, or ca = -1 and F_c > 0 on the box. Valid because F(+-1, theta) >= 0 for all
       theta in Theta (endpoint lemma, verified separately by endpoint_lemma()).
Every accepted box gives F > 0 except on the touching set {c = c0(theta)} u {c = +-1}: the equality set is contained
in {c0, -1, +1}."""
import pickle, json, sys, time
import sympy as sp, flint
from flint import arb, arb_series, fmpq
from fractions import Fraction
from multiprocessing import Pool
from arbeval import compile_expr, S_, C_, R_
flint.ctx.prec = 110
D = __import__('thm3_data').load()  # exact data from thm3_data.json
CERT = json.load(open('certs_v13.json'))
SYMS = [S_, C_, R_]
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
PI = arb.pi()
def q(x): return arb(fmpq(Fraction(x).numerator, Fraction(x).denominator))

class Vertex:
    def __init__(self, name):
        self.name = name; V = D['verts'][name]; ce = CERT[name]
        self.fg = [compile_expr(e, SYMS) for e in V['g']]
        self.fd0 = [compile_expr(e, SYMS) for e in V['d0']]
        self.fd1 = [compile_expr(e, SYMS) for e in V['d1']]
        self.k, self.L = ce['k'], ce['L']
        self.R = [[q(x) for x in row] for row in ce['R']]
        self.SM, self.SW = q(ce['SM']), q(ce['SW'])
    def data(self, x, series=False):
        """x: arb ball for theta/pi. series=True: Taylor series in x (order 1) with ball constant terms."""
        if series:
            th = arb_series([x * PI, PI], prec=2)
        else:
            th = x * PI
        s = (2 * th).sin(); C = (2 * th).cos(); r = (1 + s * s).sqrt()
        v = [s, C, r]
        g = [f(v) for f in self.fg]; d0 = [f(v) for f in self.fd0]; d1 = [f(v) for f in self.fd1]
        c0 = (1 - s * s) / (1 + s * s); sig = (s - self.SM) / self.SW
        return dict(g=g, d0=d0, d1=d1, c0=c0, sig=sig)
    def F(self, dat, c):
        g, d0, d1, c0, sig = dat['g'], dat['d0'], dat['d1'], dat['c0'], dat['sig']
        t = c - c0
        poly = [0, 0]
        for j in range(self.k + 1):
            for l in range(self.L + 1):
                w = (c ** j if j else 1) * (sig ** l if l else 1)
                row = self.R[j * (self.L + 1) + l]
                poly[0] = poly[0] + row[0] * w; poly[1] = poly[1] + row[1] * w
        dl = [d0[i] + d1[i] * t + t * t * poly[i] for i in range(2)]
        b = [g[2] / 2, g[3] / 2]
        out = []
        for s0, s1 in COMBOS:
            u = [b[i] + dl[i] + s0 * g[4 + i] for i in range(2)]
            w = [b[i] - dl[i] + s1 * g[6 + i] for i in range(2)]
            A = u[0] * u[0] + u[1] * u[1] + 2 * c * u[0] * u[1]
            B = w[0] * w[0] + w[1] * w[1] + 2 * c * w[0] * w[1]
            Rr = 1 - s0 * g[0] - s1 * g[1]
            out.append(Rr - A.sqrt() - B.sqrt())
        return out

def ball(a, b): return arb((a + b) / 2, (b - a) / 2)

def check_box(V, ca, cb, xa, xb, cache):
    """returns None if all 4 combos accepted, else list of failing combo indices."""
    key = (xa, xb)
    if key not in cache:
        xb_ = ball(xa, xb)
        cache[key] = (V.data(xb_), V.data(xb_, series=True), V.data(arb((xa + xb) / 2)))
    dB, dS, dM = cache[key]
    cm = (ca + cb) / 2; rc = (cb - ca) / 2; rx = (xb - xa) / 2
    cball = ball(ca, cb)
    Fm = V.F(dM, arb(cm))                                  # point c, narrow theta ball (midpoint)
    Fc = V.F(dB, arb_series([cball, 1], prec=3))           # series in c over box, theta ball
    Fx = V.F(dS, cball)                                    # series in x over box, c ball
    fail = []
    for i in range(4):
        fc = Fc[i].coeffs() if not isinstance(Fc[i], arb) else [Fc[i], arb(0), arb(0)]
        fx = Fx[i].coeffs()
        d1c = fc[1] if len(fc) > 1 else arb(0); d2c = 2 * fc[2] if len(fc) > 2 else arb(0)
        dx = fx[1] if len(fx) > 1 else arb(0)
        lb = Fm[i] - abs(d1c) * arb(rc) - abs(dx) * arb(rx)
        if lb > 0: continue                                                    # T1
        if cb == 1 and d1c < 0: continue                                       # T3
        if ca == -1 and d1c > 0: continue                                      # T3
        c0 = dB['c0']
        hb = cball.union(c0)                     # ball containing hull([ca,cb] u c0(theta-box))
        if hb.rad() < 0.3:
            F2 = V.F(dB, arb_series([hb, 1], prec=3))[i].coeffs()
            if len(F2) > 2 and F2[2] > 0: continue                             # T2
        fail.append(i)
    return fail or None


def run_strip(args):
    name, xa, xb, nc0 = args
    V = Vertex(name); cache = {}
    stack = [(fmpq(-1) + fmpq(2 * j, nc0), fmpq(-1) + fmpq(2 * (j + 1), nc0), xa, xb) for j in range(nc0)]
    nacc = 0; minw = 1
    while stack:
        ca, cb, a, b = stack.pop()
        f = check_box(V, ca, cb, a, b, cache)
        if f is None:
            nacc += 1; continue
        wc, wx = cb - ca, b - a
        if wc < fmpq(1, 2**22) and wx < fmpq(1, 2**26):
            return dict(ok=False, box=(float(ca), float(cb), float(a), float(b)), fail=f)
        # split the dimension with the larger estimated contribution
        if wc * 3 > wx * 64 or wx < fmpq(1, 2**26):
            m = (ca + cb) / 2; stack += [(ca, m, a, b), (m, cb, a, b)]
        else:
            m = (a + b) / 2; stack += [(ca, cb, a, m), (ca, cb, m, b)]
        if len(cache) > 400: cache.clear()
    return dict(ok=True, nboxes=nacc, strip=(float(xa), float(xb)))

def endpoint_lemma(name, xa, xb, n=256):
    """For all theta in [xa,xb]*pi: at c = +-1 the point D = delta0 +- delta1 lies where the x=0 term of
    Phi_{+-1} and the x=1 term have opposite slopes (+-1, strictly, on a neighbourhood), hence D is in a flat piece of Phi_{+-1}.
    On such a piece Phi_{+-1} = s0 a0 + s1 a1 + sig (b0 + e b1 + s0 (c00 + e c01) + s1 (c10 + e c11))  (e = +-1, sig = common
    sign, s0, s1 = active branches) = g.L for the deterministic L = (s0, s1, sig, sig e); and g.L <= 1 for all deterministic L
    (hexagon_m2.py).  Hence Phi_{+-1}(delta(+-1)) <= 1, i.e. F_{s0 s1}(+-1, theta) >= 0 for every sign pair."""
    V = Vertex(name); stack = [(xa + (xb - xa) * fmpq(j, n), xa + (xb - xa) * fmpq(j + 1, n)) for j in range(n)]
    cnt = 0
    while stack:
        a, b = stack.pop(); dat = V.data(ball(a, b)); ok = True
        for sgn in (1, -1):
            c = arb(sgn); g = dat['g']; t = c - dat['c0']
            poly = [0, 0]
            for j in range(V.k + 1):
                for l in range(V.L + 1):
                    w = (c ** j if j else 1) * (dat['sig'] ** l if l else 1); row = V.R[j * (V.L + 1) + l]
                    poly[0] += row[0] * w; poly[1] += row[1] * w
            dl = [dat['d0'][i] + dat['d1'][i] * t + t * t * poly[i] for i in range(2)]
            Dv = dl[0] + sgn * dl[1]; bp = (g[2] + sgn * g[3]) / 2
            # x=0 branches: s0 a0 + |bp + D + s0 gam0| ; need slope +1: for every possibly-active branch, bp + D + s0 gam0 > 0
            # x=1 branches: s1 a1 + |bp - D + s1 gam1| ; need slope -1: for every possibly-active branch, bp - D + s1 gam1 > 0
            # Phi_{+-1}(D) = psi0(D) + psi1(D), psi0 slope = sign(inner0 of active branch), psi1 slope = -sign(inner1 of active
            # branch).  Flat (=> minimal) iff all possibly-active inner values of both terms have one common strict sign.
            act = []
            for x, sl in ((0, 1), (1, -1)):
                gam = g[4 + 2 * x] + sgn * g[5 + 2 * x]
                inner = [bp + sl * Dv + s * gam for s in (1, -1)]
                val = [s * g[x] + abs(inner[k]) for k, s in enumerate((1, -1))]
                for k_ in range(2):
                    if val[k_] < val[1 - k_]: continue      # branch certainly inactive
                    act.append(inner[k_])
            if not (all(a > 0 for a in act) or all(a < 0 for a in act)): ok = False
        if not ok:
            if b - a < fmpq(1, 10**8): return False, (float(a), float(b))
            m = (a + b) / 2; stack += [(a, m), (m, b)]
        cnt += 1
    return True, cnt

def sanity(name, x):
    V = Vertex(name); dat = V.data(arb(x))
    c0 = dat['c0']; F = V.F(dat, arb_series([c0, 1], prec=3))
    return [(f.coeffs()[0], f.coeffs()[1]) for f in F]

if __name__ == "__main__":
    xa, xb = fmpq(1, 8), fmpq(3, 16)
    names = sys.argv[1:] or ['V1', 'V3']
    for name in names:
        print("== %s" % name, flush=True)
        print("   sanity at theta=pi/8 (F and F_c at c0 must enclose 0):", [(str(a), str(b)) for a, b in sanity(name, fmpq(1, 8))], flush=True)
        t0 = time.time(); ok, cnt = endpoint_lemma(name, xa, xb)
        print("   endpoint lemma on Theta: %s (%s theta-boxes, %.1fs)" % (ok, cnt, time.time() - t0), flush=True)
        NS = 64; t0 = time.time()
        strips = [(name, xa + (xb - xa) * fmpq(j, NS), xa + (xb - xa) * fmpq(j + 1, NS), 32) for j in range(NS)]
        with Pool(8) as P:
            res = P.map(run_strip, strips)
        allok = all(r['ok'] for r in res); nb = sum(r.get('nboxes', 0) for r in res)
        print("   B&B on [-1,1] x Theta: ALL ACCEPTED = %s ; %d accepted boxes ; %.1fs" % (allok, nb, time.time() - t0), flush=True)
        for r in res:
            if not r['ok']: print("   FAIL", r, flush=True)
        json.dump(dict(name=name, ok=allok, nboxes=nb, res=res), open('bnb_%s.json' % name, 'w'))
