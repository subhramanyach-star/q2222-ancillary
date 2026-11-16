#!/usr/bin/env python3
"""Standalone deterministic verifier for the certificates of Theorem 2 (vertices V1, V3 of the tilted-CHSH hexagon).

INPUT (text only, no pickles):  thm3_data.json  (exact vertex data g(theta), d0(theta), d1(theta), c0(theta))
                                certs_v13.json  (rational coefficients R_jl, sigma = (s - SM)/SW)
DEPENDENCIES: python3, sympy (parsing only), python-flint (Arb ball arithmetic).  Independent of the search code.

CLAIM VERIFIED.  For V in {V1, V3}, theta = pi*x, x in X = [1/8, 3/16], c in [-1, 1], and each sign pair (s0, s1):
    F_{s0 s1}(c, x) := T_{s0 s1} - sqrt(A_{s0}) - sqrt(B_{s1}) >= 0,
    A_{s0} = ||b/2 + delta + s0 e_0||^2_{G_c},  B_{s1} = ||b/2 - delta + s1 e_1||^2_{G_c},  T_{s0 s1} = 1 - s0 a0 - s1 a1,
    delta(c, theta) = d0 + d1 (c - c0) + (c - c0)^2 sum_{j<=k, l<=L} R_jl c^j sigma^l.
Since Phi_c(delta) = max_{s0,s1} (sqrt A + sqrt B + s0 a0 + s1 a1), this is Phi_c(delta(c,theta)) <= 1 = g.p_theta.

COVERING.  X is cut into 64 equal strips; in each strip [-1,1] is cut into 32 equal c-intervals.  Each box is either
ACCEPTED or bisected (in c if 3*wc > 64*wx, else in x).  Boxes are exact rationals.  The accepted boxes therefore tile
the strip; the verifier checks this by summing their exact areas (fmpq) and comparing with the strip area.

ACCEPTANCE TESTS (per box [ca,cb] x [xa,xb], all four sign pairs; every quantity is an Arb enclosure, and a test passes
only if the stated inequality holds for the whole enclosure; an enclosure that is not finite -- e.g. because a radicand
A or B is not bounded away from 0 on the box -- never passes, so every accepted box has A, B > 0 on it):
  (T1) mean value:  F(cm, xm) - sup|dF/dc|_box * rc - sup|dF/dx|_box * rx > 0          => F > 0 on the box.
  (T2) touching curve:  d^2F/dc^2 > 0 on H x [xa,xb], H = a ball containing [ca,cb] and c0([xa,xb]) (radius < 0.3).
       Valid because F(c0(theta), theta) = dF/dc(c0(theta), theta) = 0 identically on Theta (exact identities,
       exact_c0.py + exact_c0_signs.py).  Then F(c,theta) >= (c-c0)^2 * min F_cc / 2 >= 0, with equality iff c = c0.
  (T3) endpoint monotonicity:  cb = 1 and dF/dc < 0 on the box, or ca = -1 and dF/dc > 0 on the box.
       Valid because F(+-1, theta) >= 0 on Theta (endpoint lemma, verified below), so F > 0 for |c| < 1 on the box.
EQUALITY SET.  T1 gives strict positivity; T2 gives positivity except on c = c0(theta); T3 except at c = +-1.  Hence
{c : Phi_c(delta(c,theta)) = 1} is contained in {-1, c0(theta), 1}.
ENDPOINT LEMMA (paper: Lemma "Endpoint lemma", hypotheses H1-H2).  At c = +-1 the code checks, for THIS particular
certificate delta(+-1, theta) and all theta in Theta, that on each maximising branch the inner values have a common
sign (i.e. D = delta^(0) +- delta^(1) lies in the flat set of the piecewise-linear Phi_{+-1}); for this delta this
gives Phi_{+-1}(delta) = g.L for a deterministic L.  This is NOT a general property of Phi_{+-1}: without a tight
local point Phi_{+-1}(delta) can exceed h_{K_{+-1}}(g) (e.g. g = 0, delta = (1/4, 0)).  Combined with the check
g.L <= 1 for all 16 deterministic L (exact identity g.L = 1 by computer algebra, or g.L < 1 in Arb) it gives
F(+-1, theta) >= 0, which is what tests (T3) use.
[Header comment reworded 2026-10-21 for precision; code unchanged.]

USAGE:  python3 verify_thm3.py V1 [--strips j0 j1] [--prec 110] [--jobs N] [--dump boxes.txt]
        (default: all 64 strips; --strips selects strips j0..j1-1; --dump writes every accepted box with its test)
"""
import json, sys, time, argparse
import sympy as sp
import flint
from flint import arb, arb_series, fmpq

ap = argparse.ArgumentParser()
ap.add_argument('vertex', choices=['V1', 'V3'])
ap.add_argument('--strips', nargs=2, type=int, default=[0, 64])
ap.add_argument('--prec', type=int, default=110)
ap.add_argument('--dump', default=None)
ap.add_argument('--jobs', type=int, default=1, help='parallel worker processes (no effect on results)')
args = ap.parse_args()
flint.ctx.prec = args.prec
PI = arb.pi()

# ---------------- exact data -> Arb evaluators ----------------
s_, C_, r_ = sp.symbols('s C r', positive=True)
J = json.load(open('thm3_data.json')); CE = json.load(open('certs_v13.json'))[args.vertex]
V = J['vertices'][args.vertex]
def comp(expr_str):
    e = sp.sympify(expr_str, locals={'s': s_, 'C': C_, 'r': r_})
    def b(e):
        if e.is_Symbol:
            i = {s_: 0, C_: 1, r_: 2}[e]; return lambda v: v[i]
        if e.is_Rational: q = fmpq(int(e.p), int(e.q)); return lambda v: arb(q)
        if e.is_Add:
            fs = [b(a) for a in e.args]
            def f(v):
                acc = fs[0](v)
                for g in fs[1:]: acc = acc + g(v)
                return acc
            return f
        if e.is_Mul:
            fs = [b(a) for a in e.args]
            def f(v):
                acc = fs[0](v)
                for g in fs[1:]: acc = acc * g(v)
                return acc
            return f
        if e.is_Pow and e.exp.is_Integer:
            g = b(e.base); n = int(e.exp)
            return (lambda v: g(v) ** n) if n >= 0 else (lambda v: 1 / g(v) ** (-n))
        raise ValueError('unsupported: %s' % e)
    return b(e)
fg = [comp(x) for x in V['g']]; fd0 = [comp(x) for x in V['d0']]; fd1 = [comp(x) for x in V['d1']]
fc0 = comp(J['c0'])
q = lambda x: fmpq(*map(int, x.split('/'))) if '/' in x else fmpq(int(x))
K, L = CE['k'], CE['L']; RT = [[q(a) for a in row] for row in CE['R']]; SM, SW = q(CE['SM']), q(CE['SW'])
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

def theta_data(x):
    """x: arb or arb_series in x = theta/pi."""
    th = x * PI; s = (2 * th).sin(); C = (2 * th).cos(); r = (1 + s * s).sqrt(); v = [s, C, r]
    return dict(g=[f(v) for f in fg], d0=[f(v) for f in fd0], d1=[f(v) for f in fd1], c0=fc0(v), sig=(s - SM) / SW)

def F(dat, c):
    g, d0, d1, c0, sig = dat['g'], dat['d0'], dat['d1'], dat['c0'], dat['sig']
    t = c - c0; poly = [0, 0]
    for j in range(K + 1):
        for l in range(L + 1):
            w = (c ** j if j else 1) * (sig ** l if l else 1); row = RT[j * (L + 1) + l]
            poly[0] = poly[0] + row[0] * w; poly[1] = poly[1] + row[1] * w
    dl = [d0[i] + d1[i] * t + t * t * poly[i] for i in range(2)]
    out = []
    for s0, s1 in COMBOS:
        u = [g[2 + i] / 2 + dl[i] + s0 * g[4 + i] for i in range(2)]
        w = [g[2 + i] / 2 - dl[i] + s1 * g[6 + i] for i in range(2)]
        A = u[0] * u[0] + u[1] * u[1] + 2 * c * u[0] * u[1]
        B = w[0] * w[0] + w[1] * w[1] + 2 * c * w[0] * w[1]
        out.append(1 - s0 * g[0] - s1 * g[1] - A.sqrt() - B.sqrt())
    return out

def ball(a, b): return arb((a + b) / 2, (b - a) / 2)
def coeffs(z, n):
    cs = z.coeffs() if not isinstance(z, arb) else [z]
    return list(cs) + [arb(0)] * (n - len(cs))
def finite(z): return z.is_finite()

def accept(ca, cb, xa, xb, cache):
    key = (xa, xb)
    if key not in cache:
        cache[key] = (theta_data(ball(xa, xb)), theta_data(arb_series([ball(xa, xb), 1], prec=2)),
                      theta_data(arb((xa + xb) / 2)))
    dB, dS, dM = cache[key]
    cm, rc, rx = (ca + cb) / 2, (cb - ca) / 2, (xb - xa) / 2; cb_ = ball(ca, cb)
    Fm = F(dM, arb(cm)); Fc = F(dB, arb_series([cb_, 1], prec=3)); Fx = F(dS, cb_)
    tests = []
    for i in range(4):
        fc = coeffs(Fc[i], 3); fx = coeffs(Fx[i], 2); d1c, d2c, dx = fc[1], 2 * fc[2], fx[1]
        if all(finite(z) for z in (Fm[i], d1c, dx)):
            lb = Fm[i] - abs(d1c) * arb(rc) - abs(dx) * arb(rx)
            if lb > 0: tests.append('T1'); continue
        if finite(d1c) and cb == 1 and d1c < 0: tests.append('T3+'); continue
        if finite(d1c) and ca == -1 and d1c > 0: tests.append('T3-'); continue
        H = cb_.union(dB['c0'])
        if H.rad() < 0.3:
            f2 = coeffs(F(dB, arb_series([H, 1], prec=3))[i], 3)[2]
            if finite(f2) and f2 > 0: tests.append('T2'); continue
        return None
    return tests

def run_strip(j, dump):
    xa0 = fmpq(1, 8) + fmpq(1, 16) * fmpq(j, 64); xb0 = fmpq(1, 8) + fmpq(1, 16) * fmpq(j + 1, 64)
    stack = [(fmpq(-1) + fmpq(2 * i, 32), fmpq(-1) + fmpq(2 * (i + 1), 32), xa0, xb0) for i in range(32)]
    cache = {}; n = 0; area = fmpq(0); cnt = {}
    while stack:
        ca, cb, a, b = stack.pop()
        t = accept(ca, cb, a, b, cache)
        if t is not None:
            n += 1; area += (cb - ca) * (b - a)
            for x in t: cnt[x] = cnt.get(x, 0) + 1
            if dump: dump.write('%s %s %s %s %s\n' % (ca, cb, a, b, ','.join(t)))
            continue
        wc, wx = cb - ca, b - a
        if wc < fmpq(1, 2 ** 22) and wx < fmpq(1, 2 ** 26): return dict(ok=False, box=[str(ca), str(cb), str(a), str(b)])
        if wc * 3 > wx * 64 or wx < fmpq(1, 2 ** 26):
            m = (ca + cb) / 2; stack += [(ca, m, a, b), (m, cb, a, b)]
        else:
            m = (a + b) / 2; stack += [(ca, cb, a, m), (ca, cb, m, b)]
        if len(cache) > 400: cache.clear()
    full = fmpq(2) * (xb0 - xa0)
    return dict(ok=(area == full), nboxes=n, area_ok=(area == full), tests=cnt)

def run_strip_nodump(j): return run_strip(j, None)

# ---------------- endpoint lemma ----------------
DET = [(a0, a1, b0, b1) for a0 in (1, -1) for a1 in (1, -1) for b0 in (1, -1) for b1 in (1, -1)]
def zL(L_):
    a0, a1, b0, b1 = L_; return [a0, a1, b0, b1, a0 * b0, a0 * b1, a1 * b0, a1 * b1]
def gL_check():
    """g.L <= 1 on Theta for all 16 deterministic L: exact identity g.L == 1, or Arb-certified g.L < 1."""
    G = [sp.sympify(x, locals={'s': s_, 'C': C_, 'r': r_}) for x in V['g']]
    def red(e):
        n, d = sp.fraction(sp.together(sp.expand(e)))
        n = sp.expand(n)
        for _ in range(12):
            P = sp.Poly(n, C_, r_); n2 = 0
            for (i, k), cf in P.terms():
                n2 += cf * (1 - s_ ** 2) ** (i // 2) * C_ ** (i % 2) * (1 + s_ ** 2) ** (k // 2) * r_ ** (k % 2)
            n2 = sp.expand(n2)
            if n2 == n: break
            n = n2
        return n
    res = {}
    for L_ in DET:
        e = sum(gi * zi for gi, zi in zip(G, zL(L_))) - 1
        if red(e) == 0: res[L_] = 'exact =1'; continue
        f = comp(sp.sstr(e)); stack = [(fmpq(1, 8) + fmpq(1, 16) * fmpq(i, 64), fmpq(1, 8) + fmpq(1, 16) * fmpq(i + 1, 64)) for i in range(64)]
        while stack:
            a, b = stack.pop(); th = ball(a, b) * PI; sv = (2 * th).sin()
            val = f([sv, (2 * th).cos(), (1 + sv * sv).sqrt()])
            if val < 0: continue
            if b - a < fmpq(1, 10 ** 9): return False, L_
            m = (a + b) / 2; stack += [(a, m), (m, b)]
        res[L_] = '<1'
    return True, res

def endpoint_lemma(n=256):
    stack = [(fmpq(1, 8) + fmpq(1, 16) * fmpq(i, n), fmpq(1, 8) + fmpq(1, 16) * fmpq(i + 1, n)) for i in range(n)]; cnt = 0
    while stack:
        a, b = stack.pop(); dat = theta_data(ball(a, b)); ok = True; g = dat['g']
        for sgn in (1, -1):
            c = arb(sgn); t = c - dat['c0']; poly = [0, 0]
            for j in range(K + 1):
                for l in range(L + 1):
                    w = (c ** j if j else 1) * (dat['sig'] ** l if l else 1); row = RT[j * (L + 1) + l]
                    poly[0] += row[0] * w; poly[1] += row[1] * w
            dl = [dat['d0'][i] + dat['d1'][i] * t + t * t * poly[i] for i in range(2)]
            Dv = dl[0] + sgn * dl[1]; bp = (g[2] + sgn * g[3]) / 2; act = []
            for x, sl in ((0, 1), (1, -1)):
                gam = g[4 + 2 * x] + sgn * g[5 + 2 * x]
                inner = [bp + sl * Dv + s * gam for s in (1, -1)]
                val = [s * g[x] + abs(inner[k]) for k, s in enumerate((1, -1))]
                for k_ in range(2):
                    if val[k_] < val[1 - k_]: continue
                    act.append(inner[k_])
            if not (all(z > 0 for z in act) or all(z < 0 for z in act)): ok = False
        if not ok:
            if b - a < fmpq(1, 10 ** 8): return False, cnt
            m = (a + b) / 2; stack += [(a, m), (m, b)]
        cnt += 1
    return True, cnt

if __name__ == '__main__':
    t0 = time.time(); name = args.vertex
    print('verify_thm3  vertex=%s  prec=%d bits  python-flint %s  sympy %s' % (name, args.prec, flint.__version__, sp.__version__), flush=True)
    ok_gL, info = gL_check(); print('g.L <= 1 for all 16 deterministic L on Theta:', ok_gL, flush=True)
    ok_ep, n_ep = endpoint_lemma(); print('endpoint lemma on Theta:', ok_ep, '(%d theta-boxes)' % n_ep, flush=True)
    dump = open(args.dump, 'w') if args.dump else None
    allok = ok_gL and ok_ep; tot = 0; tests = {}
    js = list(range(*args.strips))
    if args.jobs > 1 and not dump:
        from multiprocessing import Pool
        with Pool(args.jobs) as P: results = P.map(run_strip_nodump, js)
    else:
        results = [run_strip(j, dump) for j in js]
    for j, r in zip(js, results):
        if not r['ok']: print('strip', j, 'FAILED', r, flush=True); allok = False; continue
        tot += r['nboxes']
        for k_, v in r['tests'].items(): tests[k_] = tests.get(k_, 0) + v
        print('strip %2d: %7d accepted boxes, exact area check %s, tests %s' % (j, r['nboxes'], r['area_ok'], r['tests']), flush=True)
    print('RESULT %s strips %d..%d: ALL OK = %s ; %d accepted boxes ; test counts %s ; %.1fs'
          % (name, args.strips[0], args.strips[1] - 1, allok, tot, tests, time.time() - t0), flush=True)
