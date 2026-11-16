#!/usr/bin/env python3
"""General per-strip verifier for global certificates of tilted-CHSH vertex functionals (fullrange project).

CLAIM.  For the vertex g = g(theta) named in the certificate file, theta = pi*x, every x in each strip [xa, xb], every
c in [-1, 1] and every sign pair (s0, s1):   F_{s0s1}(c, x) = T - sqrt(A_{s0}) - sqrt(B_{s1}) >= 0,  i.e.
Phi_c(delta(c, theta)) <= 1 = g.p_theta, with equality only for c in {-1, c0(theta), 1}.

CERTIFICATE (per strip):  delta(c) = d0 + d1 (c - c0) + (c - c0)^2 [ sum_j R_j T_j(c) + alpha + beta c ],
R_j in Q^2 (T_j Chebyshev), d0, d1, c0 exact functions of theta; alpha = beta = 0 unless the endpoint c = eps is PINNED
to the branch (x, a): then alpha, beta (exact functions of theta) make delta(eps) = delta_pin, where
delta_pin = -(b/2 + a e_0) (x = 0, so u_{a|0}(eps) = 0) or b/2 + a e_1 (x = 1, so v_{a|1}(eps) = 0).
For a pinned branch the vector is (c - eps) rho(c) with rho = (delta - delta_pin)/(c - eps) a polynomial (exact synthetic
division), and its G_c-norm is evaluated as (1 - eps c) ||rho||_{G_c}: this is analytic up to c = eps.

TESTS on a box [ca,cb] x [xa,xb] (Arb, outward rounding; non-finite enclosures never pass):
 T0 direct enclosure of F over the box > 0; T1 mean value; T2 F_cc > 0 on hull(box, c0(box)) (valid by the exact identities F = F_c = 0 at c0, exact_c0_all.py, given
 mu > 0, comb.py); T3 F_c < 0 on a box with cb = 1 (F_c > 0 with ca = -1); T3h the same on hull(box, +-1), valid since F(+-1, theta) >= 0 (endpoint check).
ENDPOINT CHECK (per x-box, eps = +-1): F_{s0s1}(eps) = min_{eta,eta'} (T - eta u.f - eta' v.f), f = (1, eps).  Terms with
eta' = eta, and all terms of a pair with a pinned (zero) vector, equal 1 - g.L for a deterministic L: these are >= 0 by
(H1), proved exactly in comb.py for the claimed theta-range.  Every other (mixed) term must be > 0 in Arb.
NOTE: H1 and mu > 0 come from comb.py for the CLAIMED range; strips may extend slightly beyond it (e.g. past theta*),
and only the claimed range is asserted.

USAGE: python3 vgen.py cert.json [--strips j0 j1] [--jobs N] [--prec 120]
"""
import json, sys, time, argparse
import sympy as sp, flint
from flint import arb, arb_series, fmpq
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('--strips', nargs=2, type=int, default=None)
ap.add_argument('--jobs', type=int, default=1); ap.add_argument('--prec', type=int, default=120)
ap.add_argument('--nc', type=int, default=32); ap.add_argument('--maxbox', type=int, default=3000000)
ap.add_argument('--verbose', action='store_true')
args = ap.parse_args(); VERB = args.verbose; flint.ctx.prec = args.prec; PI = arb.pi()
CERT = json.load(open(args.cert)); FILES = {'thm3': 'thm3_data.json', 'oct': 'oct_data.json'}
J = json.load(open(FILES[CERT['data']])); V = J['vertices'][CERT['vertex']]
s_, C_, r_ = sp.symbols('s C r', positive=True)
def comp(expr_str):
    e = sp.sympify(expr_str, locals={'s': s_, 'C': C_, 'r': r_})
    def b(e):
        if e.is_Symbol: i = {s_: 0, C_: 1, r_: 2}[e]; return lambda v: v[i]
        if e.is_Rational: q = fmpq(int(e.p), int(e.q)); return lambda v: arb(q)
        if e.is_Add or e.is_Mul:
            fs = [b(a) for a in e.args]; add = e.is_Add
            def f(v):
                acc = fs[0](v)
                for g in fs[1:]: acc = acc + g(v) if add else acc * g(v)
                return acc
            return f
        if e.is_Pow and e.exp.is_Integer:
            g = b(e.base); n = int(e.exp); return (lambda v: g(v) ** n) if n >= 0 else (lambda v: 1 / g(v) ** (-n))
        raise ValueError(e)
    return b(e)
FG = [comp(x) for x in V['g']]; FD0 = [comp(x) for x in V['d0']]; FD1 = [comp(x) for x in V['d1']]; FC0 = comp(J['c0'])
PINS = {int(k): v for k, v in CERT.get('pins', {}).items() if v and v[0] != 'm'}
MPINS = {int(k): v for k, v in CERT.get('pins', {}).items() if v and v[0] == 'm'}   # mid-pin: D(delta(eps)) := (K_i+K_j)/2
assert not (PINS and MPINS)
from midpin import cands
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
def q(x): x = str(x); return fmpq(*map(int, x.split('/'))) if '/' in x else fmpq(int(x))
def cheb_power(n):
    T = [[1], [0, 1]]
    while len(T) < n + 1:
        a, b = T[-1], T[-2]; t = [0] + [2 * v for v in a]
        for i, v in enumerate(b): t[i] -= v
        T.append(t)
    return T[:n + 1]
def theta_data(x):
    th = x * PI; s = (2 * th).sin(); C = (2 * th).cos(); r = (1 + s * s).sqrt(); v = [s, C, r]
    return dict(g=[f(v) for f in FG], d0=[f(v) for f in FD0], d1=[f(v) for f in FD1], c0=FC0(v))
def _mapd(f, d): return {k: ([f(z) for z in v] if isinstance(v, list) else f(v)) for k, v in d.items()}
def theta_data_mv(xa, xb):
    """(value enclosure over the x-box, series [value, derivative] over the x-box, midpoint data) using mean-value forms."""
    xm = (xa + xb) / 2; rx = (xb - xa) / 2
    dM = theta_data(arb(xm)); dSer = theta_data(arb_series([ball(xa, xb), 1], prec=2))
    der = _mapd(lambda z: coeffs(z, 2)[1], dSer)
    dB = {k: ([m + dd * arb(0, rx) for m, dd in zip(dM[k], der[k])] if isinstance(dM[k], list) else dM[k] + der[k] * arb(0, rx)) for k in dM}
    dS = {k: ([arb_series([b_, dd], prec=2) for b_, dd in zip(dB[k], der[k])] if isinstance(dB[k], list) else arb_series([dB[k], der[k]], prec=2)) for k in dB}
    return dB, dS, dM
def peval(coefs, c):
    acc = coefs[-1]
    for a in reversed(coefs[:-1]): acc = acc * c + a
    return acc
def build(dat, R):
    """theta-dependent pieces of the certificate: pin corrections alpha, beta and pinned vectors (exact formulas)."""
    g, d0, d1, c0 = dat['g'], dat['d0'], dat['d1'], dat['c0']
    pinv = {}
    for eps, (bx, ba) in PINS.items():
        if bx == 0: pinv[eps] = [-(g[2 + i] / 2 + ba * g[4 + i]) for i in range(2)]
        else: pinv[eps] = [g[2 + i] / 2 + ba * g[6 + i] for i in range(2)]
    al = [0, 0]; be = [0, 0]
    if MPINS:
        (e, (_, i, j)), = MPINS.items(); K = cands(g, e); tgt = (K[i] + K[j]) / 2
        Qe = [sum((R[k][m] * (e ** k) for k in range(len(R))), fmpq(0)) for m in range(2)]
        Dal = (tgt - (d0[0] + d1[0] * (e - c0)) - e * (d0[1] + d1[1] * (e - c0))) / ((e - c0) ** 2) - (Qe[0] + e * Qe[1])
        al = [Dal / 2, e * Dal / 2]
    if PINS:
        t = {}
        for eps in PINS:
            Qe = [sum((R[j][i] * (eps ** j) for j in range(len(R))), fmpq(0)) for i in range(2)]   # T_j(+-1) = (+-1)^j
            t[eps] = [(pinv[eps][i] - d0[i] - d1[i] * (eps - c0)) / ((eps - c0) ** 2) - Qe[i] for i in range(2)]
        if len(PINS) == 1: (eps, tv), = t.items(); al = tv
        else: al = [(t[1][i] + t[-1][i]) / 2 for i in range(2)]; be = [(t[1][i] - t[-1][i]) / 2 for i in range(2)]
    return dict(R=R, al=al, be=be, pinv=pinv)
def cheb_vals(c, n):
    T = [1, c]
    while len(T) < n + 1: T.append(2 * c * T[-1] - T[-2])
    return T[:n + 1]
def dd_vals(c, n, eps):
    """D_j(c) = (T_j(c) - T_j(eps)) / (c - eps), eps = +-1:  D_0 = 0, D_1 = 1, D_{j+1} = 2c D_j + 2 eps^j - D_{j-1}."""
    Dv = [0, 1]
    while len(Dv) < n + 1: j = len(Dv) - 1; Dv.append(2 * c * Dv[-1] + 2 * (eps ** j) - Dv[-2])
    return Dv[:n + 1]
def deltas(dat, P, c):
    d0, d1, c0, R = dat['d0'], dat['d1'], dat['c0'], P['R']; n = len(R) - 1; T = cheb_vals(c, n)
    Pc = [sum((R[j][i] * T[j] for j in range(n + 1)), 0) + P['al'][i] + P['be'][i] * c for i in range(2)]
    t = c - c0
    dl = [d0[i] + d1[i] * t + t * t * Pc[i] for i in range(2)]
    rho = {}
    for eps in PINS:
        Dv = dd_vals(c, n, eps)
        DP = [sum((R[j][i] * Dv[j] for j in range(n + 1)), 0) + P['be'][i] for i in range(2)]
        rho[eps] = [d1[i] + (c + eps - 2 * c0) * Pc[i] + (eps - c0) ** 2 * DP[i] for i in range(2)]
    return dl, rho
FORCE_MANUAL = False   # testing only
def sqrt_nn(x):
    """enclosure of sqrt(x) for a quantity known to be >= 0."""
    if x > 0: return x.sqrt()
    up = x.upper(); return arb(0).union(up.sqrt()) if up >= 0 else arb('nan')
def pinned_term(rho, c, eps):
    """(1 - eps c) ||rho||_{G_c}.  For a c-series whose direct enclosure is non-finite (||rho|| ~ 0 at the pinned end),
    the value and first derivative are enclosed by hand (2nd derivative left as NaN, so T2 is never used there):
    with N = ||rho||, p = rho0 + rho1, m = rho0 - rho1, N^2 = (1+c)/2 p^2 + (1-c)/2 m^2, so
    N >= sqrt((1+c)/2)|p| and N >= sqrt((1-c)/2)|m|;  d/dc[(1-eps c)N] = -eps N + (1-eps c) N' with
    N N' = (p^2 - m^2)/4 + (1+c)/2 p p' + (1-c)/2 m m'.  For eps = +1:
    (1-c)N' = X [(p^2-m^2)/4 + (1-c)/2 m m'] + (1-c) Y p',  X = (1-c)/N in [0, sqrt(2(1-c))/|m|],
    Y = (1+c)/2 p/N in [-sqrt((1+c)/2), sqrt((1+c)/2)].  (eps = -1 symmetric: swap the roles of p and m.)"""
    direct = (1 - eps * c) * gnorm(rho, c)
    if not isinstance(c, arb_series): return direct
    cs = coeffs(direct, 3)
    cv = coeffs(c, 1)[0]
    r0 = [coeffs(z, 2) for z in rho]
    p, dp = r0[0][0] + r0[1][0], r0[0][1] + r0[1][1]; m, dm = r0[0][0] - r0[1][0], r0[0][1] - r0[1][1]
    # work in the frame c' = eps c (so the pinned end is c' = 1): N^2 = (1+c')/2 P^2 + (1-c')/2 M^2 with
    # (P, M) = (p, m) for eps = 1 and (m, p) for eps = -1; derivatives d/dc' = eps d/dc.
    if eps == 1: P, dP, M, dM, cp_ = p, dp, m, dm, cv
    else: P, dP, M, dM, cp_ = m, -dm, p, -dp, -cv
    N = sqrt_nn((1 + cp_) * P * P / 2 + (1 - cp_) * M * M / 2)
    if not (abs(M) > 0): return direct
    X = arb(0).union(sqrt_nn(2 * (1 - cp_)) / abs(M)); h = sqrt_nn((1 + cp_) / 2); Y = (-h).union(h)
    dNp = X * ((P * P - M * M) / 4 + (1 - cp_) / 2 * M * dM) + (1 - cp_) * Y * dP     # = (1 - c') dN/dc'
    val = (1 - cp_) * N; d1 = eps * (-N + dNp)                                           # d/dc = eps d/dc'
    pick = lambda u, v: v if FORCE_MANUAL or not u.is_finite() or (v.is_finite() and v.rad() < u.rad()) else u
    return arb_series([pick(cs[0], val), pick(cs[1], d1), cs[2]], prec=3)
def gnorm(w, c):
    p, m = w[0] + w[1], w[0] - w[1]   # ||w||^2 = (1+c)/2 (w0+w1)^2 + (1-c)/2 (w0-w1)^2 : both terms >= 0 on [-1,1]
    qf = (1 + c) * p * p / 2 + (1 - c) * m * m / 2
    if isinstance(qf, arb) and not qf > 0:   # true value >= 0 (sum of squares, |c| <= 1): enclose sqrt by hull(0, sqrt(upper))
        up = qf.upper()
        return arb(0).union(up.sqrt()) if up >= 0 else qf.sqrt()
    return qf.sqrt()
NTAY = None
class TP:
    """Polynomial in t = c - cm with generic coefficients (arb, arb_series in x, fmpq, int).  Exact algebra: delta and rho
    are polynomials in c, so their expansion about the exact point cm has no remainder.  (Using arb_series for c would
    wrongly identify the c-variable with the x-variable of x-series data.)"""
    def __init__(self, cs): self.c = list(cs)
    def _w(o): return o if isinstance(o, TP) else TP([o])
    def __add__(a, b):
        b = TP._w(b); n = max(len(a.c), len(b.c))
        return TP([(a.c[i] if i < len(a.c) else 0) + (b.c[i] if i < len(b.c) else 0) for i in range(n)])
    __radd__ = __add__
    def __neg__(a): return TP([-x for x in a.c])
    def __sub__(a, b): return a + (-TP._w(b))
    def __rsub__(a, b): return TP._w(b) + (-a)
    def __mul__(a, b):
        b = TP._w(b); out = [0] * (len(a.c) + len(b.c) - 1)
        for i, x in enumerate(a.c):
            for j, y in enumerate(b.c): out[i + j] = out[i + j] + x * y
        return TP(out)
    __rmul__ = __mul__
def deltas_centered(dat, P, cm, c):
    """delta and rho expanded exactly about the exact point cm (TP algebra), evaluated by Horner in t = c - cm."""
    dl0, rho0 = deltas(dat, P, TP([arb(cm), 1]))
    t = c - arb(cm)
    def ev(z):
        cs = z.c if isinstance(z, TP) else [z]; acc = cs[-1]
        for a in reversed(cs[:-1]): acc = acc * t + a
        return acc
    return [ev(z) for z in dl0], {e: [ev(z) for z in v] for e, v in rho0.items()}
def F(dat, P, c, cm=None):
    g = dat['g']
    dl, rho = deltas(dat, P, c) if cm is None else deltas_centered(dat, P, cm, c)
    out = []
    for s0, s1 in COMBOS:
        u = [g[2 + i] / 2 + dl[i] + s0 * g[4 + i] for i in range(2)]
        w = [g[2 + i] / 2 - dl[i] + s1 * g[6 + i] for i in range(2)]
        sA = sB = None
        for eps, (bx, ba) in PINS.items():
            if bx == 0 and ba == s0: sA = pinned_term(rho[eps], c, eps)
            if bx == 1 and ba == s1: sB = pinned_term(rho[eps], c, eps)
        if sA is None: sA = gnorm(u, c)
        if sB is None: sB = gnorm(w, c)
        out.append(1 - s0 * g[0] - s1 * g[1] - sA - sB)
    return out
def _midq(H):
    m = H.mid(); return fmpq(int((m * 2 ** 60).floor().unique_fmpz()), 2 ** 60)
def ball(a, b): return arb((a + b) / 2, (b - a) / 2)
def coeffs(z, n):
    cs = z.coeffs() if not isinstance(z, arb) else [z]
    return list(cs) + [arb(0)] * (n - len(cs))
def endpoint_ok(dat, poly, eps):
    g = dat['g']; e = arb(eps); dl, _ = deltas(dat, poly, e)
    for s0, s1 in COMBOS:
        pinned = any((bx == 0 and ba == s0) or (bx == 1 and ba == s1) for ep, (bx, ba) in PINS.items() if ep == eps)
        if pinned: continue
        u = [g[2 + i] / 2 + dl[i] + s0 * g[4 + i] for i in range(2)]; w = [g[2 + i] / 2 - dl[i] + s1 * g[6 + i] for i in range(2)]
        uf = u[0] + eps * u[1]; wf = w[0] + eps * w[1]; T = 1 - s0 * g[0] - s1 * g[1]
        for et in (1, -1):
            if not (T - et * uf + et * wf > 0): return False
    return True
def accept(ca, cb, xa, xb, cache, R):
    key = (xa, xb)
    if key not in cache:
        dB, dS, dM = theta_data_mv(xa, xb)
        cache[key] = ((dB, build(dB, R)), (dS, build(dS, R)), (dM, build(dM, R)))
    (dB, pB), (dS, pS), (dM, pM) = cache[key]
    cm, rc, rx = (ca + cb) / 2, (cb - ca) / 2, (xb - xa) / 2; cb_ = ball(ca, cb)
    Fm = F(dM, pM, arb(cm)); Fc = F(dB, pB, arb_series([cb_, 1], prec=3), cm); Fx = F(dS, pS, cb_, cm)
    F0 = F(dB, pB, cb_, cm)   # plain (non-series) enclosure: finite even where a norm vanishes
    tests = []; ec = ex = 0.0; hc = {}
    for i in range(4):
        fc = coeffs(Fc[i], 3); fx = coeffs(Fx[i], 2); d1c, dx = fc[1], fx[1]
        if F0[i].is_finite() and F0[i] > 0: tests.append('T0'); continue
        if all(z.is_finite() for z in (Fm[i], d1c, dx)):
            if Fm[i] - abs(d1c) * arb(rc) - abs(dx) * arb(rx) > 0: tests.append('T1'); continue
            ec = max(ec, float(abs(d1c).upper()) * float(rc)); ex = max(ex, float(abs(dx).upper()) * float(rx))
        else: ec = float('nan')
        if d1c.is_finite() and cb == 1 and d1c < 0: tests.append('T3+'); continue
        if d1c.is_finite() and ca == -1 and d1c > 0: tests.append('T3-'); continue
        if cb > fmpq(3, 4) and ca < 1:   # T3h: F_c < 0 on hull(box, 1) and F(1) >= 0  =>  F >= 0 on the box
            if '_h1' not in hc: hc['_h1'] = F(dB, pB, arb_series([ball(ca, fmpq(1)), 1], prec=3), _midq(ball(ca, fmpq(1))))
            z = coeffs(hc['_h1'][i], 3)[1]
            if z.is_finite() and z < 0: tests.append('T3h+'); continue
        if ca < -fmpq(3, 4) and cb > -1:
            if '_h-1' not in hc: hc['_h-1'] = F(dB, pB, arb_series([ball(fmpq(-1), cb), 1], prec=3), _midq(ball(fmpq(-1), cb)))
            z = coeffs(hc['_h-1'][i], 3)[1]
            if z.is_finite() and z > 0: tests.append('T3h-'); continue
        H = cb_.union(dB['c0']); cmH = fmpq(int(H.mid().man_exp()[0]) if False else 0) if False else None
        cmH = _midq(H)
        if H.rad() < 0.3:
            f2 = coeffs(F(dB, pB, arb_series([H, 1], prec=3), cmH)[i], 3)[2]
            if f2.is_finite() and f2 > 0: tests.append('T2'); continue
        if ec != ec: return ('split', 'r')
        return ('split', 'x' if ex > ec else 'c')
    return tests
def run_strip(j):
    S = CERT['strips'][j]; xa0, xb0 = q(S['xa']), q(S['xb']); R = [[q(a), q(b)] for a, b in S['R']]
    # endpoint check on x-subboxes
    stack = [(xa0 + (xb0 - xa0) * fmpq(i, 8), xa0 + (xb0 - xa0) * fmpq(i + 1, 8)) for i in range(8)]; nep = 0
    while stack:
        a, b = stack.pop(); dat = theta_data(ball(a, b)); poly = build(dat, R)
        if endpoint_ok(dat, poly, 1) and endpoint_ok(dat, poly, -1): nep += 1; continue
        if b - a < (xb0 - xa0) / 2 ** 16: return dict(ok=False, why='endpoint', x=[str(a), str(b)])
        m = (a + b) / 2; stack += [(a, m), (m, b)]
    nc = args.nc; stack = [(fmpq(-1) + fmpq(2 * i, nc), fmpq(-1) + fmpq(2 * (i + 1), nc), xa0, xb0) for i in range(nc)]
    cache = {}; n = 0; area = fmpq(0); cnt = {}
    while stack:
        ca, cb, a, b = stack.pop(); t = accept(ca, cb, a, b, cache, R)
        if t is not None and t[0] != 'split':
            n += 1; area += (cb - ca) * (b - a)
            for x in t: cnt[x] = cnt.get(x, 0) + 1
            if n > args.maxbox: return dict(ok=False, why='too many boxes', n=n)
            continue
        wc, wx = cb - ca, b - a
        nrej = cache.get('_nrej', 0) + 1; cache['_nrej'] = nrej
        if VERB and nrej % 2000 == 0: print('  strip %d: %d accepted, %d split; now at c=[%.6g,%.6g] x=[%.8g,%.8g]' % (j, n, nrej, float(ca), float(cb), float(a), float(b)), flush=True)
        if wc < fmpq(1, 2 ** 30) and wx < (xb0 - xa0) / 2 ** 24: return dict(ok=False, why='box', box=[str(ca), str(cb), str(a), str(b)])
        if t[1] == 'r': sc = wc * (xb0 - xa0) > wx * 2
        else: sc = t[1] == 'c'
        if (sc and wc > fmpq(1, 2 ** 30)) or wx < (xb0 - xa0) / 2 ** 24: m = (ca + cb) / 2; stack += [(ca, m, a, b), (m, cb, a, b)]
        else: m = (a + b) / 2; stack += [(ca, cb, a, m), (ca, cb, m, b)]
        if len(cache) > 300: k_ = cache.get('_nrej', 0); cache.clear(); cache['_nrej'] = k_
    return dict(ok=(area == 2 * (xb0 - xa0)), nboxes=n, endpoint_boxes=nep, tests=cnt)
def lambda_run(kj): return kj[0], run_strip(kj[1])
if __name__ == '__main__':
    t0 = time.time(); js = list(range(*(args.strips or [0, len(CERT['strips'])])))
    print('vgen %s vertex=%s pins=%s prec=%d strips %d..%d' % (args.cert, CERT['vertex'], PINS, args.prec, js[0], js[-1]), flush=True)
    if args.jobs > 1:
        from multiprocessing import Pool
        res = [None] * len(js)
        with Pool(args.jobs) as P:
            for k, r in P.imap_unordered(lambda_run, list(enumerate(js))):
                res[k] = r; print('  done strip %d: ok=%s %s' % (js[k], r['ok'], {a: b for a, b in r.items() if a != 'tests'}), flush=True)
    else: res = [run_strip(j) for j in js]
    allok = True; tot = 0
    for j, r in zip(js, res):
        S = CERT['strips'][j]; allok &= r['ok']; tot += r.get('nboxes', 0)
        print('strip %3d x in [%s, %s] (theta/pi ~ [%.6f, %.6f]): %s' % (j, S['xa'], S['xb'], float(q(S['xa'])), float(q(S['xb'])), r), flush=True)
    print('RESULT %s: ALL OK = %s ; %d boxes ; %.1fs' % (args.cert, allok, tot, time.time() - t0), flush=True)
