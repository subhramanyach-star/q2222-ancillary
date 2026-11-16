#!/usr/bin/env python3
"""tgen3.py (= tgen.py + Taylor endpoint check near t=0 + mean-value edge bound in E-) -- small-theta verifier for the V3 global certificate (smalltheta project).

CLAIM. For t = s^2 = sin^2(2 theta) in [0, t0] (cert 'tmax') and every c in [-1,1] and every sign pair (s0,s1):
   G_{s0 s1}(c,t) >= 0,   where for t > 0   G = F/t,  F = 1 - s0 g0 - s1 g1 - ||u_{s0}||_c - ||w_{s1}||_c  (vgen.py's F),
so Phi_c(delta(c,t)) <= 1 = g.p_theta for all theta with 0 < sin^2 2theta <= t0.
DATA. Exact hats from prep.py: X(t) = X(0) + t X^(t), X^ = q(t,gam,rho)/D(t,gam,rho), gam = 1/(1+sqrt(1-t)), rho = 1/(1+sqrt(1+t))
(exact identities sqrt(1-t) = 1 - t gam, sqrt(1+t) = 1 + t rho). V3 facts used (asserted exactly below):
g_a(0) = 0 (a = Alice marginals), d1(0) = 0, u_{s0}(c,0) = b(0)/2 + d0(0) + s0 e0(0) = 0, w00_{s1} = b(0)/2 - d0(0) + s1 e1(0) in {(1,0),(0,1)}.
CERTIFICATE. delta = d0(0) + t dhat, dhat = d0^ + d1^ (c-c0) + (c-c0)^2 [R(c) + al],  R(c) = R0 + (1+c) sum_j Rt_j T_j(c),
al = (Dal/2, e Dal/2) the exact mid-pin at e = -1 (D(delta(e)) = (K_i+K_j)/2, K = midpin.cands, linear in g).
Then  G = -(s0 g0^ + s1 g1^) - ||u^||_c - Wq,  u^ = b^/2 + dhat + s0 e0^,  w^ = b^/2 - dhat + s1 e1^,
Wq = (||w||_c - 1)/t = (2<w00,w^>_c + t ||w^||_c^2)/(||w||_c + 1),  w = w00 + t w^.   (Algebraic identities; no cancellation.)
TESTS per box [ca,cb]x[ta,tb] (Arb): T0 direct; T1 mean value (G_c, G_t enclosures); T2 G_cc > 0 on hull(box, c0(tbox))
[valid: F(c0)=F_c(c0)=0 exactly for t>0 (fullrange/exact_c0.log), hence G too; at t=0 (c0=1) by continuity];
T3/T3h monotone in c towards an endpoint where G >= 0 (endpoint check); CORNER (box touching (c,t)=(-1,0)): see corner_dc().
ENDPOINT CHECK at c = +-1 per t-box: G(eps) = min_{eta,eta'} [1 - eta' w00.f]/t + [-(s0 g0^+s1 g1^) - eta u^.f - eta' w^.f];
terms with eta = eta' equal (1 - g.L)/t >= 0 (comb.py, exact, 0 < s < s*); mixed terms must be > 0 in Arb.
"""
import json, sys, time, argparse
import sympy as sp, flint
from flint import arb, arb_series, fmpq
sys.path.insert(0, '../fullrange')
from midpin import cands
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('--prec', type=int, default=100)
ap.add_argument('--nc', type=int, default=64); ap.add_argument('--nt', type=int, default=32); ap.add_argument('--maxbox', type=int, default=4000000)
ap.add_argument('--tstrips', nargs=2, type=int, default=None); ap.add_argument('--jobs', type=int, default=1)
args = ap.parse_args(); flint.ctx.prec = args.prec
CERT = json.load(open(args.cert)); import os; H = json.load(open(os.path.basename(CERT['hats'])))[CERT['vertex']]  # bundle: hats file resolved in the working directory (smalltheta/)
t_, gam_, rho_ = sp.symbols('t gam rho')
def comp(expr_str):
    e = sp.sympify(expr_str, locals={'t': t_, 'gam': gam_, 'rho': rho_})
    def b(e):
        if e.is_Symbol: i = {t_: 0, gam_: 1, rho_: 2}[e]; return lambda v: v[i]
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
def Q(x): x = str(x); return fmpq(*map(int, x.split('/'))) if '/' in x else fmpq(int(x))
X0 = {k: [Q(a) for a, _, _ in H[k]] for k in ('g', 'd0', 'd1')}
XH = {k: [(comp(q), comp(D)) for _, q, D in H[k]] for k in ('g', 'd0', 'd1')}
# exact structural facts (V3)
g0, d00 = X0['g'], X0['d0']
assert g0[0] == 0 and g0[1] == 0 and X0['d1'] == [0, 0]
for s0 in (1, -1): assert all(g0[2 + i] / 2 + d00[i] + s0 * g0[4 + i] == 0 for i in range(2))
W00 = {s1: [g0[2 + i] / 2 - d00[i] + s1 * g0[6 + i] for i in range(2)] for s1 in (1, -1)}
assert W00[1] == [1, 0] and W00[-1] == [0, 1]
E, PI_, PJ_ = CERT['mpin']; assert E == -1
K00 = cands(g0, E); assert (K00[PI_] + K00[PJ_]) / 2 == d00[0] + E * d00[1]
R0 = [Q(a) for a in CERT['R0']]; RT = [[Q(a), Q(b)] for a, b in CERT['Rt']]
def hats(t):
    """t: arb or arb_series (in t). returns hatted data and c0."""
    if isinstance(t, arb_series): sq1 = (1 - t).sqrt() if False else None
    one = t * 0 + 1
    gam = 1 / (1 + (1 - t).sqrt()); rho = 1 / (1 + (1 + t).sqrt()); v = (t, gam, rho)
    d = {k: [q(v) / D(v) for q, D in XH[k]] for k in XH}
    d['c0'] = (1 - t) / (1 + t); d['t'] = t
    return d
def hats_mv(ta, tb):
    """mean-value enclosure of the hatted data over [ta,tb]: X(tm) + X'([ta,tb]) [-rt, rt]  (also returns the t-series)."""
    tB = ball(ta, tb); rt = arb(0, (tb - ta) / 2)
    dM = hats(arb((ta + tb) / 2)); dS = hats(arb_series([tB, 1], prec=2)); out = {}
    for k, v in dM.items():
        if k == 't': out[k] = tB; continue
        if isinstance(v, list): out[k] = [m + coeffs(z, 2)[1] * rt for m, z in zip(v, dS[k])]
        else: out[k] = v + coeffs(dS[k], 2)[1] * rt
    return out, dS, dM
def cheb(c, n):
    T = [c * 0 + 1, c]
    while len(T) < n + 1: T.append(2 * c * T[-1] - T[-2])
    return T[:n + 1]
# R(c) = R0 + (1+c) sum_j Rt_j T_j(c) converted EXACTLY to monomials in c (fmpq), then evaluated about an exact
# centre cm by an exact Taylor shift (no wrapping from the Chebyshev recurrence).
def _cheb_mono(n):
    T = [[fmpq(1)], [fmpq(0), fmpq(1)]]
    while len(T) < n + 1:
        a, b = T[-1], T[-2]; t = [fmpq(0)] + [2 * v for v in a]
        for i, v in enumerate(b): t[i] -= v
        T.append(t)
    return T[:n + 1]
_TM = _cheb_mono(len(RT) - 1)
RMONO = []
for i in range(2):
    m = [fmpq(0)] * (len(RT) + 1)
    for j in range(len(RT)):
        for k, v in enumerate(_TM[j]): m[k] += RT[j][i] * v; m[k + 1] += RT[j][i] * v   # (1+c) T_j
    m[0] += R0[i]; RMONO.append(m)
def _shift(m, cm):
    """coefficients of p(cm + x) in x, exact."""
    n = len(m); out = [fmpq(0)] * n
    from math import comb
    for k in range(n):
        for j in range(k, n): out[k] += m[j] * comb(j, k) * cm ** (j - k)
    return out
_SH = {}
def Rpoly(c, cm=None):
    if cm is None:
        return [sum((RMONO[i][k] * c ** k for k in range(len(RMONO[i]))), c * 0) for i in range(2)]
    if cm not in _SH: _SH[cm] = [_shift(RMONO[i], cm) for i in range(2)]
    x = c - cm; out = []
    for cs in _SH[cm]:
        acc = c * 0 + cs[-1]
        for a in reversed(cs[:-1]): acc = acc * x + a
        out.append(acc)
    return out
def pin(d):
    gh, d0h, d1h, c0 = d['g'], d['d0'], d['d1'], d['c0']; e = E
    K = cands(gh, e); tgt = (K[PI_] + K[PJ_]) / 2
    Re = Rpoly(arb(e))
    Dal = (tgt - (d0h[0] + d1h[0] * (e - c0)) - e * (d0h[1] + d1h[1] * (e - c0))) / ((e - c0) ** 2) - (Re[0] + e * Re[1])
    return [Dal / 2, e * Dal / 2]
def ip(a, b, c): return a[0] * b[0] + a[1] * b[1] + c * (a[0] * b[1] + a[1] * b[0])
def nsq(a, c):
    p, m = a[0] + a[1], a[0] - a[1]; return (1 + c) * p * p / 2 + (1 - c) * m * m / 2
def sqrt_nn(x):
    if isinstance(x, arb) and not x > 0:
        up = x.upper(); return arb(0).union(up.sqrt()) if up >= 0 else arb('nan')
    return x.sqrt()
def vecs(d, c, al, cm=None):
    gh, d0h, d1h, c0 = d['g'], d['d0'], d['d1'], d['c0']
    Rc = Rpoly(c, cm); tc = c - c0 if cm is None else (c - cm) + (cm - c0)
    dh = [d0h[i] + d1h[i] * tc + tc * tc * (Rc[i] + al[i]) for i in range(2)]
    U = {s0: [gh[2 + i] / 2 + dh[i] + s0 * gh[4 + i] for i in range(2)] for s0 in (1, -1)}
    Wh = {s1: [gh[2 + i] / 2 - dh[i] + s1 * gh[6 + i] for i in range(2)] for s1 in (1, -1)}
    return U, Wh
def G(d, c, al=None, cm=None, no_u=False):
    if al is None: al = pin(d)
    gh, t = d['g'], d['t']; U, Wh = vecs(d, c, al, cm); out = []
    for s0 in (1, -1):
        nu = 0 if no_u else sqrt_nn(nsq(U[s0], c))
        for s1 in (1, -1):
            w00 = W00[s1]; wh = Wh[s1]
            w = [w00[i] + t * wh[i] for i in range(2)]
            Wq = (2 * ip(w00, wh, c) + t * nsq(wh, c)) / (sqrt_nn(nsq(w, c)) + 1)
            out.append(-(s0 * gh[0] + s1 * gh[1]) - nu - Wq)
    return out
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
# ---------------------------------------------------------------------------------------------------------------
# ENDPOINT CHECK (tgen2).  At c = eps the slice norm is |v0 + eps v1|, so F(eps) = 1 - s0 a0 - s1 a1 - |u_f| - |w_f|.
# u_f = t uf, w_f = wfull = w00f + t wf (u00 = 0 for V3).  Sign choices (eta, eta') for |u_f|, |w_f| are dropped when the
# sign of uf / wfull is certified on the t-box (for t > 0); same-sign terms equal (1 - g.L)/t >= 0 by comb.py
# (deterministic L, all 16 margins >= 0 on (0, pi/4]); mixed terms val = const + t br (const = 1 - eta' w00f) must be > 0.
# Each scalar X in {uf, wfull, val} is an analytic function of t on [0, 1/63] (hats are polynomials in t, gam, rho over
# nonvanishing denominators).  Its sign on a box [ta, tb] (t > 0) is certified either
#   (direct) by an Arb enclosure over the t-ball, or
#   (Taylor) by exact rational Taylor coefficients x_0..x_{M-1} at t = 0 (fmpq_series, exact arithmetic) and the Lagrange
#            remainder: X(t) = sum_{j<M} x_j t^j + C_M(xi) t^M, xi in [0, tb], C_M enclosed by the M-th coefficient of the
#            Arb Taylor series of X over the ball [0, tb].  If x_k is the first nonzero coefficient (k < M), then
#            X(t)/t^k = sum_{k<=j<M} x_j t^(j-k) + C_M t^(M-k), enclosed over t in [ta, tb]; its sign is the sign of X on (0,tb].
from flint import fmpq_series
MT = 5
def _comp_generic(expr_str, const):
    e = sp.sympify(expr_str, locals={'t': t_, 'gam': gam_, 'rho': rho_})
    def b(e):
        if e.is_Symbol: i = {t_: 0, gam_: 1, rho_: 2}[e]; return lambda v: v[i]
        if e.is_Rational: q = fmpq(int(e.p), int(e.q)); return lambda v: const(q)
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
_XQ = {k: [(_comp_generic(q, lambda x: fmpq_series([x])), _comp_generic(D, lambda x: fmpq_series([x]))) for _, q, D in H[k]] for k in ('g', 'd0', 'd1')}
def hats_exact():
    flint.ctx.cap = MT + 2
    t = fmpq_series([0, 1]); gam = 1 / (1 + (1 - t).sqrt()); rho = 1 / (1 + (1 + t).sqrt()); v = (t, gam, rho)
    d = {k: [q(v) / D(v) for q, D in _XQ[k]] for k in _XQ}
    d['c0'] = (1 - t) / (1 + t); d['t'] = t
    return d
def pin_generic(d, Re):
    gh, d0h, d1h, c0 = d['g'], d['d0'], d['d1'], d['c0']; e = E
    K = cands(gh, e); tgt = (K[PI_] + K[PJ_]) / 2
    Dal = (tgt - (d0h[0] + d1h[0] * (e - c0)) - e * (d0h[1] + d1h[1] * (e - c0))) / ((e - c0) ** 2) - (Re[0] + e * Re[1])
    return [Dal / 2, e * Dal / 2]
def ep_quantities(d, al, cconst):
    """dict (eps, k, name) -> scalar; name in 'uf', 'wfull', ('val', eta, etp)."""
    gh, t = d['g'], d['t']; out = {}
    for eps in (1, -1):
        U, Wh = vecs(d, cconst(eps), al)
        for k, (s0, s1) in enumerate(COMBOS):
            uf = U[s0][0] + eps * U[s0][1]; wf = Wh[s1][0] + eps * Wh[s1][1]; w00f = W00[s1][0] + eps * W00[s1][1]
            out[(eps, k, 'uf')] = uf; out[(eps, k, 'wfull')] = w00f + t * wf
            for eta in (1, -1):
                for etp in (1, -1):
                    if eta == etp: continue
                    const = 1 - etp * w00f; br = -(s0 * gh[0] + s1 * gh[1]) - eta * uf - etp * wf
                    out[(eps, k, ('val', eta, etp))] = br if const == 0 else const + t * br
    return out
def _exact_coeffs():
    d = hats_exact(); Re = Rpoly(fmpq(E)); al = pin_generic(d, Re)
    q = ep_quantities(d, al, lambda e: fmpq(e)); out = {}
    for key, x in q.items():
        if not isinstance(x, fmpq_series): x = fmpq_series([x])
        cs = list(x.coeffs()); cs = cs + [fmpq(0)] * (MT - len(cs)); out[key] = cs[:MT]
    flint.ctx.cap = 64
    return out
EXC = _exact_coeffs()
_REMC = {}
def _remainders(tb):
    if tb not in _REMC:
        dS = hats(arb_series([ball(fmpq(0), tb), 1], prec=MT + 1)); al = pin(dS)
        q = ep_quantities(dS, al, lambda e: arb(e)); out = {}
        for key, x in q.items():
            cs = coeffs(x, MT + 1) if isinstance(x, arb_series) else [x] + [arb(0)] * MT
            out[key] = cs[MT]
        _REMC[tb] = out
        if len(_REMC) > 64: _REMC.clear()
    return _REMC[tb]
def _sign(key, qB, ta, tb, tB):
    """+1 / -1 if the sign of X on the box (t > 0) is certified, else 0."""
    x = qB[key]
    if x > 0: return 1
    if x < 0: return -1
    if tb > fmpq(1, 64): return 0
    cs = EXC[key]; k = next((j for j in range(MT) if cs[j] != 0), None)
    if k is None: return 0
    acc = arb(0); tp = arb(1)
    for j in range(k, MT): acc += arb(cs[j]) * tp; tp = tp * tB
    acc += _remainders(tb)[key] * tp
    if acc > 0: return 1
    if acc < 0: return -1
    return 0
def endpoint_ok(d, eps, ta=None, tb=None):
    al = pin(d); qB = ep_quantities(d, al, lambda e: arb(e)); tB = d['t']
    for k, (s0, s1) in enumerate(COMBOS):
        su = _sign((eps, k, 'uf'), qB, ta, tb, tB); sw = _sign((eps, k, 'wfull'), qB, ta, tb, tB)
        etas = [su] if su else [1, -1]; etps = [sw] if sw else [1, -1]
        for eta in etas:
            for etp in etps:
                if etp == eta: continue
                if _sign((eps, k, ('val', eta, etp)), qB, ta, tb, tB) != 1: return False
    return True
def ball(a, b): return arb((a + b) / 2, (b - a) / 2)
def coeffs(z, n):
    cs = z.coeffs() if isinstance(z, arb_series) else [z]
    cs = list(cs); return cs + [arb(0)] * (n - len(cs))
def midq(x):
    m = x.mid(); return fmpq(int((m * 2 ** 60).floor().unique_fmpz()), 2 ** 60)
from math import comb as _comb
KT = 4
def ctaylor(dat, al, ca, cb, K=KT, jmax=2):
    """Enclosures, valid for all c in [ca,cb] (and all t in the data's t-range), of the Taylor coefficients
    G^(j)(c)/j!, j = 0,1,2, per pair:  sum_{k=j}^{K-1} C(k,j) a_k(cm) x^(k-j) + C(K,j) A_K x^(K-j),  x = c - cm in [-r,r],
    with a_k the Taylor coefficients at the exact centre cm and A_K = G^(K)(xi)/K! enclosed over the ball (Taylor's theorem)."""
    cm = (ca + cb) / 2; r = (cb - ca) / 2
    P = G(dat, arb_series([arb(cm), 1], prec=K + 1), al, cm)
    B = G(dat, arb_series([ball(ca, cb), 1], prec=K + 1), al, cm)
    x = arb(0, r); out = []
    xp = [arb(1)]
    for _ in range(K): xp.append(xp[-1] * x)
    for i in range(4):
        a_ = coeffs(P[i], K + 1); AK = coeffs(B[i], K + 1)[K]
        enc = []
        for j in range(jmax + 1):
            acc = _comb(K, j) * AK * xp[K - j]
            for k in range(j, K): acc = acc + _comb(k, j) * a_[k] * xp[k - j]
            enc.append(acc)
        out.append(enc)
    return out
def accept(ca, cb, ta, tb, cache):
    key = (ta, tb)
    if key not in cache:
        dB, dS, dM = hats_mv(ta, tb)
        cache[key] = (dB, pin(dB), dS, pin(dS), dM, pin(dM))
    dB, aB, dS, aS, dM, aM = cache[key]
    cm, rc, rt = (ca + cb) / 2, (cb - ca) / 2, (tb - ta) / 2; cB = ball(ca, cb)
    EB = ctaylor(dB, aB, ca, cb)              # over the whole box
    tests = []; ec = et = 0.0; EM = None; Gt = None; Hc = {}
    for i in range(4):
        v0 = EB[i][0]
        if v0.is_finite() and v0 > 0: tests.append('T0'); continue
        if 'direct' not in Hc: Hc['direct'] = G(dB, cB, aB, cm)    # plain ball evaluation (no c-derivatives)
        if Hc['direct'][i].is_finite() and Hc['direct'][i] > 0: tests.append('T0d'); continue
        if EM is None: EM = ctaylor(dM, aM, ca, cb); Gt = G(dS, cB, aS, cm)
        vm = EM[i][0]; dt = coeffs(Gt[i], 2)[1]   # G(c,t) >= G(c,tm) - |G_t| rt
        if vm.is_finite() and dt.is_finite():
            if vm - abs(dt) * arb(rt) > 0: tests.append('T1'); continue
            et = max(et, float(abs(dt).upper()) * float(rt))
        d1c = EB[i][1]
        ec = max(ec, float(abs(d1c).upper()) * float(rc)) if d1c.is_finite() else float('nan')
        if not (dt.is_finite()): et = float('inf')
        if d1c.is_finite() and cb == 1 and d1c < 0: tests.append('T3+'); continue
        if d1c.is_finite() and ca == -1 and d1c > 0: tests.append('T3-'); continue
        if cb > fmpq(3, 4) and ca < 1:
            if 'h+' not in Hc: Hc['h+'] = ctaylor(dB, aB, ca, fmpq(1))
            z = Hc['h+'][i][1]
            if z.is_finite() and z < 0: tests.append('T3h+'); continue
        if ca < -fmpq(3, 4) and cb > -1:
            if 'h-' not in Hc: Hc['h-'] = ctaylor(dB, aB, fmpq(-1), cb)
            z = Hc['h-'][i][1]
            if z.is_finite() and z > 0: tests.append('T3h-'); continue
        if ca < -fmpq(3, 4) and COMBOS[i][0] == -1:   # G_c > 0 on [-1, cb] x tbox and G(-1) >= 0  =>  G >= 0 on the box
            if 'edge' not in Hc: Hc['edge'] = edge_dc(dB, aB, dS, aS, fmpq(-1), cb, ta, tb)
            z = Hc['edge'][i]
            if z is not None and z.is_finite() and z > 0: tests.append('E-'); continue
        if cb > fmpq(1, 2):
            lo = min(ca, fmpq(int((dB['c0'].lower() * 2 ** 40).floor().unique_fmpz()), 2 ** 40))
            if 'T4' not in Hc: Hc['T4'] = ctaylor(dB, aB, lo, fmpq(1), K=5, jmax=3)
            z = Hc['T4'][i][3]
            if z.is_finite() and z < 0: tests.append('T4'); continue
        Hh = cB.union(dB['c0'])
        if Hh.rad() < 0.3:
            if 'T2' not in Hc:
                lo = min(ca, fmpq(int((Hh.lower() * 2 ** 40).floor().unique_fmpz()), 2 ** 40)); hi = max(cb, fmpq(int((Hh.upper() * 2 ** 40).ceil().unique_fmpz()), 2 ** 40))
                Hc['T2'] = ctaylor(dB, aB, lo, min(hi, fmpq(1)))
            f2 = Hc['T2'][i][2]
            if f2.is_finite() and f2 > 0: tests.append('T2'); continue
        if ec != ec: return ('split', 'r', i)
        return ('split', 't' if et > ec else 'c', i)
    return tests
def edge_dc(dB, aB, dS, aS, ca, cb, ta, tb):
    """Lower bound of G_c on the box [ca=-1, cb] x [ta, tb] for the pairs with s0 = -1, robust where u^_- is small.
    G = A - N with N = ||v||_c, v = u^_{-}; N^2 = (1+c)/2 P^2 + (1-c)/2 M^2 (P = v0+v1, M = v0-v1), so
    N' = (P^2 - M^2)/(4N) + [(1+c) P P_c + (1-c) M M_c]/(2N)
       <= P^2/(4 |M| sqrt((1-c)/2)) + sqrt((1+c)/2)|P_c| + sqrt((1-c)/2)|M_c|      (N >= sqrt((1+-c)/2)|P|,|M|).
    P^2/|M|: if ta > 0 directly (needs M bounded away from 0); if ta = 0, by Hadamard at the corner (P = M = 0 at (-1,0),
    asserted): P = (1+c) Pc~ + t Pt~, M = (1+c) Mc~ + t Mt~ with ~ = derivative ranges over the box (the segment from the
    corner lies in the box); if Mc~, Mt~ > 0 then P^2/|M| <= ((1+c)+t) max(|Pc~|,|Pt~|)^2 / min(Mc~, Mt~).
    Returns per pair a lower bound of A_c - N'_upper (None for s0 = +1)."""
    cm = (ca + cb) / 2; cB = ball(ca, cb)
    A = G(dB, arb_series([cB, 1], prec=2), aB, cm, no_u=True)
    U, _ = vecs(dB, arb_series([cB, 1], prec=2), aB, cm); v = U[-1]
    P, M = v[0] + v[1], v[0] - v[1]; Pv, Pc = coeffs(P, 2)[:2]; Mv, Mc = coeffs(M, 2)[:2]
    onep = cB + 1; onem = 1 - cB
    if ta > 0:
        # (a) direct: |M| bounded away from 0 on the box;  (b) tgen3 edge bound: mean value from c = -1 (same t):
        #   P(c) = P(-1) + (1+c) P_c(xi),  M(c) >= M(-1) + (1+c) min M_c  (needs M(-1) > 0, M_c > 0 on [-1,cb] x tbox),
        #   P^2/M <= 2 P(-1)^2/M(-1) + 2 (1+cb) P_c^2 / min M_c   ((a+b)^2 <= 2a^2 + 2b^2;  a^2/(x+y) <= a^2/x).
        Qs = []
        if Mv > 0 or Mv < 0: Qs.append(Pv * Pv / abs(Mv))
        Ue, _ = vecs(dB, arb(-1), aB); Pe = Ue[-1][0] + Ue[-1][1]; Me = Ue[-1][0] - Ue[-1][1]
        if Me > 0 and Mc > 0: Qs.append(2 * Pe * Pe / Me + 2 * (1 + arb(cb)) * Pc * Pc / Mc)
        if not Qs: return [None] * 4
        Q = Qs[0] if len(Qs) == 1 else (Qs[0] if Qs[0].upper() < Qs[1].upper() else Qs[1])
    else:
        Ut, _ = vecs(dS, cB, aS, cm); Pt = coeffs(Ut[-1][0] + Ut[-1][1], 2)[1]; Mt = coeffs(Ut[-1][0] - Ut[-1][1], 2)[1]
        if not (Mc > 0 and Mt > 0): return [None] * 4
        mx = abs(Pc).union(abs(Pt)); mn = Mc.union(Mt)
        Q = (onep + ball(ta, tb)) * mx * mx / mn
    Nup = Q / (4 * sqrt_nn(onem / 2)) + sqrt_nn(onep / 2) * abs(Pc) + sqrt_nn(onem / 2) * abs(Mc)
    out = []
    for k, (s0, s1) in enumerate(COMBOS):
        out.append(coeffs(A[k], 2)[1] - Nup.upper() if s0 == -1 else None)
    return out
CORNER_ZERO_CHECKED = False
def run_tstrip(j):
    S = CERT['tstrips'][j]; ta0, tb0 = Q(S[0]), Q(S[1])
    stack = [(ta0 + (tb0 - ta0) * fmpq(i, 8), ta0 + (tb0 - ta0) * fmpq(i + 1, 8)) for i in range(8)]; nep = 0
    while stack:
        a, b = stack.pop(); d = hats(ball(a, b))
        if endpoint_ok(d, 1, a, b) and endpoint_ok(d, -1, a, b): nep += 1; continue
        if b - a < (tb0 - ta0) / 2 ** 20: return dict(ok=False, why='endpoint', t=[str(a), str(b)])
        m = (a + b) / 2; stack += [(a, m), (m, b)]
    nc = args.nc; NT = args.nt
    stack = [(fmpq(-1) + fmpq(2 * i, nc), fmpq(-1) + fmpq(2 * (i + 1), nc), ta0 + (tb0 - ta0) * fmpq(k, NT), ta0 + (tb0 - ta0) * fmpq(k + 1, NT)) for i in range(nc) for k in range(NT)]
    cache = {}; n = 0; area = fmpq(0); cnt = {}; nsplit = 0; fails = {}
    while stack:
        ca, cb, a, b = stack.pop(); r = accept(ca, cb, a, b, cache)
        if r and r[0] == 'split':
            nsplit += 1; wc, wt = cb - ca, b - a
            if nsplit % 5000 == 0: print('   strip %d: %d acc %d split; at c=[%.8g,%.8g] t=[%.6g,%.6g] pair %d why %s' % (j, n, nsplit, float(ca), float(cb), float(a), float(b), r[2], r[1]), flush=True)
            if wc < fmpq(1, 2 ** 34) and wt < (tb0 - ta0) / 2 ** 30: return dict(ok=False, why='box', box=[str(ca), str(cb), str(a), str(b)], pair=r[2])
            sc = (r[1] == 'c') or (r[1] == 'r' and wc * (tb0 - ta0) > wt * 2)
            if (sc and wc > fmpq(1, 2 ** 34)) or wt < (tb0 - ta0) / 2 ** 30: m = (ca + cb) / 2; stack += [(ca, m, a, b), (m, cb, a, b)]
            else: m = (a + b) / 2; stack += [(ca, cb, a, m), (ca, cb, m, b)]
            if len(cache) > 400: cache.clear()
            continue
        n += 1; area += (cb - ca) * (b - a)
        for x in r: cnt[x] = cnt.get(x, 0) + 1
        if n > args.maxbox: return dict(ok=False, why='too many boxes', n=n)
    return dict(ok=(area == 2 * (tb0 - ta0)), nboxes=n, nsplit=nsplit, endpoint_boxes=nep, tests=cnt)
def lrun(kj): return kj[0], run_tstrip(kj[1])
if __name__ == '__main__':
    t0 = time.time(); js = list(range(*(args.tstrips or [0, len(CERT['tstrips'])])))
    print('tgen3 %s vertex=%s prec=%d tstrips %s' % (args.cert, CERT['vertex'], args.prec, js), flush=True)
    res = {}
    if args.jobs > 1:
        from multiprocessing import Pool
        with Pool(args.jobs) as P:
            for k, r in P.imap_unordered(lrun, list(enumerate(js))): res[js[k]] = r; print('  done tstrip %d: %s' % (js[k], r), flush=True)
    else:
        for j in js: res[j] = run_tstrip(j); print('  done tstrip %d: %s (%.0fs)' % (j, res[j], time.time() - t0), flush=True)
    allok = all(r['ok'] for r in res.values())
    print('RESULT %s: ALL OK = %s ; %d boxes ; %.1fs' % (args.cert, allok, sum(r.get('nboxes', 0) for r in res.values()), time.time() - t0), flush=True)
