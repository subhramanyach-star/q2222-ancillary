"""sgen.py -- small-theta certificate for vertex V1 via ALICE-overlap slices (= Bob slices of the swapped functional
g' = P g_V1 at p' = P p_theta, target overlap c0' = 0).  Variable t = s^2 = sin^2(2 theta); all data are analytic in t
(even in s; C = sqrt(1-t), r = sqrt(1+t)), exact formulas from swap_data.json (swapdata.py, verified exactly).

Certificate: delta(c,t) = d0(t) + d1(t) c + t c^2 R(c), R polynomial with rational coefficients (cert json).
F_{s0 s1}(c,t) = 1 - s0 g0 - s1 g1 - ||u_{s0}||_c - ||w_{s1}||_c,  u = b/2 + delta + s0 e0,  w = b/2 - delta + s1 e1.
h_{K'_c}(g) <= 1  <=  F_{s0s1}(c,t) >= 0 for all four pairs (slice duality, as in the PRL proof).

Structure (asserted exactly below, from Taylor coefficients at t = 0 computed in exact arithmetic):
  u, w at t = 0 are the constant axis vectors U00, W00 (norm independent of c);  F(c,0) == 0 for all c and all pairs.
  Pairs (+-),(-+),(--):  G = F/t  = -(s0 g^0 + s1 g^1) - (2<u00,u^> + t||u^||^2)/(||u||+||u00||) - (same for w)
  Pair (++) (u00 = w00 = a e_0, a = 2/5):  G = F/t^2 =
     Lhh - ||u^||^2/(||u||+a) - ||w^||^2/(||w||+a) + (2<u00,u^>+t||u^||^2)<u00,u^>/(a(||u||+a)^2) + (same for w),
     Lhh = -(g^^0 + g^^1) - <e_0, g^^_b + g^^_e0 + g^^_e1>_c   (delta cancels; requires L(c,0) == 0, asserted).
  Here X^ = (X(t)-X(0))/t and X^^ = (X(t)-X(0)-t X'(0))/t^2, enclosed over a t-box [ta,tb] by the Taylor polynomial at
  0 (exact rational coefficients) plus the Lagrange remainder (Arb Taylor coefficient over [0,tb]).
Tests on a box [ca,cb] x [ta,tb] (G >= 0 on the box):
  T0: enclosure of G > 0.
  T2: hull H = conv(box, 0) small: G_cc > 0 on H x tbox.  Valid because G(0,t) = G_c(0,t) = 0 for all t (exact
      c0'-identities, swapdata.py, with mu' > 0).
  T3+/T3-: G_c < 0 on [ca,1] (resp. > 0 on [-1,cb]) and G(+-1,t) >= 0 (endpoint lemma, below).
Endpoint lemma: at c = +-1, ||v|| = |v0 +- v1|.  For each t-box and pair either G(e) encloses > 0, or the signs
  sigma_u, sigma_w of u0+-u1, w0+-w1 are certified equal, so F(e) = 1 - g'.L' with L' = (s0, s1, sigma, +-sigma)
  deterministic, and 1 - g'.L' = 1 - g.L >= 0 for all theta in (0, pi/4] by fullrange/comb.py (polygon combinatorics).
"""
import sys, json, argparse, itertools, sympy as sp, flint
from flint import arb, fmpq, arb_series
ap = argparse.ArgumentParser(); ap.add_argument('cert'); ap.add_argument('--prec', type=int, default=120)
ap.add_argument('--tstrips', type=int, nargs=2, default=None); ap.add_argument('--nc', type=int, default=64)
ap.add_argument('--maxbox', type=int, default=3000000); ap.add_argument('--K', type=int, default=10)
args = ap.parse_args(); flint.ctx.prec = args.prec; flint.ctx.cap = 64   # series length cap (default 10 would truncate!)
CERT = json.load(open(args.cert)); V = CERT['vertex']; KT = args.K
S_, C_, R_, T_ = sp.symbols('s C r t', positive=True)
SW = json.load(open(CERT.get('data', 'swap_data.json')))['vertices'][V]
def Q(x): x = str(x); return fmpq(*map(int, x.split('/'))) if '/' in x else fmpq(int(x))
def tosym(e):
    x = sp.sympify(e, locals={'s': S_, 'C': C_, 'r': R_}).subs(S_, sp.sqrt(T_))
    assert not x.has(S_); return x
EX = {k: [tosym(e) for e in SW[k]] for k in ('g', 'd0', 'd1', 'mu')}
# exact Taylor coefficients at t = 0
def taylor0(x, n):
    y = sp.series(x.subs({C_: sp.sqrt(1 - T_), R_: sp.sqrt(1 + T_)}), T_, 0, n).removeO()
    P = sp.Poly(sp.expand(y), T_); cs = [sp.Rational(P.coeff_monomial(T_ ** k)) for k in range(n)]
    return cs
TC = {k: [taylor0(x, KT) for x in EX[k]] for k in ('g', 'd0', 'd1', 'mu')}
def ev(x, env):
    if x.is_Symbol: return env[x.name]
    if x.is_Integer: return int(x)
    if x.is_Rational: return env['_q'](fmpq(int(x.p), int(x.q)))
    if x.is_Add:
        it = iter(x.args); a = ev(next(it), env)
        for y in it: a = a + ev(y, env)
        return a
    if x.is_Mul:
        it = iter(x.args); a = ev(next(it), env)
        for y in it: a = a * ev(y, env)
        return a
    if x.is_Pow:
        b, e = x.args; assert e.is_Integer; bb = ev(b, env); n = int(e)
        if n < 0: bb = 1 / bb; n = -n
        a = bb
        for _ in range(n - 1): a = a * bb
        return a
    raise ValueError(x)
def remainder(tb):
    """KT-th Taylor coefficient of each datum over the ball [0, tb] (Lagrange remainder)."""
    Tb = arb_series([arb((tb) / 2, (tb) / 2), 1], prec=KT + 1)
    env = {'t': Tb, 'C': (1 - Tb).sqrt(), 'r': (1 + Tb).sqrt(), '_q': lambda q: arb(q)}
    out = {}
    for k in TC:
        out[k] = []
        for x in EX[k]:
            sx = ev(x, env); cs = sx.coeffs()
            assert sx.prec > KT, 'series truncated'          # trailing zero coefficients may be dropped by coeffs()
            out[k].append(cs[KT] if len(cs) > KT else arb(0))
    return out
def hat(cs, Rk, t, m):
    """(X(t) - sum_{k<m} c_k t^k)/t^m for t in ball t, via Taylor coefficients cs (exact) + remainder Rk."""
    acc = arb(0); tp = arb(1)
    for k in range(m, KT): acc += arb(fmpq(int(cs[k].p), int(cs[k].q))) * tp; tp = tp * t
    return acc + Rk * tp
def tdata(ta, tb):
    t = arb((ta + tb) / 2, (tb - ta) / 2); Rm = remainder(tb)
    d = {'t': t}
    for k in ('g', 'd0', 'd1'):
        d[k + '0'] = [fmpq(int(c[0].p), int(c[0].q)) for c in TC[k]]
        d[k + 'h'] = [hat(TC[k][i], Rm[k][i], t, 1) for i in range(len(TC[k]))]
    d['ghh'] = [hat(TC['g'][i], Rm['g'][i], t, 2) for i in range(8)]
    d['mu'] = [hat(TC['mu'][i], Rm['mu'][i], t, 0) for i in range(4)]
    return d
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
# ---- exact structural assertions at t = 0
g0 = [c[0] for c in TC['g']]; g1 = [c[1] for c in TC['g']]; d00 = [c[0] for c in TC['d0']]
assert all(c[0] == 0 for c in TC['d1']), 'd1(0) != 0'
U00 = {s0: [g0[2 + i] / 2 + d00[i] + s0 * g0[4 + i] for i in range(2)] for s0 in (1, -1)}
W00 = {s1: [g0[2 + i] / 2 - d00[i] + s1 * g0[6 + i] for i in range(2)] for s1 in (1, -1)}
NRM = {}
for key, v in list(U00.items()) + [((9, k), v) for k, v in W00.items()]:
    assert (v[0] == 0) != (v[1] == 0), 'not a nonzero axis vector'   # ||v||_c = |v0| + |v1| for all c
for s0, s1 in COMBOS:
    T0 = 1 - s0 * g0[0] - s1 * g0[1]; nu = abs(U00[s0][0]) + abs(U00[s0][1]); nw = abs(W00[s1][0]) + abs(W00[s1][1])
    assert T0 == nu + nw, 'F(c,0) != 0'
assert U00[1] == W00[1] and U00[1][1] == 0 and U00[1][0] > 0
A_PP = U00[1][0]
# L(c,0) == 0 : -(g1_0+g1_1) - (g1_2+g1_4+g1_6) == 0 and g1_3+g1_5+g1_7 == 0
assert -(g1[0] + g1[1]) - (g1[2] + g1[4] + g1[6]) == 0 and g1[3] + g1[5] + g1[7] == 0, 'L(c,0) != 0'
Rco = [[Q(x) for x in CERT['R'][i]] for i in range(2)]      # R_i(c) = sum_k Rco[i][k] c^k
def Rpoly(c):
    out = []
    for i in range(2):
        acc = 0 * c
        for a in reversed(Rco[i]): acc = acc * c + arb(a)
        out.append(acc)
    return out
def A(x): x = sp.Rational(x); return arb(fmpq(int(x.p), int(x.q)))
def ip(a, b, c): return ((1 + c) * (a[0] + a[1]) * (b[0] + b[1]) + (1 - c) * (a[0] - a[1]) * (b[0] - b[1])) / 2
def Gall(d, c):
    """c: arb or arb_series in c.  Returns [G_++, G_+-, G_-+, G_--]."""
    t = d['t']; gh = d['gh']; Rc = Rpoly(c)
    dh = [d['d0h'][i] + d['d1h'][i] * c + c * c * Rc[i] for i in range(2)]
    Uh = {s0: [gh[2 + i] / 2 + dh[i] + s0 * gh[4 + i] for i in range(2)] for s0 in (1, -1)}
    Wh = {s1: [gh[2 + i] / 2 - dh[i] + s1 * gh[6 + i] for i in range(2)] for s1 in (1, -1)}
    out = []
    for s0, s1 in COMBOS:
        u00 = [A(x) for x in U00[s0]]; w00 = [A(x) for x in W00[s1]]; uh, wh = Uh[s0], Wh[s1]
        u = [u00[i] + t * uh[i] for i in range(2)]; w = [w00[i] + t * wh[i] for i in range(2)]
        nu = ip(u, u, c).sqrt(); nw = ip(w, w, c).sqrt()
        au = A(abs(U00[s0][0]) + abs(U00[s0][1])); aw = A(abs(W00[s1][0]) + abs(W00[s1][1]))
        qu = 2 * ip(u00, uh, c) + t * ip(uh, uh, c); qw = 2 * ip(w00, wh, c) + t * ip(wh, wh, c)
        if (s0, s1) != (1, 1):
            out.append(-(s0 * gh[0] + s1 * gh[1]) - qu / (nu + au) - qw / (nw + aw))
        else:
            a = au; hh = d['ghh']
            L = -(hh[0] + hh[1]) - ((hh[2] + hh[4] + hh[6]) + c * (hh[3] + hh[5] + hh[7]))
            out.append(L - ip(uh, uh, c) / (nu + a) - ip(wh, wh, c) / (nw + a)
                       + qu * ip(u00, uh, c) / (a * (nu + a) ** 2) + qw * ip(w00, wh, c) / (a * (nw + a) ** 2))
    return out
def coeffs(x, n):
    cs = x.coeffs(); return [cs[k] if k < len(cs) else arb(0) for k in range(n)]
def ctaylor(d, ca, cb, K=4, jmax=2):
    """Taylor-form enclosures over c in [ca,cb] of G^(j)/j!, j = 0..jmax, per pair."""
    cm = (ca + cb) / 2; h = arb((cb - ca) / 2) * arb(0, 1)  # c - cm in [-h, h]
    hball = arb(0, float((cb - ca) / 2) * (1 + 2 ** -40))
    hball = arb(0, arb((cb - ca) / 2).upper())
    P = Gall(d, arb_series([arb(cm), 1], prec=K + 1)); B = Gall(d, arb_series([arb((ca + cb) / 2, (cb - ca) / 2), 1], prec=K + 1))
    out = []
    from math import comb
    for i in range(4):
        pc = coeffs(P[i], K + 1); bc = coeffs(B[i], K + 1); res = []
        for j in range(jmax + 1):
            acc = arb(0); hp = arb(1)
            for k in range(j, K):
                acc += comb(k, j) * pc[k] * hp; hp = hp * hball
            acc += comb(K, j) * bc[K] * hp
            res.append(acc)
        out.append(res)
    return out
def endpoint_ok(d, e):
    """G(e, t) >= 0 for all four pairs over the t-box d (endpoint lemma)."""
    c = arb(e); Gv = Gall(d, c); t = d['t']; gh = d['gh']; Rc = Rpoly(c)
    dh = [d['d0h'][i] + d['d1h'][i] * c + c * c * Rc[i] for i in range(2)]
    res = []
    for k, (s0, s1) in enumerate(COMBOS):
        if Gv[k] > 0: res.append('pos'); continue
        u = [A(U00[s0][i]) + t * (gh[2 + i] / 2 + dh[i] + s0 * gh[4 + i]) for i in range(2)]
        w = [A(W00[s1][i]) + t * (gh[2 + i] / 2 - dh[i] + s1 * gh[6 + i]) for i in range(2)]
        su = u[0] + e * u[1]; sw = w[0] + e * w[1]
        if (su > 0 and sw > 0) or (su < 0 and sw < 0):
            sg = 1 if su > 0 else -1; L = (s0, s1, sg, e * sg)
            if L in DET_OK: res.append('det' + str(L)); continue
        return None
    return res
# deterministic points L' (swapped labelling) with 1 - g'.L' >= 0 proven by comb.py for all theta in (0, pi/4]
DET_OK = set(itertools.product((1, -1), repeat=4))   # comb.py: all 16 margins of the vertex are >= 0 (H1)
def accept(ca, cb, d, dcache):
    E = ctaylor(d, ca, cb); tests = []; Hc = {}
    for i in range(4):
        if E[i][0] > 0: tests.append('T0'); continue
        if cb > fmpq(1, 2):
            if 'h+' not in Hc: Hc['h+'] = ctaylor(d, ca, fmpq(1))
            if Hc['h+'][i][1] < 0 and dcache['ep+']: tests.append('T3+'); continue
        if ca < -fmpq(1, 2):
            if 'h-' not in Hc: Hc['h-'] = ctaylor(d, fmpq(-1), cb)
            if Hc['h-'][i][1] > 0 and dcache['ep-']: tests.append('T3-'); continue
        lo, hi = min(ca, fmpq(0)), max(cb, fmpq(0))
        if hi - lo <= fmpq(1, 2):
            if 'T2' not in Hc: Hc['T2'] = ctaylor(d, lo, hi)
            if Hc['T2'][i][2] > 0: tests.append('T2'); continue
        return ('split', i)
    return tests
def run_strip(ta0, tb0):
    stack = [(ta0, tb0)]; n = 0; nsplit = 0; area = fmpq(0); cnt = {}; nt = 0
    while stack:
        ta, tb = stack.pop(); d = tdata(ta, tb)
        if not all(m > 0 for m in d['mu']): return dict(ok=False, why='mu not certified positive', t=[str(ta), str(tb)])
        dc = {'ep+': endpoint_ok(d, 1) is not None, 'ep-': endpoint_ok(d, -1) is not None}
        cst = [(fmpq(-1) + fmpq(2 * i, args.nc), fmpq(-1) + fmpq(2 * (i + 1), args.nc)) for i in range(args.nc)]
        ok = True; loc = 0; acc = []
        while cst:
            ca, cb = cst.pop(); r = accept(ca, cb, d, dc)
            if isinstance(r, tuple):
                if cb - ca < fmpq(1, 2 ** 16): ok = False; break
                m = (ca + cb) / 2; cst += [(ca, m), (m, cb)]; nsplit += 1; continue
            loc += 1; acc.append(r)
        if not ok:
            if tb - ta < (tb0 - ta0) / 2 ** 16: return dict(ok=False, why='box', c=[str(ca), str(cb)], t=[str(ta), str(tb)], pair=r[1])
            m = (ta + tb) / 2; stack += [(ta, m), (m, tb)]; continue
        nt += 1; n += loc; area += tb - ta
        for r in acc:
            for x in r: cnt[x] = cnt.get(x, 0) + 1
        if nt % 50 == 0: print('   t-boxes done %d, area %.6g/%.6g, c-boxes %d' % (nt, float(area), float(tb0 - ta0), n), flush=True)
    return dict(ok=(area == tb0 - ta0), tboxes=nt, cboxes=n, nsplit=nsplit, tests=cnt)
if __name__ == '__main__':
    print('sgen %s vertex=%s prec=%d K=%d' % (args.cert, V, args.prec, KT), flush=True)
    js = range(len(CERT['tstrips'])) if args.tstrips is None else range(*args.tstrips)
    for j in js:
        a, b = map(Q, CERT['tstrips'][j]); res = run_strip(a, b)
        print('strip %d [%s, %s]: %s' % (j, a, b, json.dumps(res)), flush=True)
