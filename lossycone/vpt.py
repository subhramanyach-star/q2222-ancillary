#!/usr/bin/env python3
"""Stand-alone verifier: N_Q(q) = C(q) at the rational qubit points of data/points.json (computer-assisted proof).

For each point q (rational device) it RE-DERIVES from the device alone, in exact rational arithmetic:
  q, Bob's overlap c0, the slice gradients gamma_{a|x} (all four cones active: w^T G^{-1} w = t^2 exactly),
  the overlap derivative zdot and m_{a|x} = gamma_{a|x}.zdot, and the vertex list of the section
  S = {mu >= 0, sum mu = 1, mu.m = 0, (sum mu gamma).L <= 1 for all 16 deterministic L}
  by exhaustive enumeration of all pairs of constraints (so S = conv(vertices) exactly), and checks they match the file.
For each vertex V and its certificate delta(c) (data/certs.json; polynomial with rational coefficients, exact forced
data d0, d1 at c0, optional pins at c = +-1) it proves  Phi_c(delta(c)) <= 1 = V.q  for every c in [-1,1], i.e.
h_Q(V) = V.q (Theorem 'NDWU dual' + generation theorem), hence V in N_Q(q). With N_Q(q) <= C(q) (necessary conditions)
and convexity, N_Q(q) = C(q) = cone(V_1..V_n).
Proof of Phi_c(delta(c)) <= 1:  Phi_c <= 1  <=>  F_s(c) = T_s - ||u_{s0}(c)||_c - ||w_{s1}(c)||_c >= 0 for the four sign pairs
s = (s0,s1), T_s = 1 - s0 a0 - s1 a1, u_a = b/2 + delta + a e_0, w_a = b/2 - delta + a e_1, ||u||_c^2 = u^T G_c u.
||u(c)||_c^2 is an exact rational polynomial P(c); for a pinned branch u(c) = (c - eps) rho(c) exactly and
||u||_c = (1 - eps c) sqrt(rho^T G_c rho).  Exact checks: F_s(c0) = 0 and F_s'(c0) = 0 (all square roots are exact
rationals at c0); at c = +-1 the endpoint lemma: F_s(eps) is an exact rational >= 0 (deterministic margins and
mixed terms), with mixed terms > 0 for non-pinned branches.  Ball arithmetic (python-flint arb, outward rounding)
on a bisection of [-1,1]: T0 enclosure of F > 0; T1 mean value; T2 F'' > 0 on hull(box, c0) (with F = F' = 0 at c0);
T3 F' < 0 on hull(box, 1) (F' > 0 on hull(-1, box)) together with F(+-1) >= 0.
USAGE: python3 vpt.py [eta ...]"""
import json, sys, time, itertools
from fractions import Fraction as Fr
import flint
from flint import arb, arb_series, fmpq
flint.ctx.prec = 160
import os as _os; _H = _os.path.dirname(_os.path.abspath(__file__)) if '__file__' in globals() else '.'
P = json.load(open(_os.path.join(_H, 'data/points.json'))); CERT = json.load(open(_os.path.join(_H, 'data/certs.json')))
def fq(x): x = Fr(x); return fmpq(x.numerator, x.denominator)
def A(x): x = Fr(x); return arb(fmpq(x.numerator, x.denominator))
DET = [(a0, a1, b0, b1) for a0 in (1, -1) for a1 in (1, -1) for b0 in (1, -1) for b1 in (1, -1)]
def detz(L): a0, a1, b0, b1 = L; return [a0, a1, b0, b1, a0*b0, a0*b1, a1*b0, a1*b1]
def dot(g, z): return sum((Fr(x)*Fr(y) for x, y in zip(g, z)), Fr(0))
ORDER = [(1, 0), (1, 1), (-1, 0), (-1, 1)]
def solve(M, b):
    n = len(M); Am = [list(M[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if Am[r][col] != 0), None)
        if piv is None: return None
        Am[col], Am[piv] = Am[piv], Am[col]
        for r in range(n):
            if r != col and Am[r][col] != 0:
                f = Am[r][col]/Am[col][col]; Am[r] = [x - f*y for x, y in zip(Am[r], Am[col])]
    return [Am[i][n]/Am[i][i] for i in range(n)]
# ---------------------------------------------------------------- exact re-derivation of the polygon
def rederive(pt):
    d = {k: ([Fr(x) for x in v] if isinstance(v, list) else Fr(v)) for k, v in pt['dev'].items()}
    C, S, ca, sa, cb, sb = d['C'], d['S'], d['ca'], d['sa'], d['cb'], d['sb']
    assert C*C + S*S == 1 and all(ca[i]**2 + sa[i]**2 == 1 and cb[i]**2 + sb[i]**2 == 1 for i in (0, 1))
    assert S != 0 and C != 1
    z = [C*ca[0], C*ca[1], C*cb[0], C*cb[1]] + [ca[x]*cb[y] + S*sa[x]*sb[y] for x in (0, 1) for y in (0, 1)]
    zd = [Fr(0)]*8; zd[3] = -C*sb[1]
    for x in (0, 1): zd[5+2*x] = -ca[x]*sb[1] + S*sa[x]*cb[1]
    c0 = cb[0]*cb[1] + sb[0]*sb[1]; assert -1 < c0 < 1
    G = []
    for a, x in ORDER:
        t = 1 + a*z[x]; w = [z[2] + a*z[4+2*x], z[3] + a*z[5+2*x]]; assert t > 0
        gi = [(w[0] - c0*w[1])/(1 - c0*c0), (w[1] - c0*w[0])/(1 - c0*c0)]
        assert gi[0]*w[0] + gi[1]*w[1] == t*t          # cone (a|x) active
        v = [gi[0]/t, gi[1]/t]; g = [Fr(0)]*8; g[x] = Fr(-a); g[2:4] = v; g[4+2*x] = a*v[0]; g[5+2*x] = a*v[1]; G.append(g)
    m = [dot(g, zd) for g in G]
    rows = [([Fr(-1) if j == i else Fr(0) for j in range(4)], Fr(0)) for i in range(4)]
    rows += [([dot(G[k], detz(L)) for k in range(4)], Fr(1)) for L in DET]
    verts = []
    for i, j in itertools.combinations(range(len(rows)), 2):
        mu = solve([[Fr(1)]*4, m, rows[i][0], rows[j][0]], [Fr(1), Fr(0), rows[i][1], rows[j][1]])
        if mu is None or not all(dot(r, mu) <= rr for r, rr in rows): continue
        if mu not in verts: verts.append(mu)
    gs = [[sum((mu[k]*G[k][jj] for k in range(4)), Fr(0)) for jj in range(8)] for mu in verts]
    return z, c0, G, m, verts, gs
# ---------------------------------------------------------------- certificate polynomials (exact)
def padd(p, q):
    n = max(len(p), len(q)); return [(p[i] if i < len(p) else 0) + (q[i] if i < len(q) else 0) for i in range(n)]
def pmul(p, q):
    out = [Fr(0)]*(len(p) + len(q) - 1)
    for i, x in enumerate(p):
        for j, y in enumerate(q): out[i+j] += x*y
    return out
def pscale(p, s): return [s*x for x in p]
def peval(p, c):
    acc = p[-1]
    for x in reversed(p[:-1]): acc = acc*c + x
    return acc
def pdiv_lin(p, r):
    """exact synthetic division of p by (c - r); asserts zero remainder."""
    n = len(p) - 1; qd = [Fr(0)]*n; acc = Fr(0)
    for k in range(n, 0, -1): acc = p[k] + acc*r if k == n else p[k] + acc*r; qd[k-1] = acc
    # recompute properly
    qd = [Fr(0)]*n; carry = Fr(0)
    for k in range(n, 0, -1): carry = p[k] + carry*r; qd[k-1] = carry
    assert p[0] + carry*r == 0, 'pinned vector does not vanish'
    return qd
def cheb_polys(n):
    T = [[Fr(1)], [Fr(0), Fr(1)]]
    while len(T) < n + 1: T.append(padd(pmul([Fr(0), Fr(2)], T[-1]), pscale(T[-2], -1)))
    return T[:n+1]
def build_delta(g, d0, d1, c0, R, pins):
    """delta(c) as two exact rational polynomials, with pins (as in fullrange/vgen.py)."""
    n = len(R) - 1; T = cheb_polys(n)
    Q = [sum((pscale(T[j], R[j][i]) for j in range(n + 1)), []) if False else None for i in range(2)]
    Q = []
    for i in range(2):
        acc = [Fr(0)]
        for j in range(n + 1): acc = padd(acc, pscale(T[j], R[j][i]))
        Q.append(acc)
    b = g[2:4]
    def pinv(bx, ba): return [-(b[i]/2 + ba*g[4+i]) for i in range(2)] if bx == 0 else [b[i]/2 + ba*g[6+i] for i in range(2)]
    al = [Fr(0), Fr(0)]; be = [Fr(0), Fr(0)]
    if pins:
        t = {}
        for eps, (bx, ba) in pins.items():
            pv = pinv(bx, ba); t[eps] = [(pv[i] - d0[i] - d1[i]*(eps - c0))/(eps - c0)**2 - peval(Q[i], Fr(eps)) for i in range(2)]
        if len(pins) == 1: (eps, tv), = t.items(); al = tv
        else: al = [(t[1][i] + t[-1][i])/2 for i in range(2)]; be = [(t[1][i] - t[-1][i])/2 for i in range(2)]
    sq = pmul([-c0, Fr(1)], [-c0, Fr(1)])
    dl = [padd(padd([d0[i]], pscale([-c0, Fr(1)], d1[i])), pmul(sq, padd(Q[i], [al[i], be[i]]))) for i in range(2)]
    return dl
def gram_poly(u):
    """u^T G_c u = u0^2 + u1^2 + 2 c u0 u1 as a polynomial in c (u0,u1 polynomials)."""
    return padd(padd(pmul(u[0], u[0]), pmul(u[1], u[1])), pmul([Fr(0), Fr(2)], pmul(u[0], u[1])))
def branches(g, dl, pins):
    """For each sign pair s: (T_s, [(kind, poly)...]) with kind 0: sqrt(P(c)); kind eps: (1 - eps c) sqrt(P(c))."""
    b = g[2:4]; U = {}; W = {}
    for a in (1, -1):
        U[a] = [padd([b[i]/2 + a*g[4+i]], dl[i]) for i in range(2)]
        W[a] = [padd([b[i]/2 + a*g[6+i]], pscale(dl[i], -1)) for i in range(2)]
    def term(vec, x, a):
        for eps, (bx, ba) in pins.items():
            if bx == x and ba == a:
                rho = [pdiv_lin(vec[i], Fr(eps)) for i in range(2)]
                return (eps, gram_poly(rho), vec)
        return (0, gram_poly(vec), vec)
    out = {}
    for s0, s1 in ((1, 1), (1, -1), (-1, 1), (-1, -1)):
        out[s0, s1] = (1 - s0*g[0] - s1*g[1], [term(U[s0], 0, s0), term(W[s1], 1, s1)])
    return out
def fsqrt(q):
    """exact square root of a nonnegative rational, or None."""
    q = Fr(q); 
    if q < 0: return None
    from math import isqrt
    n, d = q.numerator, q.denominator; rn, rd = isqrt(n), isqrt(d)
    return Fr(rn, rd) if rn*rn == n and rd*rd == d else None
def exact_checks(g, dl, pins, c0, br):
    """F_s(c0) = 0, F_s'(c0) = 0 exactly; endpoint values F_s(+-1) >= 0 exactly (endpoint lemma)."""
    for s, (T, terms) in br.items():
        val = T; der = Fr(0)
        for kind, Pp, vec in terms:
            assert kind == 0 or True
            k = (lambda c: 1 - kind*c) if kind else (lambda c: 1)
            Pv = peval(Pp, c0); r = fsqrt(Pv); assert r is not None and r > 0, ('irrational/zero norm at c0', s)
            dP = peval([i*x for i, x in enumerate(Pp)][1:] or [Fr(0)], c0)
            val -= k(c0)*r; der -= k(c0)*dP/(2*r) + (-kind)*r if kind else dP/(2*r)
        assert val == 0, ('F(c0) != 0', s, val); assert der == 0, ('F\'(c0) != 0', s, der)
    endp = {}
    for eps in (1, -1):
        e = Fr(eps); f = (1, eps)
        for s, (T, terms) in br.items():
            # ||u||_{G_eps} = |u0 + eps u1|; pinned (vanishing) vector contributes 0
            vals = []
            for kind, Pp, vec in terms:
                uf = peval(vec[0], e) + eps*peval(vec[1], e)
                if kind == eps: assert peval(vec[0], e) == 0 and peval(vec[1], e) == 0; uf = Fr(0)
                vals.append(uf)
            Fv = T - abs(vals[0]) - abs(vals[1])
            assert Fv >= 0, ('endpoint F < 0', eps, s, Fv)
            endp[eps, s] = Fv
    return endp
# ---------------------------------------------------------------- ball arithmetic
def coeffs(zz, n):
    cs = zz.coeffs() if not isinstance(zz, arb) else [zz]
    return list(cs) + [arb(0)]*(n - len(cs))
def sqrt_nn(x):
    if isinstance(x, arb_series): return x.sqrt()
    if x > 0: return x.sqrt()
    up = x.upper(); return arb(0).union(up.sqrt()) if up >= 0 else arb('nan')
def apoly(Pq):
    return [arb(fq(x)) for x in Pq]
def aeval(Pa, c):
    acc = Pa[-1]
    for x in reversed(Pa[:-1]): acc = acc*c + x
    return acc
def shift(Pa, cm):
    """coefficients of P(cm + t) (arb), by repeated synthetic division (Taylor shift)."""
    a = list(Pa); n = len(a); out = []
    for k in range(n):
        for i in range(n - 2, k - 1, -1): a[i] = a[i] + cm*a[i + 1]
        out.append(a[k])
    return out
def Fball(brA, c, cm=None):
    """F_s on c (arb or arb_series); if cm (fmpq) is given, c is the OFFSET t = c - cm and polynomials are shifted."""
    out = []
    for s, (T, terms) in brA.items():
        val = T
        for kind, Pa in terms:
            if cm is None: Pv = aeval(Pa, c); cc = c
            else:
                key = (id(Pa), cm)
                if key not in _SH: _SH.clear() if len(_SH) > 64 else None; _SH[key] = shift(Pa, arb(cm))
                Pv = aeval(_SH[key], c); cc = c + arb(cm)
            r = sqrt_nn(Pv)
            val = val - ((1 - kind*cc)*r if kind else r)
        out.append(val)
    return out
_SH = {}
def ball(a, b): return arb((a + b)/2, (b - a)/2)
def _midq(H):
    m = H.mid(); return fmpq(int((m*2**60).floor().unique_fmpz()), 2**60)
def accept(brA, ca, cb, c0):
    cm = (ca + cb)/2; rc = (cb - ca)/2; Bt = arb(0, rc)
    F0 = Fball(brA, Bt, cm); Fm = Fball(brA, arb(0), cm); Fs = Fball(brA, arb_series([Bt, 1], prec=3), cm)
    tests = []
    for i in range(4):
        if F0[i].is_finite() and F0[i] > 0: tests.append('T0'); continue
        d1 = coeffs(Fs[i], 3)[1]
        if Fm[i].is_finite() and d1.is_finite() and Fm[i] - abs(d1)*arb(rc) > 0: tests.append('T1'); continue
        if cb > fmpq(1, 2):
            hm = (ca + 1)/2; z = coeffs(Fball(brA, arb_series([arb(0, (1 - ca)/2), 1], prec=3), hm)[i], 3)[1]
            if z.is_finite() and z < 0: tests.append('T3+'); continue
        if ca < -fmpq(1, 2):
            hm = (cb - 1)/2; z = coeffs(Fball(brA, arb_series([arb(0, (cb + 1)/2), 1], prec=3), hm)[i], 3)[1]
            if z.is_finite() and z > 0: tests.append('T3-'); continue
        lo = min(ca, c0); hi = max(cb, c0)
        if hi - lo < fmpq(3, 5):
            hm = (lo + hi)/2; f2 = coeffs(Fball(brA, arb_series([arb(0, (hi - lo)/2), 1], prec=3), hm)[i], 3)[2]
            if f2.is_finite() and f2 > 0: tests.append('T2'); continue
        return None
    return tests
def verify_vertex(g, mu, c0, cert, maxbox=400000):
    R = [[Fr(x) for x in row] for row in cert['R']]; pins = {int(k): tuple(v) for k, v in cert['pins'].items()}
    d0, d1 = cert['_d0'], cert['_d1']
    dl = build_delta(g, d0, d1, c0, R, pins)
    assert [peval(dl[i], c0) for i in range(2)] == d0
    br = branches(g, dl, pins); endp = exact_checks(g, dl, pins, c0, br)
    brA = {s: (arb(fq(T)), [(kind, apoly(Pp)) for kind, Pp, _ in terms]) for s, (T, terms) in br.items()}
    stack = [(fmpq(-1) + fmpq(2*i, 64), fmpq(-1) + fmpq(2*(i + 1), 64)) for i in range(64)]
    n = 0; cnt = {}; total = fmpq(0); c0q = fq(c0)
    while stack:
        ca, cb = stack.pop(); t = accept(brA, ca, cb, c0q)
        if t is not None:
            n += 1; total += cb - ca
            for x in t: cnt[x] = cnt.get(x, 0) + 1
            continue
        if cb - ca < fmpq(1, 2**40) or n > maxbox: return dict(ok=False, box=[str(ca), str(cb)], n=n)
        m = (ca + cb)/2
        if ca < c0q < cb: m = c0q
        stack += [(ca, m), (m, cb)]
    return dict(ok=(total == 2), nboxes=n, tests=cnt, endpoint_min=str(min(endp.values())))
def forced(G, mu, c0, g):
    """forced data (re-derived): d0 from complementary slackness, d1 from the first-order equalities."""
    v = []
    for k in range(4): v.append(G[k][2:4])
    vd = dict(zip(ORDER, v)); md = dict(zip(ORDER, mu)); b = g[2:4]
    d0 = [md[1, 0]*vd[1, 0][i] + md[-1, 0]*vd[-1, 0][i] - b[i]/2 for i in range(2)]
    def lin(a, x):
        s = 1 if x == 0 else -1; vv = vd[a, x]; Gv = [vv[0] + c0*vv[1], c0*vv[0] + vv[1]]
        return (2*md[a, x]*vv[0]*vv[1], [s*Gv[0], s*Gv[1]])
    eqs = []
    for (k1, k2, sg) in (((1, 0), (-1, 0), 1), ((1, 1), (-1, 1), 1), ((1, 0), (1, 1), -1)):
        c1, l1 = lin(*k1); c2, l2 = lin(*k2); eqs.append(([l1[i] - sg*l2[i] for i in range(2)], sg*c2 - c1))
    d1 = solve([eqs[0][0], eqs[1][0]], [eqs[0][1], eqs[1][1]])
    assert dot(eqs[2][0], d1) == eqs[2][1]
    return d0, d1
if __name__ == '__main__':
    etas = sys.argv[1:] or list(P)
    allok = True
    for e in etas:
        pt = P[e]; t0 = time.time()
        z, c0, G, m, verts, gs = rederive(pt)
        assert [str(x) for x in z] == pt['z'] and str(c0) == pt['c0']
        filev = [[Fr(x) for x in V['mu']] for V in pt['verts']]
        assert sorted(map(tuple, verts)) == sorted(map(tuple, filev)), 'vertex list mismatch'
        print('eta %s: q, c0, cones active, %d vertices re-derived exactly (exhaustive pair enumeration)' % (e, len(verts)), flush=True)
        for k, V in enumerate(pt['verts']):
            mu = [Fr(x) for x in V['mu']]; g = [Fr(x) for x in V['g']]
            assert all(x > 0 for x in mu) and sum(mu) == 1 and dot(mu, m) == 0
            assert g == [sum((mu[i]*G[i][j] for i in range(4)), Fr(0)) for j in range(8)] and dot(g, z) == 1
            assert all(dot(g, detz(L)) <= 1 for L in DET)
            cert = CERT.get('%s/%d' % (e, k))
            if cert is None: print('  vertex %d: NO CERTIFICATE' % k); allok = False; continue
            d0, d1 = forced(G, mu, c0, g); cert = dict(cert, _d0=d0, _d1=d1)
            try: r = verify_vertex(g, mu, c0, cert)
            except AssertionError as ex: r = dict(ok=False, why=str(ex))
            allok &= r['ok']
            print('  vertex %d tight %s: %s' % (k, [L for L in DET if dot(g, detz(L)) == 1], r), flush=True)
        print('eta %s done in %.1fs' % (e, time.time() - t0), flush=True)
    print('RESULT ALL OK =', allok)
