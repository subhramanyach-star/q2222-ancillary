"""Exact polygon combinatorics of the candidate section C_theta for theta in (0, pi/4] (s = sin 2theta in (0, 1]).
Each condition (mu_{a|x} > 0 at a vertex; 1 - g.L > 0 for non-facet L; 1 - g.L == 0 for facet L) is an element x of
K = Q(s)(C, r).  Its sign on a closed s-interval is decided exactly: the real roots of the norm N(x) (and of all
denominators) in [0, 1] are isolated with rational intervals; x is evaluated in Arb on every isolating interval that
lies inside the range (if the enclosure excludes 0, x has no zero there) and at a rational point of every gap."""
import json, sys, itertools, sympy as sp
from flint import arb, fmpq, ctx
from kfield import K, from_sympy, s_sym, FS
from kfun import KF, norm_poly, denoms_poly
ctx.prec = 200
S_, C_, R_ = sp.symbols('s C r', positive=True)
def load(path):
    J = json.load(open(path)); loc = {'s': S_, 'C': C_, 'r': R_}
    return {n: {k: [sp.sympify(e, locals=loc) for e in V[k]] for k in ('g', 'mu')} for n, V in J['vertices'].items()}
D = load('thm3_data.json'); D.update(load('oct_data.json'))
def Lv(L): A0, A1, B0, B1 = L; return [A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1]
def lab(L): return ''.join('+' if v > 0 else '-' for v in L)
FACETS = {'V1': ('++++', '+--+'), 'V2': ('+++-', '+--+'), 'V1s': ('+++-', '+-++'), 'V3': ('+-++', '-++-'),
          'V4': ('-++-', '---+'), 'V3s': ('++++', '---+'), 'Wa': ('++++', '-+--'), 'Wb': ('+--+', '-+--'),
          'Wa_s': ('+++-', '----'), 'Wb_s': ('+-++', '----')}
DET = list(itertools.product((1, -1), repeat=4))
def sval(s):  # rational s -> arb triple
    sv = arb(s); return sv
def sign_on(x, lo, hi, allow_zero_at=()):
    """sign of x on [lo,hi] (sympy Rationals); zeros of x allowed only at points listed in allow_zero_at (exact
    endpoints given as ('rat', q) or ('root', poly, (a,b)) ).  Returns (sign, info)."""
    if x.iszero(): return 'zero', 'identically zero'
    n, d = norm_poly(x); den = denoms_poly(x) * d
    roots = []
    for P in (n, den):
        P = sp.Poly(P, s_sym)
        if P.degree() <= 0: continue
        for iv, mult in P.intervals(inf=lo, sup=hi, eps=sp.Rational(1, 10**12)):
            roots.append(iv)
    f = KF(x); bad = []
    for (a, b) in roots:
        # exclude roots at an allowed endpoint
        I = arb((fmpq(int(a.p), int(a.q)) + fmpq(int(b.p), int(b.q))) / 2, float(b - a) / 2 + 1e-40)
        v = f(I)
        if v > 0 or v < 0: continue
        bad.append((a, b))
    pts = sorted(set([lo, hi] + [a for a, b in roots] + [b for a, b in roots]))
    signs = set()
    for a, b in zip(pts[:-1], pts[1:]):
        m = (a + b) / 2
        if any(aa <= m <= bb for aa, bb in bad): continue
        v = f(arb(fmpq(int(m.p), int(m.q))))
        signs.add(1 if v > 0 else -1 if v < 0 else 0)
    return signs, bad

def rat2fmpq(q): q = sp.Rational(q); return fmpq(int(q.p), int(q.q))
def sign_open(x, lo, hi, endpoly=None):
    """Exact: constant strict sign of x on the OPEN interval (lo, hi) (lo, hi rational, 0 <= lo < hi <= 1).
    endpoly: minimal polynomial of an algebraic endpoint approximated by lo or hi (its factor is divided out of
    the norm, so that roots AT that endpoint are not mistaken for interior roots).
    Returns +1/-1, or raises if x has (or may have) a zero / pole inside."""
    if x.iszero(): return 0
    n, d = norm_poly(x); polys = [sp.Poly(n, s_sym), sp.Poly(denoms_poly(x) * d, s_sym)]
    if endpoly is not None:
        for k in range(2):
            while polys[k].degree() > 0 and sp.rem(polys[k], endpoly).is_zero: polys[k] = sp.quo(polys[k], endpoly)
            if polys[k].degree() > 0 and endpoly.gcd(polys[k]).degree() > 0: raise ValueError('partial common factor with endpoint poly')
    cuts = []
    for k, P in enumerate(polys):
        if P.degree() <= 0: continue
        for (a, b), mult in P.intervals(eps=sp.Rational(1, 10**15)):
            if b <= lo or a >= hi: continue
            if a < lo or b > hi:
                # isolating interval straddles an endpoint: shrink until it does not, unless the root IS the endpoint
                for e in range(16, 60, 4):
                    ivs = [iv for iv, m in P.intervals(eps=sp.Rational(1, 10**e)) if not (iv[1] <= lo or iv[0] >= hi)]
                    if all(lo <= iv[0] and iv[1] <= hi for iv in ivs): break
                    if any(P.eval(e_) == 0 for e_ in (lo, hi)): break
                ivs = [iv for iv, m in P.intervals(eps=sp.Rational(1, 10**e)) if not (iv[1] <= lo or iv[0] >= hi)]
                for (a2, b2) in ivs:
                    if P.eval(lo) == 0 and a2 <= lo <= b2 and a2 == b2: continue
                    if P.eval(hi) == 0 and a2 <= hi <= b2 and a2 == b2: continue
                    if a2 == b2 and (a2 == lo or a2 == hi): continue
                    cuts.append((max(a2, lo), min(b2, hi), k))
                continue
            if a == b and (a == lo or a == hi): continue
            cuts.append((a, b, k))
    f = KF(x); signs = set()
    for a, b, k in cuts:
        if k == 1: raise ValueError('pole of x (denominator root) in (%s, %s)' % (a, b))
        ok = False
        for e in range(0, 8):
            I = arb((rat2fmpq(a) + rat2fmpq(b)) / 2, float(b - a) / 2 * 1.000001 + 1e-60)
            v = f(I)
            if v > 0 or v < 0: ok = True; signs.add(1 if v > 0 else -1); break
            raise ValueError('possible zero of x in [%s, %s] (enclosure %s)' % (float(a), float(b), v.str(5)))
    pts = sorted(set([lo, hi] + [a for a, b, k in cuts] + [b for a, b, k in cuts]))
    for a, b in zip(pts[:-1], pts[1:]):
        m = (a + b) / 2
        if any(aa <= m <= bb for aa, bb, k in cuts): continue
        v = f(arb(rat2fmpq(m)))
        if not (v > 0 or v < 0): raise ValueError('undecided sign at %s' % m)
        signs.add(1 if v > 0 else -1)
    if len(signs) != 1: raise ValueError('sign change: %s' % signs)
    return signs.pop()
def at_point(x, sball):
    """Arb sign of x at s in sball (C, r from s)."""
    v = KF(x)(sball)
    return 1 if v > 0 else -1 if v < 0 else 0, v
def conditions(name):
    V = D[name]; g = [from_sympy(e) for e in V['g']]; mu = [from_sympy(e) for e in V['mu']]
    out = [('mu%d' % i, mu[i]) for i in range(4)]
    for L in DET:
        x = K(1) - sum((gi * K(FS(li)) for gi, li in zip(g, Lv(L))), K(0))
        out.append(('L' + lab(L), x))
    return out
if __name__ == '__main__':
    import time
    t0 = time.time()
    # s* : root of s^6 - 8 s^4 + 24 s^2 - 16 in (0.95, 0.96)
    Pst = sp.Poly(s_sym**6 - 8*s_sym**4 + 24*s_sym**2 - 16, s_sym)
    ivs = [iv for iv, m in Pst.intervals(eps=sp.Rational(1, 10**40)) if 0 < iv[0] < 1]
    assert len(ivs) == 1; sa, sb = ivs[0]
    sstar = arb((rat2fmpq(sa) + rat2fmpq(sb)) / 2, float(sb - sa))
    print('s* in [%s, %s] ~ %.15f ; theta*/pi = %s' % (float(sa), float(sb), float(sa), (sstar.asin() / 2 / arb.pi()).str(15)))
    REG = {'hex': (['V1', 'V2', 'V1s', 'V3', 'V4', 'V3s'], sp.Integer(0), sa), 'oct': (['Wa', 'Wb', 'V2', 'Wa_s', 'Wb_s', 'V3', 'V4', 'V3s'], sb, sp.Integer(1))}
    allok = True
    for reg, (names, lo, hi) in REG.items():
        print('== regime %s: open s-interval (%s, %s); exact facet identities; endpoints by Arb' % (reg, float(lo), float(hi)), flush=True)
        for name in names:
            res = {}
            for cl, x in conditions(name):
                facet = cl[0] == 'L' and cl[1:] in FACETS[name]
                if facet:
                    res[cl] = 'facet: exact zero' if x.iszero() else 'FACET NOT ZERO'; allok &= x.iszero(); continue
                try:
                    sg = sign_open(x, lo, hi, Pst)
                except ValueError as e:
                    res[cl] = 'FAIL ' + str(e); allok = False; continue
                res[cl] = sg
                if sg != 1: allok = False
            # endpoint values
            ends = {}
            for ename, sb_ in (('s=0', arb(0)), ('s*', sstar), ('s=1', arb(1))):
                if (reg == 'hex' and ename == 's=1') or (reg == 'oct' and ename == 's=0'): continue
                vals = {cl: at_point(x, sb_)[0] for cl, x in conditions(name) if not (cl[0] == 'L' and cl[1:] in FACETS[name])}
                ends[ename] = {k: v for k, v in vals.items() if v != 1}
            bad = {k: v for k, v in res.items() if v != 1 and not str(v).startswith('facet')}
            print('  %-5s facets %s: open-interval conditions all > 0: %s %s ; endpoint non-positive: %s' % (name, FACETS[name], not bad, bad if bad else '', ends), flush=True)
    print('ALL OK (open intervals):', allok, '%.1fs' % (time.time() - t0))
