"""Interval (Arb) evaluation of K-elements as functions of s, and exact sign decisions."""
import sympy as sp
from flint import arb, fmpq, ctx
from kfield import K, FS, S, s_sym
ctx.prec = 128

def _coeffs(poly_fs):
    # FS numerator/denominator are polynomial ring elements; get dense coeffs low->high
    P = sp.Poly(FS.ring.to_sympy(poly_fs) if hasattr(FS, 'ring') else poly_fs, s_sym)
    return [fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in reversed(P.all_coeffs())]

def _horner(cs, x):
    acc = arb(0)
    for c in reversed(cs): acc = acc * x + c
    return acc

class KF:
    """compiled evaluator of a K element"""
    def __init__(self, k: K):
        self.k = k
        self.parts = []
        for coef in k.c:
            e = FS.to_sympy(coef)
            n, d = sp.fraction(sp.together(e))
            self.parts.append((_coeffs(sp.Poly(n, s_sym).as_expr()), _coeffs(sp.Poly(d, s_sym).as_expr())))
    def __call__(self, s, C=None, r=None):
        if C is None: C = (1 - s*s).sqrt()
        if r is None: r = (1 + s*s).sqrt()
        vals = [_horner(n, s) / _horner(d, s) for n, d in self.parts]
        return vals[0] + vals[1]*C + vals[2]*r + vals[3]*C*r

def kderiv(x: K) -> K:
    a, b, c, d = x.c
    Cp = K(0, -S / (1 - S**2))          # C' = -s/C = -s C/(1-s^2)
    rp = K(0, 0, S / (1 + S**2))        # r' = s r/(1+s^2)
    Cv, rv = K(0, 1), K(0, 0, 1)
    da, db, dc, dd = (y.diff(FS.gens[0]) if hasattr(y, 'diff') else FS(0) for y in (a, b, c, d))
    return (K(da) + K(db)*Cv + K(b)*Cp + K(dc)*rv + K(c)*rp + K(dd)*Cv*rv + K(d)*(Cp*rv + Cv*rp))

def ball(lo, hi):
    """arb enclosing the real interval [lo,hi] (lo,hi fmpq or arb)."""
    lo = arb(lo); hi = arb(hi)
    lo_l = lo.lower(); hi_u = hi.upper()
    m = (lo_l + hi_u) / 2
    rad = (hi_u - lo_l) / 2
    return arb(m.mid(), (rad + m.rad()).upper())

def norm_poly(x: K):
    N = x.norm()
    e = sp.together(FS.to_sympy(N)); n, d = sp.fraction(e)
    return sp.Poly(n, s_sym), sp.Poly(d, s_sym)

def denoms_poly(x: K):
    out = 1
    for coef in x.c:
        out = out * sp.fraction(sp.together(FS.to_sympy(coef)))[1]
    return sp.Poly(out, s_sym)

def exact_sign_on(x: K, lo: sp.Rational, hi: sp.Rational):
    """Return (+1/-1, info) if x has constant strict sign on [lo,hi], decided exactly
    (no root of the norm / of denominators in [lo,hi], sign at midpoint by Arb);
    return (0, info) if undecided (norm has a root in [lo, hi]); 'zero' if x == 0."""
    if x.iszero(): return 'zero', 'identically zero'
    n, d = norm_poly(x)
    nr = n.count_roots(lo, hi); dr = denoms_poly(x).count_roots(lo, hi) + d.count_roots(lo, hi)
    if dr: return 0, 'denominator root in interval'
    if nr: return 0, 'norm has %d roots in interval' % nr
    mid = (lo + hi) / 2
    v = KF(x)(arb(fmpq(int(mid.p), int(mid.q))))
    if v > 0: return 1, 'norm root-free; value at mid %s' % v.str(5)
    if v < 0: return -1, 'norm root-free; value at mid %s' % v.str(5)
    return 0, 'midpoint sign undecided'

def arb_min_on(x: K, lo, hi, n=256):
    """rigorous lower bound of x on [lo,hi] by n subintervals with centered (mean-value) form"""
    f = KF(x); fp = KF(kderiv(x))
    lo = fmpq(int(lo.p), int(lo.q)); hi = fmpq(int(hi.p), int(hi.q))
    best_lo = None; best_hi = None
    for i in range(n):
        a = lo + (hi - lo) * i / n; b = lo + (hi - lo) * (i + 1) / n
        I = ball(a, b); m = arb((a + b) / 2)
        enc = f(m) + fp(I) * (I - m)
        nat = f(I)
        lb = max(enc.lower(), nat.lower()); ub = min(enc.upper(), nat.upper())
        best_lo = lb if best_lo is None else min(best_lo, lb)
        best_hi = ub if best_hi is None else min(best_hi, ub)
    return best_lo, best_hi
