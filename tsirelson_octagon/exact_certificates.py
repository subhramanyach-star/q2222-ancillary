"""Exact (symbolic, Q(sqrt2)) verification of explicit global dual certificates delta(c) for the
Barizien-Bancal octagon vertices at the Tsirelson point:   Phi_c(delta(c)) <= 1 = g.P_T  for all c in [-1,1].

Phi_c(delta) = sum_x max_a [ a a_x + ||B_x + a c_x||_{G_c} ],  B_0 = b/2 + delta, B_1 = b/2 - delta.
Phi_c(delta) <= 1  <=>  for all (a,a'):  sqrt(A_a) + sqrt(B_a') <= R_{aa'} := 1 - a a_0 - a' a_1
                   <=>  P1 := R^2 - A - B >= 0  and  P2 := P1^2 - 4AB >= 0   (given R > 0, A, B >= 0).
Each P is a polynomial in c with coefficients in Q(sqrt2): P = U + sqrt2 V.  Its real zeros are zeros of
N = U^2 - 2V^2 in Q[c]; we isolate those exactly (sympy, rational isolating intervals), and evaluate sign(P)
exactly at a rational point of every gap and at the isolating-interval endpoints. This proves P >= 0 on [-1,1]
(by continuity: P has no sign change except at zeros of N), in exact arithmetic.
"""
import sympy as sp
c = sp.symbols('c', real=True)
r2 = sp.sqrt(2)
RHO = 1 - 1/r2
k = 1/(2*r2)

def gnorm2(v):
    return sp.expand(v[0]**2 + v[1]**2 + 2*c*v[0]*v[1])

def split_sqrt2(P):
    P = sp.expand(P)
    V = sp.expand(P.coeff(r2))
    U = sp.expand(P - r2*V)
    assert not U.has(r2) and not V.has(r2), (U, V)
    return sp.Poly(U, c, domain='QQ'), sp.Poly(V, c, domain='QQ')

def sign_exact(U, V, x):
    u = U.eval(x); v = V.eval(x)
    if u >= 0 and v >= 0: return 0 if (u == 0 and v == 0) else 1
    if u <= 0 and v <= 0: return -1
    if u > 0:  return 1 if u*u > 2*v*v else (-1 if u*u < 2*v*v else 0)
    return 1 if 2*v*v > u*u else (-1 if 2*v*v < u*u else 0)

def nonneg_on(P, lo=-1, hi=1):
    U, V = split_sqrt2(P)
    if U.is_zero and V.is_zero: return True, "identically 0"
    N = U**2 - 2*V**2
    if N.is_zero: raise ValueError("degenerate norm")
    ivs = sorted([(sp.Rational(a_), sp.Rational(b_)) for (a_, b_), m in N.intervals(eps=sp.Rational(1, 10**8))])
    lo, hi = sp.Rational(lo), sp.Rational(hi)
    edges = sorted(set([lo, hi] + [x for iv in ivs for x in iv if lo <= x <= hi]))
    pts = list(edges) + [(a_ + b_)/2 for a_, b_ in zip(edges[:-1], edges[1:])]
    sgns = [sign_exact(U, V, x) for x in pts]
    zeros_in = [iv for iv in ivs if iv[1] >= lo and iv[0] <= hi]
    return all(s >= 0 for s in sgns), "deg N=%d, %d real zeros of N in [lo,hi], %d exact sample points" % (N.degree(), len(zeros_in), len(pts))

def verify(a, b, C, delta, name):
    print("=== %s ===" % name)
    print("   a =", a, " b =", b)
    print("   delta(c) =", delta)
    B0 = [b[0]/2 + delta[0], b[1]/2 + delta[1]]
    B1 = [b[0]/2 - delta[0], b[1]/2 - delta[1]]
    ok_all = True
    for s0 in (1, -1):
        for s1 in (1, -1):
            A = gnorm2([B0[0] + s0*C[0][0], B0[1] + s0*C[0][1]])
            B = gnorm2([B1[0] + s1*C[1][0], B1[1] + s1*C[1][1]])
            R = sp.radsimp(1 - s0*a[0] - s1*a[1])
            assert float(R) > 0
            P1 = sp.expand(R**2 - A - B)
            P2 = sp.expand(P1**2 - 4*A*B)
            ok1, i1 = nonneg_on(P1); ok2, i2 = nonneg_on(P2)
            print("   (a,a')=(%+d,%+d): R=%s | P1>=0: %s [%s] | P2>=0: %s [%s]" % (s0, s1, R, ok1, i1, ok2, i2))
            ok_all &= ok1 and ok2
    print("   ==> Phi_c(delta(c)) <= 1 for all c in [-1,1] (exact):", ok_all)
    return ok_all

def vertex(kk):
    th = kk*sp.pi/4
    r0, r1 = sp.radsimp(RHO*sp.cos(th)), sp.radsimp(RHO*sp.sin(th))
    a = [sp.radsimp(sp.expand((r0 + r1)/r2)), sp.radsimp(sp.expand((r0 - r1)/r2))]
    b = [-r0, -r1]
    C = [[k, k], [k, -k]]
    return a, b, C

if __name__ == "__main__":
    res = {}
    kap = sp.Rational(1, 10)
    d0 = [RHO*c/2, -RHO/2 - kap*c**2]
    certs = {0: d0, 2: [d0[1], d0[0]], 4: [-d0[0], -d0[1]], 6: [-d0[1], -d0[0]]}
    for kk, dl in certs.items():
        a, b, C = vertex(kk)
        res[kk] = verify(a, b, C, [sp.expand(x) for x in dl], "octagon vertex %d (angle %d*pi/4)" % (kk, kk))
    print("SUMMARY even vertices:", res)
