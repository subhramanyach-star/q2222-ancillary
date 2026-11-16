"""Exact verification for the odd octagon vertices (k = 1,3,5,7), whose certificate is piecewise-algebraic:
   vertex 1: delta_0 = delta_1 = d(c),  d(c) = rho/(2 sqrt2) - rho/sqrt(2+2c)  for c >= c1 = 5 - 4 sqrt2,
                                         d(c) = -1/4                          for c <= c1.
We substitute u = sqrt(1 + sigma c) (sigma = +1 for k=1,5 and -1 for k=3,7) on the algebraic piece, multiply every
dual vector by u > 0, and reduce to polynomial nonnegativity in u over Q(sqrt2) (same exact method as for the even
vertices). The breakpoint is moved to the rational point sigma*c = -16/25 (u = 3/5 > 2 - sqrt2): the algebraic piece
is checked on u in [3/5, 14143/10000] (contains [3/5, sqrt2]), the constant piece on sigma*c in [-1, -16/25]."""
import sympy as sp
from exact_certificates import nonneg_on, vertex, RHO, r2, k
t = sp.symbols('t', real=True)
c_sym = sp.symbols('c', real=True)

def gnorm2(v, cc):
    return sp.expand(v[0]**2 + v[1]**2 + 2*cc*v[0]*v[1])

def verify_piece(a, b, C, delta, cc, var, lo, hi, mult, label):
    """check sqrt||B0+s0 c0||^2 + sqrt||B1+s1 c1||^2 <= R for all signs, on var in [lo,hi];
       delta, cc expressions in var; every vector multiplied by mult(var)>0 (also R)."""
    import exact_certificates as ec
    ok_all = True
    B0 = [b[0]/2 + delta[0], b[1]/2 + delta[1]]
    B1 = [b[0]/2 - delta[0], b[1]/2 - delta[1]]
    for s0 in (1, -1):
        for s1 in (1, -1):
            v0 = [sp.expand(sp.cancel(mult*(B0[i] + s0*C[0][i]))) for i in (0, 1)]
            v1 = [sp.expand(sp.cancel(mult*(B1[i] + s1*C[1][i]))) for i in (0, 1)]
            A = gnorm2(v0, cc); B = gnorm2(v1, cc)
            R = sp.expand(mult*sp.radsimp(1 - s0*a[0] - s1*a[1]))
            P1 = sp.expand(R**2 - A - B); P2 = sp.expand(P1**2 - 4*A*B)
            # nonneg_on works in the module's symbol c: substitute
            sub = {var: ec.c}
            okR, _ = nonneg_on(R.subs(sub), lo, hi)
            ok1, i1 = nonneg_on(P1.subs(sub), lo, hi); ok2, i2 = nonneg_on(P2.subs(sub), lo, hi)
            print("   %s (a,a')=(%+d,%+d): R>=0 %s | P1>=0 %s [%s] | P2>=0 %s [%s]" % (label, s0, s1, okR, ok1, i1, ok2, i2))
            ok_all &= okR and ok1 and ok2
    return ok_all

if __name__ == "__main__":
    res = {}
    lo_u, hi_u = sp.Rational(3, 5), sp.Rational(14143, 10000)
    for kk, sigma, sgn in [(1, 1, 1), (5, 1, -1), (7, -1, None), (3, -1, None)]:
        a, b, C = vertex(kk)
        print("=== octagon vertex %d: a=%s b=%s ===" % (kk, a, b))
        # algebraic piece: sigma*c = u^2 - 1, d(sigma c) = rho/(2 sqrt2) - rho/(sqrt2 u)
        dalg = RHO/(2*r2) - RHO/(r2*t)
        cc = sigma*(t**2 - 1)
        if kk == 1: delta = [dalg, dalg]
        if kk == 5: delta = [-dalg, -dalg]
        if kk == 7: delta = [-dalg, dalg]
        if kk == 3: delta = [dalg, -dalg]
        ok1 = verify_piece(a, b, C, delta, cc, t, lo_u, hi_u, t, "alg  u in [3/5,1.4143]")
        # constant piece: d = -1/4 on sigma*c in [-1,-16/25]; parametrise by c directly
        dcon = sp.Rational(-1, 4)
        if kk == 1: delta = [dcon, dcon]
        if kk == 5: delta = [-dcon, -dcon]
        if kk == 7: delta = [-dcon, dcon]
        if kk == 3: delta = [dcon, -dcon]
        if sigma == 1: lo, hi = -1, sp.Rational(-16, 25)
        else:          lo, hi = sp.Rational(16, 25), 1
        ok2 = verify_piece(a, b, C, delta, t, t, lo, hi, 1, "const c in [%s,%s]" % (lo, hi))
        res[kk] = ok1 and ok2
        print("   ==> vertex %d: Phi_c(delta(c)) <= 1 on all of [-1,1] (exact): %s" % (kk, res[kk]))
    print("SUMMARY odd vertices:", res)
