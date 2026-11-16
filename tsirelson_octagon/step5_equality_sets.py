"""Step 5 (exact part 1): equality set E_k = {c in [-1,1] : Phi_c(delta_k(c)) = 1} of the even-vertex certificates.
Phi_c(delta) = 1  <=>  for some sign pair (a,a'): P2_{aa'}(c) = 0 (P1 >= 0 is already proved), so
E_k = union over (a,a') of the real zeros in [-1,1] of P2_{aa'}.  We factor each P2 into irreducibles over Q(sqrt2)
(sympy, extension=sqrt2).  Irreducible factors are squarefree, so all their real zeros are simple; a factor F has a zero
in an isolating interval I of N_F = F*conj(F) (in Q[c]) iff F changes sign on I (or vanishes exactly at an endpoint).
Everything is exact rational / Q(sqrt2) arithmetic."""
import sympy as sp
from exact_certificates import c, r2, RHO, k, gnorm2, split_sqrt2, sign_exact, vertex

def conj(P):
    U, V = split_sqrt2(P)
    return U, V

def zeros_in(F, lo=-1, hi=1):
    """exact list of real zeros of F (squarefree, coeffs in Q(sqrt2)) in [lo,hi], as isolating rational intervals."""
    U, V = split_sqrt2(F)
    if V.is_zero:
        N = U
    else:
        N = U**2 - 2*V**2
    out = []
    lo, hi = sp.Rational(lo), sp.Rational(hi)
    for (a_, b_), m in N.intervals(eps=sp.Rational(1, 10**12)):
        a_, b_ = sp.Rational(a_), sp.Rational(b_)
        if b_ < lo or a_ > hi: continue
        if a_ == b_:
            if sign_exact(U, V, a_) == 0: out.append((a_, a_))
            continue
        sa, sb = sign_exact(U, V, a_), sign_exact(U, V, b_)
        if sa == 0: out.append((a_, a_)); continue
        if sb == 0: out.append((b_, b_)); continue
        if sa != sb:
            if a_ >= lo and b_ <= hi: out.append((a_, b_))
            else: raise ValueError("zero isolating interval straddles the domain boundary; refine")
    return sorted(set(out))

def equality_set(a, b, C, delta, label):
    B0 = [b[0]/2 + delta[0], b[1]/2 + delta[1]]
    B1 = [b[0]/2 - delta[0], b[1]/2 - delta[1]]
    allz = set()
    for s0 in (1, -1):
        for s1 in (1, -1):
            A = gnorm2([B0[0] + s0*C[0][0], B0[1] + s0*C[0][1]])
            B = gnorm2([B1[0] + s1*C[1][0], B1[1] + s1*C[1][1]])
            R = sp.radsimp(1 - s0*a[0] - s1*a[1])
            P1 = sp.expand(R**2 - A - B); P2 = sp.expand(P1**2 - 4*A*B)
            fac = sp.factor_list(P2, c, extension=r2)
            zs = []
            for F, mult in fac[1]:
                for iv in zeros_in(sp.expand(F)):
                    zs.append((iv, mult, sp.degree(F, c)))
                    allz.add(iv)
            print("   %s (a,a')=(%+d,%+d): P2 = const * prod of %d irreducible factors over Q(sqrt2) (degrees %s);"
                  " zeros in [-1,1]: %s" % (label, s0, s1, len(fac[1]), [sp.degree(F, c) for F, m in fac[1]],
                  [("c=%s" % iv[0] if iv[0] == iv[1] else "c in (%s,%s)" % (float(iv[0]), float(iv[1])), "mult %d" % m) for iv, m, d in zs]))
    return allz

if __name__ == "__main__":
    kap = sp.Rational(1, 10)
    d0 = [RHO*c/2, -RHO/2 - kap*c**2]
    certs = {0: d0, 2: [d0[1], d0[0]], 4: [-d0[0], -d0[1]], 6: [-d0[1], -d0[0]]}
    res = {}
    for kk, dl in certs.items():
        a, b, C = vertex(kk)
        print("=== vertex %d ===" % kk)
        z = equality_set(a, b, C, [sp.expand(x) for x in dl], "v%d" % kk)
        res[kk] = sorted(z)
        print("   ==> E_%d =" % kk, [str(iv[0]) if iv[0] == iv[1] else iv for iv in res[kk]])
    ok = all(all(iv[0] == iv[1] and iv[0] in (-1, 0, 1) for iv in v) and len(v) == 3 for v in res.values())
    print("SUMMARY: every even-vertex equality set is exactly {-1,0,1}:", ok)
