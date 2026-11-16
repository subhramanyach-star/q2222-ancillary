#!/usr/bin/env python3
"""Short, self-contained exact check of the even-vertex certificate of Theorem 3(a) (Tsirelson point P_T).

Vertex g_0 = BB functional rho((A0+A1)/sqrt2 - B0) + CHSH/(2 sqrt2), rho = 1 - 1/sqrt2:
   a = (rho/sqrt2, rho/sqrt2), b = (-rho, 0), e_0 = k(1,1), e_1 = k(1,-1), k = 1/(2 sqrt2).
Certificate delta_0(c) = (rho c/2, -rho/2 - c^2/10).  (Vertices g_2, g_4, g_6 follow by the symmetries in the paper and
are checked as well.)  Claim: Phi_c(delta_0(c)) <= 1 for all c in [-1,1], with equality exactly at c in {-1,0,1}.

Phi_c <= 1  <=>  for all a,a' = +-1:  sqrt(A_a) + sqrt(B_a') <= T_aa' := 1 - a a0 - a' a1   (T > 0)
            <=>  P1 := T^2 - A - B >= 0  and  P2 := P1^2 - 4AB >= 0  on [-1,1],
with A_a = ||B_0 + a e_0||^2_{G_c}, B_a' = ||B_1 + a' e_1||^2_{G_c}, B_{0,1} = b/2 +- delta.  P1, P2 lie in Q(sqrt2)[c].
Writing P = U + sqrt2 V (U, V in Q[c]), every real zero of P is a zero of N = U^2 - 2V^2 in Q[c].  We isolate the real
zeros of N in disjoint rational intervals (sympy, exact), and decide sign(P) EXACTLY at all interval endpoints and one
rational point per gap by comparing U^2 with 2V^2 (sign of U + sqrt2 V for rational U, V).  Since P is continuous and
changes sign only at zeros of N, and the isolating intervals contain exactly one zero of N each, nonnegativity of P at
all these points proves P >= 0 on [-1,1].  The equality set is read off from the exact factorisation of P2 over
Q(sqrt2) (step5_equality_sets.py).

Output: a PASS/FAIL line per condition and the file thm2a_even_polys.txt with the exact coefficients of U and V for all
32 polynomials, for independent auditing.   Run: python3 verify_thm2a_even.py   (< 2 s; sympy only)."""
import sympy as sp
from exact_certificates import c, r2, RHO, k, gnorm2, split_sqrt2, nonneg_on, vertex
kap = sp.Rational(1, 10)
d0 = [RHO * c / 2, -RHO / 2 - kap * c ** 2]
certs = {0: d0, 2: [d0[1], d0[0]], 4: [-d0[0], -d0[1]], 6: [-d0[1], -d0[0]]}
allok = True
with open('thm2a_even_polys.txt', 'w') as f:
    f.write('# vertex kk, sign pair (a,a\'), polynomial P1 or P2 = U + sqrt2*V; coefficients of U and V in c, highest degree first\n')
    for kk, dl in certs.items():
        a, b, Cm = vertex(kk); dl = [sp.expand(x) for x in dl]
        B0 = [b[0] / 2 + dl[0], b[1] / 2 + dl[1]]; B1 = [b[0] / 2 - dl[0], b[1] / 2 - dl[1]]
        for s0 in (1, -1):
            for s1 in (1, -1):
                A = gnorm2([B0[0] + s0 * Cm[0][0], B0[1] + s0 * Cm[0][1]])
                B = gnorm2([B1[0] + s1 * Cm[1][0], B1[1] + s1 * Cm[1][1]])
                T = sp.radsimp(1 - s0 * a[0] - s1 * a[1]); assert float(T) > 0
                P1 = sp.expand(T ** 2 - A - B); P2 = sp.expand(P1 ** 2 - 4 * A * B)
                for nm, P in (('P1', P1), ('P2', P2)):
                    U, Vv = split_sqrt2(P); ok, info = nonneg_on(P)
                    allok &= ok
                    f.write('g_%d (%+d,%+d) %s: U = %s ; V = %s\n' % (kk, s0, s1, nm, U.all_coeffs(), Vv.all_coeffs()))
                    print('g_%d (a,a\')=(%+d,%+d) %s >= 0 on [-1,1]: %s   [%s]' % (kk, s0, s1, nm, ok, info))
print('ALL 32 CONDITIONS PASS:', allok)
