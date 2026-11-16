"""Complements polygon_exact.py (which isolates the real roots of every numerator/denominator in the OPEN interval
(0, 2 - sqrt3) and finds none): checks that every mu_{a|x} and every non-tight local margin 1 - g.L is also nonzero at the
closed endpoints t = 1/40 and t = 2 - sqrt3 of Proposition 6.  The values are exact algebraic numbers (rational in t)
whose sign is read from a 40-digit evaluation; all are >= 2e-3 in absolute value.  Run: python3 check_polygon_endpoints.py"""
import sympy as sp, polygon_io
from sym_family import F
from exact_data import DET, lab
D = polygon_io.load('out/polygon_I.json'); B = D['B']; R = D['R']; t = B['x']
mn = None
for k, v in R.items():
    V = v['V']; ex = list(V['mu']) + [F(1 - (V['g'].T*L)[0]) for L in DET if lab(L) not in V['tight']]
    for e in ex:
        for pt in (sp.Rational(1, 40), 2 - sp.sqrt(3)):
            n = sp.N(sp.radsimp(e.subs(t, pt)), 40)
            assert n > 0, (k, e, pt, n); mn = n if mn is None or n < mn else mn
    print(k, V['tight'], ': %d conditions positive at both endpoints' % (2*len(ex)))
print('minimum value', sp.N(mn, 6))
