# Exact check (added for the v13 bundle) that the V2, V4 vertex data in thm3_data.json equal the closed forms printed in
# Data used by v24_conditions.py for the hexagonal section:
#   V2 = ((r+C)/(2D), 0, bbar, bbar, kappa, kappa, nu, -nu), D = r+1, bbar = (C-1)/(2D), kappa = (C+r+2s^2)/(4D), nu = (r+2-C)/(4D);
#   V4: same form with kappa, nu from C -> -C, and a0 = (C-r)/(2D), bbar = (C+1)/(2D).
# Equality is decided exactly in K = Q(s)(C, r), C^2 = 1-s^2, r^2 = 1+s^2 (kfield.from_sympy; zero test is exact).
import json, sympy as sp
from comb import S_, C_, R_
from kfield import from_sympy
s, C, r = S_, C_, R_; D = r + 1
def form(a0, bbar, kap, nu): return [a0, 0, bbar, bbar, kap, kap, nu, -nu]
CF = {'V2': form((r + C)/(2*D), (C - 1)/(2*D), (C + r + 2*s**2)/(4*D), (r + 2 - C)/(4*D)),
      'V4': form((C - r)/(2*D), (C + 1)/(2*D), (-C + r + 2*s**2)/(4*D), (r + 2 + C)/(4*D))}
d = json.load(open('thm3_data.json')); loc = {'s': S_, 'C': C_, 'r': R_}
allok = True
for n, cf in CF.items():
    g = [sp.sympify(x, locals=loc) for x in d['vertices'][n]['g']]
    ok = all(from_sympy(sp.together(gi - ci)).iszero() for gi, ci in zip(g, cf))
    print('%s data matches closed form used in v24_conditions: %s' % (n, ok)); allok &= ok
print('V2/V4 DATA = PRINTED CLOSED FORMS (exact):', allok)
