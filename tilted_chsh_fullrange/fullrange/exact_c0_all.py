"""Exact identities behind test T2 (forced touching at c0), for V1, V3, Wa, Wb, in K = Q(s)(C, r):
for each sign pair (s0,s1): A_{s0}(c0) = 4 mu_{s0|0}^2, B_{s1}(c0) = 4 mu_{s1|1}^2, T_{s0s1} = 2 mu_{s0|0} + 2 mu_{s1|1},
and A'_{s0}(c0)/mu_{s0|0} + B'_{s1}(c0)/mu_{s1|1} = 0, where ' = d/dc along delta(c) = d0 + d1 (c - c0) + O((c-c0)^2).
With mu > 0 (comb.py, exact) these give F = T - sqrt A - sqrt B = 0 and F_c = -A'/(4mu) - B'/(4mu) = 0 at c = c0."""
import json, sympy as sp, time
from kfield import K, from_sympy
S_, C_, R_ = sp.symbols('s C r', positive=True)
loc = {'s': S_, 'C': C_, 'r': R_}
SRC = {'V1': 'thm3_data.json', 'V3': 'thm3_data.json',
       'Wa': 'oct_data.json', 'Wb': 'oct_data.json'}
KEYS = [(1, 0), (1, 1), (-1, 0), (-1, 1)]
allok = True
for name, path in SRC.items():
    J = json.load(open(path)); V = J['vertices'][name]
    E = {k: [from_sympy(sp.sympify(e, locals=loc)) for e in V[k]] for k in ('g', 'mu', 'd0', 'd1')}
    c0 = from_sympy(sp.sympify(J['c0'], locals=loc)); g, mu, d0, d1 = E['g'], E['mu'], E['d0'], E['d1']
    M = {k: mu[i] for i, k in enumerate(KEYS)}
    res = []
    for s0 in (1, -1):
        for s1 in (1, -1):
            u = [g[2 + i] * K(sp.Rational(1, 2)) + d0[i] + g[4 + i] * K(s0) for i in range(2)]
            w = [g[2 + i] * K(sp.Rational(1, 2)) - d0[i] + g[6 + i] * K(s1) for i in range(2)]
            A = u[0] * u[0] + u[1] * u[1] + K(2) * c0 * u[0] * u[1]; B = w[0] * w[0] + w[1] * w[1] + K(2) * c0 * w[0] * w[1]
            # d/dc A along delta: dA = 2 u.G u' + 2 u0 u1 with u' = d1 ; dB: w' = -d1
            dA = K(2) * (u[0] * d1[0] + u[1] * d1[1] + c0 * (u[0] * d1[1] + u[1] * d1[0])) + K(2) * u[0] * u[1]
            dB = K(-2) * (w[0] * d1[0] + w[1] * d1[1] + c0 * (w[0] * d1[1] + w[1] * d1[0])) + K(2) * w[0] * w[1]
            T = K(1) - K(s0) * g[0] - K(s1) * g[1]
            m0, m1 = M[(s0, 0)], M[(s1, 1)]
            ids = [A - K(4) * m0 * m0, B - K(4) * m1 * m1, T - K(2) * m0 - K(2) * m1, dA * m1 + dB * m0]
            ok = all(x.iszero() for x in ids); allok &= ok; res.append(((s0, s1), ok))
    print(name, res, flush=True)
print('ALL c0 IDENTITIES EXACT:', allok)
