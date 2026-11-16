"""Exact verification (computer algebra, mod C^2 = 1-s^2, r^2 = 1+s^2) of the identities used by test T2 of bnb_v13.py:
for V in {V1, V3} and every sign pair (s0,s1), with delta(c0) = d0, delta'(c0) = d1:
   I1:  (R^2 - A - B)^2 - 4AB = 0                    at c = c0   (=> sqrt A + sqrt B = R given R >= 0, R^2-A-B >= 0)
   I2:  (dA/dc)^2 B - (dB/dc)^2 A = 0                 at c = c0   (=> F_c(c0) = 0 given (dA/dc)(dB/dc) <= 0)
Sign side-conditions R > 0, R^2-A-B >= 0, A, B > 0, A'B' <= 0 are then checked in interval arithmetic on Theta."""
import pickle, sympy as sp, sys, time
from symbolic_m2 import red, s, C, r
from arbeval import compile_expr
D = __import__('thm3_data').load()  # exact data from thm3_data.json
c0 = D['c0']
c = sp.Symbol('c'); e0, e1 = sp.symbols('e0 e1')
COMBOS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]
out = {}
for name in sys.argv[1:] or ['V1', 'V3']:
    V = D['verts'][name]; g = V['g']; d0, d1 = V['d0'], V['d1']
    for s0, s1 in COMBOS:
        t0 = time.time()
        dl = [e0, e1]
        u = [g[2 + i] / 2 + dl[i] + s0 * g[4 + i] for i in range(2)]
        w = [g[2 + i] / 2 - dl[i] + s1 * g[6 + i] for i in range(2)]
        A = u[0]**2 + u[1]**2 + 2 * c * u[0] * u[1]; B = w[0]**2 + w[1]**2 + 2 * c * w[0] * w[1]
        R = 1 - s0 * g[0] - s1 * g[1]
        dA = sp.diff(A, c) + sp.diff(A, e0) * d1[0] + sp.diff(A, e1) * d1[1]
        dB = sp.diff(B, c) + sp.diff(B, e0) * d1[0] + sp.diff(B, e1) * d1[1]
        sub = {c: c0, e0: d0[0], e1: d0[1]}
        A0, B0, dA0, dB0 = [red(x.subs(sub)) for x in (A, B, dA, dB)]
        R0 = red(R)
        I1 = red((R0**2 - A0 - B0)**2 - 4 * A0 * B0)
        I2 = red(dA0**2 * B0 - dB0**2 * A0)
        print(name, (s0, s1), 'I1 =', I1, ' I2 =', I2, ' (%.1fs)' % (time.time() - t0), flush=True)
        out[(name, s0, s1)] = dict(I1=I1, I2=I2, R=R0, A=A0, B=B0, dA=dA0, dB=dB0)
import json
json.dump({'%s|%d|%d' % k: {kk: sp.sstr(vv) for kk, vv in v.items()} for k, v in out.items()},
          open('exact_c0_%s.json' % '_'.join(sys.argv[1:] or ['V1', 'V3']), 'w'), indent=1)
