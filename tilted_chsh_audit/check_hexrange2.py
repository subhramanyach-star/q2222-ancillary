"""Refine: decide the 108 inequalities on [pi/16, 0.2022pi] exactly (norm root counting) and locate true minima."""
import sympy as sp, mpmath as mp, numpy as np
from flint import arb, fmpq
from kfun import arb_min_on, exact_sign_on, KF
import model as M
from check_exact_helpers import vertices_and_mu
mp.mp.dps = 40
def sr(t0, t1): return (sp.Rational(int(mp.floor(mp.sin(2*t0) * 10**30)), 10**30), sp.Rational(int(mp.ceil(mp.sin(2*t1) * 10**30)), 10**30))
V, MU = vertices_and_mu()
for label, (lo, hi) in (('Theta', sr(mp.pi/8, 3*mp.pi/16)), ('wide', sr(mp.pi/16, mp.mpf('0.2022')*mp.pi))):
    res = []
    for name in ('V1', 'V2', 'V1s', 'V3', 'V4', 'V3s'):
        g = V[name]['g']
        items = [('mu%+d|%d' % ax, MU[name][ax]) for ax in M.AX]
        items += [('1-g.L%s' % (L,), 1 - M.dot(g, M.Lvec(*L))) for L in M.LPTS if L not in M.FACETS[name]]
        for lab, x in items:
            f = KF(x)
            ss = np.linspace(float(lo), float(hi), 4001)
            vals = [float(f(arb(float(s))).mid()) for s in ss]
            i = int(np.argmin(vals))
            sg, info = exact_sign_on(x, lo, hi)
            ok = sg == 1
            if not ok:
                lb, ub = arb_min_on(x, lo, hi, 8192); ok = lb > 0
                info += '; arb lb %.3g' % float(lb)
            res.append((vals[i], name, lab, ok, float(mp.asin(ss[i]) / 2 / mp.pi), info))
    res.sort()
    print(label, 'all 108 proved > 0:', all(r[3] for r in res))
    for r in res[:5]:
        print('   min %.3e  %-4s %-22s at theta=%.5f pi  proved=%s (%s)' % (r[0], r[1], r[2], r[4], r[3], r[5]))
