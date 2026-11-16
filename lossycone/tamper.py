"""Negative tests for vpt.py: perturbed certificates must be rejected."""
import json, sys, os
from fractions import Fraction as Fr
sys.argv = ['vpt.py']
src = open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vpt.py')).read().split("if __name__ == '__main__':")[0]; exec(src)
res = []
for e, k, j, i, dlt in [('0.8', 0, 0, 0, Fr(1, 10)), ('0.8', 2, 3, 1, Fr(1, 20)), ('0.98', 5, 1, 0, Fr(-1, 10)), ('0.95', 0, 2, 1, Fr(1, 5)), ('0.75', 3, 0, 1, Fr(1, 50))]:
    pt = P[e]; z, c0, G, m, verts, gs = rederive(pt); V = pt['verts'][k]
    mu = [Fr(x) for x in V['mu']]; g = [Fr(x) for x in V['g']]; cert = json.loads(json.dumps(CERT['%s/%d' % (e, k)]))
    cert['R'][j][i] = str(Fr(cert['R'][j][i]) + dlt)
    d0, d1 = forced(G, mu, c0, g); cert = dict(cert, _d0=d0, _d1=d1)
    try: r = verify_vertex(g, mu, c0, cert, maxbox=20000)
    except AssertionError as ex: r = dict(ok=False, why='assert: ' + str(ex)[:80])
    print(e, k, 'R[%d][%d] += %s' % (j, i, dlt), '->', 'REJECTED' if not r['ok'] else 'ACCEPTED(!)', {a: b for a, b in r.items() if a != 'tests'}, flush=True)
# perturbed functional: g scaled towards a non-cone direction (stationarity broken)
e = '0.8'; pt = P[e]; z, c0, G, m, verts, gs = rederive(pt); V = pt['verts'][0]
mu = [Fr(x) for x in V['mu']]; mu2 = [mu[0] + Fr(1, 100), mu[1], mu[2], mu[3] - Fr(1, 100)]
g2 = [sum((mu2[a]*G[a][b] for a in range(4)), Fr(0)) for b in range(8)]
try: d0, d1 = forced(G, mu2, c0, g2); print('non-stationary functional: forced data exist (!)')
except AssertionError: print('non-stationary functional: REJECTED (first-order equalities inconsistent)')
