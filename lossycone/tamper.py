"""Negative tests: a perturbed cone certificate must be rejected by the verifier.

The shipped JSON files are not modified.  The script exits with status 0
only if an unperturbed certificate is accepted and every perturbation is
rejected.
"""
import json
import os
import sys
from fractions import Fraction as Fr

_H = os.path.dirname(os.path.abspath(__file__))
if _H not in sys.path:
    sys.path.insert(0, _H)

from vpt import load_points_certs, rederive, forced, verify_vertex


def _perturbations():
    return [
        ('0.8', 0, 0, 0, Fr(1, 10)),
        ('0.8', 2, 3, 1, Fr(1, 20)),
        ('0.98', 5, 1, 0, Fr(-1, 10)),
        ('0.95', 0, 2, 1, Fr(1, 5)),
        ('0.75', 3, 1, 0, Fr(1, 5)),
    ]


def main():
    certs_path = os.path.join(_H, 'data', 'certs.json')
    P, CERT = load_points_certs(certs_path=certs_path)
    accepted_good = False
    rejected = 0
    leaked = 0

    # Sanity: an unperturbed certificate must still be accepted.
    e0, k0 = '0.8', 0
    pt = P[e0]
    z, c0, G, m, verts, gs = rederive(pt)
    V = pt['verts'][k0]
    mu = [Fr(x) for x in V['mu']]
    g = [Fr(x) for x in V['g']]
    cert = json.loads(json.dumps(CERT['%s/%d' % (e0, k0)]))
    d0, d1 = forced(G, mu, c0, g)
    cert = dict(cert, _d0=d0, _d1=d1)
    try:
        r = verify_vertex(g, mu, c0, cert, maxbox=20000)
    except AssertionError as ex:
        r = dict(ok=False, why='assert: ' + str(ex)[:80])
    if r.get('ok'):
        accepted_good = True
        print('unperturbed %s/%d: ACCEPTED' % (e0, k0), flush=True)
    else:
        print('unperturbed %s/%d: REJECTED(!)' % (e0, k0), r, flush=True)

    for e, k, j, i, dlt in _perturbations():
        pt = P[e]
        z, c0, G, m, verts, gs = rederive(pt)
        V = pt['verts'][k]
        mu = [Fr(x) for x in V['mu']]
        g = [Fr(x) for x in V['g']]
        cert = json.loads(json.dumps(CERT['%s/%d' % (e, k)]))
        cert['R'][j][i] = str(Fr(cert['R'][j][i]) + dlt)
        d0, d1 = forced(G, mu, c0, g)
        cert = dict(cert, _d0=d0, _d1=d1)
        try:
            r = verify_vertex(g, mu, c0, cert, maxbox=20000)
        except AssertionError as ex:
            r = dict(ok=False, why='assert: ' + str(ex)[:80])
        ok = bool(r.get('ok'))
        if ok:
            leaked += 1
            tag = 'ACCEPTED(!)'
        else:
            rejected += 1
            tag = 'REJECTED'
        brief = {a: b for a, b in r.items() if a != 'tests'}
        print(e, k, 'R[%d][%d] += %s' % (j, i, dlt), '->', tag, brief, flush=True)

    # Perturbed functional: stationarity is broken.
    e = '0.8'
    pt = P[e]
    z, c0, G, m, verts, gs = rederive(pt)
    V = pt['verts'][0]
    mu = [Fr(x) for x in V['mu']]
    mu2 = [mu[0] + Fr(1, 100), mu[1], mu[2], mu[3] - Fr(1, 100)]
    g2 = [sum((mu2[a] * G[a][b] for a in range(4)), Fr(0)) for b in range(8)]
    stat_ok = False
    try:
        forced(G, mu2, c0, g2)
        print('non-stationary functional: forced data exist (!)')
    except AssertionError:
        stat_ok = True
        print('non-stationary functional: REJECTED (first-order equalities inconsistent)')

    allok = accepted_good and leaked == 0 and rejected == len(_perturbations()) and stat_ok
    print('RESULT TAMPER ALL OK =', allok)
    return 0 if allok else 1


if __name__ == '__main__':
    sys.exit(main())
