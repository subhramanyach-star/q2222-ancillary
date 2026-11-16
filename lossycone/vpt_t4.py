#!/usr/bin/env python3
"""Stand-alone re-check of the supporting-cone certificates at all certified detection efficiencies:
runs the per-point checks of vpt.py (exact re-derivation of q, c0, active cones, exhaustive vertex enumeration of
C(q), forced certificate data, exact first-order/endpoint checks, Arb bisection over c in [-1,1]) on every data/t4/<eta>.json.

The T4 grid is validated before any vertex certificate is checked: exactly 112 files, no duplicate
rational efficiencies, exact agreement with EXPECTED_ETAS (same require_grid as the rate verifier),
and at least one reconstructed vertex at each point.  An incomplete nonempty directory does not
report success.  vpt.py is imported normally; this wrapper does not exec vpt.py or load data/points.json.

USAGE: python3 vpt_t4.py [datadir]
Exit status 0 if and only if the full grid is present and every vertex certificate is accepted.
"""
if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

import sys, os, glob, json
from fractions import Fraction as Fr

_H = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_H, os.path.join(_H, '..', 'rigorous_rates')]

import vpt
from verify_t4_rates import require_grid


def load_t4(datadir):
    files = sorted(glob.glob(os.path.join(datadir, '*.json')))
    P, CERT = {}, {}
    for f in files:
        e = os.path.basename(f)[:-5]
        st = json.load(open(f))
        P[e] = st['point']
        for k, c in st['certs'].items():
            CERT['%s/%s' % (e, k)] = c
    return P, CERT, files


def check_grid_and_vertices(P, nfiles):
    """Grid enforcement and a nonempty reconstructed vertex list, before certificates."""
    if nfiles != 112:
        raise ValueError('expected exactly 112 T4 files, found %d' % nfiles)
    found = [Fr(e) for e in P]
    if len(found) != nfiles:
        raise ValueError('number of loaded points does not match the number of files')
    require_grid(sorted(found))
    for e, pt in P.items():
        z, c0, G, m, verts, gs = vpt.rederive(pt)
        if not verts:
            raise ValueError('no reconstructed vertex at eta %s' % e)
        if not pt.get('verts'):
            raise ValueError('no stored vertex at eta %s' % e)


def main(datadir):
    P, CERT, files = load_t4(datadir)
    try:
        check_grid_and_vertices(P, len(files))
        print('T4 grid: 112 files, unique rationals, exact EXPECTED_ETAS, vertices at every point  OK', flush=True)
    except Exception as ex:
        print('GRID FAIL', ex)
        print('RESULT ALL OK =', False)
        return False
    etas = sorted(P, key=lambda s: Fr(s))
    return vpt.run(P, CERT, etas)


if __name__ == '__main__':
    datadir = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_H, 'data', 't4')
    allok = main(datadir)
    sys.exit(0 if allok else 1)
