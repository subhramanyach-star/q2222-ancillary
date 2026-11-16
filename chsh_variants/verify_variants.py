#!/usr/bin/env python3
"""Stand-alone check of Proposition (prop:variants)(iii).

At each of the 112 efficiencies of the continuum rate certificate, the class-2
(three-negative) NPA upper bound is strictly below a certified lower bound on the
canonical CHSH maximum S_max(eta) of Theorem smax (Gigena et al. formula).
The efficiency set is required to be exactly the rate-verification grid.

Inputs: class2_certs.json and ../lossycone/data/t4/<eta>.json.
USAGE: python3 verify_variants.py
Exit status 0 if and only if ALL OK.
"""
import os, sys, json
from fractions import Fraction as F

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

H = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(H, '..')
sys.path[:0] = [os.path.join(A, 'rigorous_rates')]
from verify_all import check_cert
from loss_exact import loss_functional_F, G_CHSH_F, class2_functional
from smax_exact import S_max_bounds
from rate_grid import EXPECTED_ETA_STRS, EXPECTED_ETAS, TABLE_ETAS, ENDPOINTS
from arb_bounds import format_lower, aq, definitely_gt


def check_variants(cert_path=None, t4_dir=None, expected=None, quiet=False):
    cert_path = cert_path or os.path.join(H, 'class2_certs.json')
    t4_dir = t4_dir or os.path.join(A, 'lossycone', 'data', 't4')
    expected = tuple(expected) if expected is not None else EXPECTED_ETA_STRS
    if not os.path.isfile(cert_path):
        return dict(ok=False, why='missing class2_certs.json', n=0)
    rows = json.load(open(cert_path))
    etas = [r['eta'] for r in rows]
    if len(etas) != len(set(etas)):
        return dict(ok=False, why='duplicate efficiencies in class2_certs.json', n=len(etas))
    if set(etas) != set(expected):
        missing = [s for s in expected if s not in set(etas)]
        extra = [s for s in etas if s not in set(expected)]
        return dict(ok=False, why='efficiency set is not the rate-verification grid (missing %s extra %s)' % (missing, extra), n=len(etas))
    # endpoints and table efficiencies are included once the sets match
    ok_all = True
    worst = None
    fails = []
    lines = []
    for r in sorted(rows, key=lambda x: F(x['eta'])):
        e = F(r['eta'])
        t4 = os.path.join(t4_dir, r['eta'] + '.json')
        if not os.path.isfile(t4):
            ok_all = False
            fails.append(r['eta'])
            lines.append('%s  FAIL missing t4 certificate' % r['eta'])
            continue
        st = json.load(open(t4))
        c0, ge = loss_functional_F(G_CHSH_F, e)
        b1, _ = check_cert(st['chsh_cert'], ge)
        S1 = c0 + b1  # NPA upper on the canonical score (not used as the cap)
        c02, c2 = class2_functional(e)
        if [F(x) for x in r['cert']['c']] != c2:
            ok_all = False
            fails.append(r['eta'])
            lines.append('%s  FAIL certificate is for a different functional' % r['eta'])
            continue
        b2, _ = check_cert(r['cert'], c2)
        S2 = c02 + b2
        SL, SU = S_max_bounds(e)
        # Honest comparison: class-2 NPA upper versus a certified lower bound
        # on S_max(eta), enforcing the manuscript gap 0.00319 in exact rationals.
        gap = SL - S2
        ok = gap >= F('0.00319')
        ok_all &= ok
        if worst is None or gap < worst[0]:
            worst = (gap, r['eta'], S2, SL, SU, S1)
        if not ok:
            fails.append(r['eta'])
        if not quiet:
            print('%s  S_max_lo = %.8f  class-2 S_up = %.8f  canonical NPA S_up = %.8f  %s' % (
                r['eta'], float(SL), float(S2), float(S1), 'OK' if ok else 'FAIL',
            ))
        lines.append(r['eta'])
    if worst is None:
        return dict(ok=False, why='no rows', n=0)
    min_gap, wet, S2w, SLw, SUw, S1w = worst
    min_gap_str = str(min_gap) if min_gap <= 0 else format_lower(aq(min_gap), 9)
    out = dict(
        ok=ok_all, n=len(rows), fails=fails,
        worst_eta=wet, min_gap=min_gap, min_gap_str=min_gap_str,
        why=None if ok_all else 'class-2 NPA upper is not at least 0.00319 below certified S_max at %s' % fails,
    )
    if not quiet:
        print('min certified gap S_max_lo - S_class2_up = %s at eta = %s' % (min_gap_str if min_gap > 0 else str(float(min_gap)), wet))
        print('ALL OK =', ok_all)
    return out


if __name__ == '__main__':
    rec = check_variants(quiet=False)
    sys.exit(0 if rec['ok'] else 1)
