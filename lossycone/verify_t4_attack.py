#!/usr/bin/env python3
"""Verify the stored T4 attacks that support the table's entropy upper endpoints.

verify_t4_rates.py checks randomness lower bounds, not these attacks.
verify_all.verify() expects a different schema (case, spec, ...).

This checker reuses:
  - verify_t4_rates.check_one(path) for the reconstructed witness g and score omega
  - verify_all.check_attack_input, zvec_arb, and guess_funcs
It recomputes the mixture zbar, certifies g · zbar >= omega, and bounds the
guessing success from below.  That upper-entropy endpoint is independent of
the separation proved by the lower bounds.

USAGE: python3 verify_t4_attack.py [file-or-datadir ...]
With no arguments, checks every data/t4/<eta>.json and requires the full 112-point grid.
Exit status 0 if and only if every requested file passes (and, for a directory of
the default grid, the collection is the complete EXPECTED_ETAS set).
"""
if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

import sys, os, glob, json
from fractions import Fraction as Fr
from flint import arb, ctx
ctx.prec = 200

_H = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_H, os.path.join(_H, '..', 'rigorous_rates')]

from verify_t4_rates import check_one, require_grid, EXPECTED_ETAS
from verify_all import check_attack_input, zvec_arb, guess_funcs
from families_arb import dot
from arb_bounds import (
    aq, strictly_below, certified_float_lower, hmin_upper_from_P_lower,
    format_upper, definitely_ge,
)

# Same Alice-input guess used by the T4 rate verifier.
GUESS = ('A', 1)

# Upper endpoints of the selected-expression intervals in Table tab:lossy
# of the accompanying paper, at the six reference points.  These are the
# paper's outward-rounded six-decimal endpoints, hard-coded from that table.
# The attack path certifies these upper ends.  Do not invent unpublished
# numeric endpoints for the other 106 grid points.
PUBLISHED_HHI = {
    Fr('0.75'): Fr('0.009847'),
    Fr('0.80'): Fr('0.037864'),
    Fr('0.85'): Fr('0.095105'),
    Fr('0.90'): Fr('0.198627'),
    Fr('0.95'): Fr('0.385566'),
    Fr('0.98'): Fr('0.597368'),
}


def hhi_is_finite(Hhi):
    return Hhi is not None and hasattr(Hhi, 'is_finite') and bool(Hhi.is_finite())


def frac_as_sixdec(q):
    """Display a six-decimal table endpoint without rounding."""
    q = Fr(q)
    k = q * 10**6
    if k.denominator != 1:
        return str(q)
    n = int(k)
    sign = '-' if n < 0 else ''
    n = abs(n)
    return '%s%d.%06d' % (sign, n // 10**6, n % 10**6)


def published_upper_ok(Hhi, eta):
    """Rigorous test that the certified H_hi ball is <= the v20 table upper, if any."""
    published = PUBLISHED_HHI.get(Fr(eta))
    if published is None:
        return None, None
    if not hhi_is_finite(Hhi):
        return False, published
    # definitely_ge(bound, Hhi) means the entire enclosure of Hhi is <= bound.
    return bool(definitely_ge(aq(published), Hhi)), published


def decide_attack_ok(score_ok, Hhi, eta):
    """Acceptance: feasible score, finite entropy upper bound, and (for the
    six table efficiencies) H_hi <= the published v20 upper rational.
    A missing or non-finite H_hi is a failure even if the score holds."""
    table_ok, published = published_upper_ok(Hhi, eta)
    attack_ok = bool(score_ok) and hhi_is_finite(Hhi)
    if table_ok is not None:
        attack_ok = attack_ok and bool(table_ok)
    return attack_ok, table_ok, published


def check_attack(path):
    rec = check_one(path)
    st = json.load(open(path))
    if 'attack' not in st:
        raise ValueError('no stored attack')
    A = st['attack']
    if 't' in A:
        raise ValueError('T4 attack unexpectedly includes a local remainder t')
    G = guess_funcs(GUESS)
    kind = 'single'
    pis = check_attack_input(A, G, kind)
    zs = [zvec_arb(st_) for st_ in A['strategies']]
    zbar = [sum((aq(pi) * z[j] for pi, z in zip(pis, zs)), arb(0)) for j in range(8)]
    g = [aq(x) for x in rec['g']]
    omega = aq(rec['omega'])
    score = dot(g, zbar)
    # Certify g · zbar >= omega by a fully nonnegative enclosure of the difference.
    gap = score - omega
    score_ok = definitely_ge(score, omega) or bool(gap >= 0)
    gi = {lab: (c0, v) for c0, v, lab in G}
    Pm = sum(
        (aq(pi) * (aq(gi[tuple(e)][0]) + dot([aq(x) for x in gi[tuple(e)][1]], z))
         for pi, e, z in zip(pis, A['guesses'], zs)),
        arb(0),
    )
    if not score_ok:
        Plo = None
        Hhi = None
    else:
        Plo_ball = strictly_below(Pm)
        Plo = certified_float_lower(Plo_ball)
        Hhi = hmin_upper_from_P_lower(Plo_ball)
    attack_ok, table_ok, published = decide_attack_ok(score_ok, Hhi, rec['eta'])
    return dict(
        eta_s=rec['eta_s'], eta=rec['eta'],
        score_ok=bool(score_ok), attack_ok=attack_ok,
        table_ok=table_ok, published_upper=str(published) if published is not None else None,
        score_minus_omega=gap,
        P_lo=Plo,
        H_hi=Hhi,
        H_hi_dec=format_upper(Hhi, 6) if hhi_is_finite(Hhi) else None,
        n_cols=len(pis),
    )


def collect_targets(args):
    """Return (files, require_full_grid)."""
    if not args:
        return sorted(glob.glob(os.path.join(_H, 'data', 't4', '*.json'))), True
    files, require_grid_flag = [], False
    for a in args:
        if os.path.isdir(a):
            files.extend(sorted(glob.glob(os.path.join(a, '*.json'))))
            require_grid_flag = True
        else:
            files.append(a)
    return files, require_grid_flag


def main(args):
    files, need_grid = collect_targets(args)
    allok = True
    rows = []
    if not files:
        print('empty requested collection')
        print('RESULT ALL OK =', False)
        return False
    if need_grid:
        try:
            found = [Fr(os.path.basename(f)[:-5]) for f in files]
            if len(files) != 112:
                raise ValueError('expected exactly 112 T4 files, found %d' % len(files))
            require_grid(sorted(found))
            print('T4 grid: 112 files, unique rationals, exact EXPECTED_ETAS  OK', flush=True)
        except Exception as ex:
            print('GRID FAIL', ex)
            print('RESULT ALL OK =', False)
            return False
    for f in files:
        es = os.path.basename(f)[:-5]
        try:
            r = check_attack(f)
            recs = 'OK' if r['score_ok'] else 'FAIL'
            hdec = r['H_hi_dec'] if r['H_hi_dec'] is not None else 'FAIL'
            line = '%s  g·zbar >= omega  %s  H_hi <= %s  ncols=%d' % (
                r['eta_s'], recs, hdec, r['n_cols'],
            )
            if r['published_upper'] is not None:
                tstat = 'OK' if r['table_ok'] else 'FAIL'
                line += '  table upper %s  %s' % (frac_as_sixdec(r['published_upper']), tstat)
            acc = 'OK' if r['attack_ok'] else 'FAIL'
            line += '  accepted %s' % acc
            print(line, flush=True)
            # Never treat H_hi <= FAIL as success, even if the score is feasible.
            if r['H_hi_dec'] is None and r['attack_ok']:
                raise RuntimeError('internal error: accepted an attack with no finite H_hi')
            rows.append(r)
            allok &= r['attack_ok']
        except Exception as ex:
            allok = False
            print(es, 'FAIL', repr(ex)[:300], flush=True)
    if need_grid and len(rows) != len(EXPECTED_ETAS):
        allok = False
        print('GRID FAIL incomplete attack results: %d of %d' % (len(rows), len(EXPECTED_ETAS)))
    table_rows = [r for r in rows if r.get('published_upper') is not None]
    if table_rows:
        print('published v20 table upper endpoints (selected-expression intervals):')
        for r in sorted(table_rows, key=lambda x: x['eta']):
            print('  eta %s  certified H_hi <= %s  published upper %s  %s' % (
                r['eta_s'],
                r['H_hi_dec'] if r['H_hi_dec'] is not None else 'FAIL',
                frac_as_sixdec(r['published_upper']),
                'OK' if r['table_ok'] else 'FAIL',
            ))
        if any(not r['table_ok'] for r in table_rows):
            allok = False
    print('RESULT ALL OK =', bool(allok and rows))
    return bool(allok and rows)


if __name__ == '__main__':
    allok = main(sys.argv[1:])
    sys.exit(0 if allok else 1)
