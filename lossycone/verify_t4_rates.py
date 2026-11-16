#!/usr/bin/env python3
"""Stand-alone re-check of the RATE side of the certified-rate / every-eta theorems.

For every data/t4/<eta>.json it re-derives in exact rationals / ball arithmetic (no SDP solver):
 (1) the honest recorded behaviour p = L_eta(q) from the stored rational device (exact loss map),
     after reconstructing q from the device trigonometry and requiring it to equal point['z'];
 (2) the witness g = (sum_i lam_i M_eta^{-T} V_i)/scale with rational convex lam over the stored
     exact polygon vertices (mu >= 0 on those vertices), and omega = g.p exactly;
 (3) for each of Eve's guesses e: the exact-rational NPA dual certificate of h_Q(Pi_e + mu*gtilde)
     (exact LDL^T, via verify_all.check_cert), the exact-rational l1 error eps = |g - gtilde|_1,
     hence P_up = max_e [c0_e + beta_e + mu*(eps - omega)] in fractions.Fraction, and
     H_lo = -log2 P_up with Arb used only for the logarithm;
 (4) the CHSH cap: exact certificate for the loss-transformed CHSH functional, S_up = c0 + beta, and
     a certified upper bound on H_cap = -log2[1/2 + 1/2 sqrt(2 - S_up^2/4)];
 (5) the gains H_lo/H_cap - 1 and the interval guarantees min_i H_lo(eta_i)/H_cap(eta_{i+1}) - 1,
     all required in the pass condition, with rational efficiencies and claimed percentages.

Acceptance keeps Arb enclosures throughout; decimals are printed only after the rounding
direction is certified.  A nonempty incomplete collection of files fails.
The cone statement N_Q(q) = C(q) at each point is checked separately by vpt_t4.py.

USAGE: python3 verify_t4_rates.py [datadir]
Exit status 0 if and only if ALL OK.
"""
import sys, os, glob, json
from fractions import Fraction as Fr

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

_H = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [_H, os.path.join(_H, '..', 'rigorous_rates'), os.path.join(_H, '..', 'chsh_variants')]
from flint import arb, ctx
ctx.prec = 200
from verify_all import check_cert, guess_funcs
import exact as EX
from loss_exact import loss_functional_F, G_CHSH_F
from rate_grid import (
    EXPECTED_ETAS, EXPECTED_ETA_STRS, EXPECTED_ETA_SET, ENDPOINTS, TABLE_ETAS,
    BANDS, GRID_GAIN_MIN, TABLE_GAIN_MIN,
)
from arb_bounds import (
    aq, certified_min_lower, hmin_lower_from_P_upper, chsh_hmin_upper, LOG2,
    strictly_below, format_lower, format_upper, definitely_gt, definitely_lt, certifies_at_least,
)

def dec2(q):
    k = q * 100
    if k.denominator != 1:
        raise ValueError('band endpoint is not a two-decimal efficiency')
    n = int(k)
    return '%d.%02d' % (n // 100, n % 100)


def device_from_point(pt):
    d = {k: ([Fr(x) for x in v] if isinstance(v, list) else Fr(v)) for k, v in pt['dev'].items()}
    needed = ('C', 'S', 'ca', 'sa', 'cb', 'sb')
    if any(k not in d for k in needed):
        raise ValueError('point is missing rational device data')
    return d


def check_one(path):
    st = json.load(open(path))
    es = os.path.basename(path)[:-5]
    e = Fr(es)
    pt = st['point']
    q_file = [Fr(x) for x in pt['z']]
    if len(q_file) != 8:
        raise ValueError('Invalid Bell-vector dimension for point[z]')
    dev = device_from_point(pt)
    an = EX.analyse(dev)  # unit-circle identities, active cones, rank-two stationarity, mu >= 0
    if an['z'] != q_file:
        raise ValueError("point['z'] does not match the rational device reconstruction")
    q = an['z']
    M, m0 = EX.lossmap(e)
    p = [sum((M[i][j] * q[j] for j in range(8)), Fr(0)) + m0[i] for i in range(8)]
    file_verts = pt['verts']
    an_mus = sorted(tuple(V['mu']) for V in an['verts'])
    file_mus = sorted(tuple(Fr(x) for x in V['mu']) for V in file_verts)
    if an_mus != file_mus:
        raise ValueError('stored polygon vertices do not match the exact enumeration')
    if len(an['gam']) != 4 or any(len(gk) != 8 for gk in an['gam']):
        raise ValueError('Invalid generator dimensions')
    Gs = []
    for V in file_verts:
        muV = [Fr(x) for x in V['mu']]
        gv = [Fr(x) for x in V['g']]
        if len(muV) != 4 or len(gv) != 8:
            raise ValueError('Invalid vertex dimensions')
        if not all(m >= 0 for m in muV):
            raise ValueError('stored vertex has a negative cone weight')
        expected_g = [
            sum((muV[k] * an['gam'][k][j] for k in range(4)), Fr(0))
            for j in range(8)
        ]
        if gv != expected_g:
            raise ValueError('Stored vertex coefficients do not match its weights')
        # Push the reconstructed supporting-family vector, not the stored g.
        Gs.append(EX.push_functional(expected_g, e))
    lam = [Fr(x) for x in st['lam']]
    if not (all(l >= 0 for l in lam) and sum(lam) == 1 and len(lam) == len(Gs)):
        raise ValueError('witness weights are not a convex combination of the pushed vertices')
    g0 = [sum((l * G_[j] for l, G_ in zip(lam, Gs)), Fr(0)) for j in range(8)]
    g = [Fr(x) for x in st['g']]
    if len(g) != 8:
        raise ValueError('Invalid Bell-vector dimension for the witness')
    sc = g0[0] / g[0] if g[0] != 0 else None
    if not (sc is not None and sc > 0 and all(g0[j] == sc * g[j] for j in range(8))):
        raise ValueError('witness not in the pushed cone')
    om = EX.dot(g, p)
    if om != Fr(st['omega']):
        raise ValueError('score mismatch')
    u = st['upper']
    mu = Fr(u['mu'])
    if mu < 0:
        raise ValueError('Negative threshold multiplier')
    gt = [Fr(x) for x in u['gtilde']]
    if len(gt) != 8:
        raise ValueError('Invalid Bell-vector dimension for gtilde')
    Gfn = guess_funcs(('A', 1))
    if len(u['per_e']) != len(Gfn):
        raise ValueError('Missing or extra guessing certificates')
    # Exact-rational P_up: eps and the dual terms stay in Fraction; Arb is used
    # only for the logarithm and for certified decimal comparisons.
    eps = sum((abs(g[j] - gt[j]) for j in range(8)), Fr(0))
    terms = []
    for (c0e, v, lab), ce in zip(Gfn, u['per_e']):
        c = [v[j] + mu * gt[j] for j in range(8)]
        beta, L = check_cert(ce['cert'], c)
        terms.append(c0e + beta + mu * (eps - om))
    Pup_exact = max(terms)
    if Pup_exact <= 0 or Pup_exact >= 1:
        raise ValueError('exact P_up is not in (0, 1)')
    Hlo = hmin_lower_from_P_upper(aq(Pup_exact))
    cj = st['chsh_cert']
    c0, ge = loss_functional_F(G_CHSH_F, e)
    if len(ge) != 8:
        raise ValueError('Invalid Bell-vector dimension for the CHSH functional')
    beta, L = check_cert(cj, ge)
    S = c0 + beta
    if 'S_up' in st and S != Fr(st['S_up']):
        raise ValueError('recomputed S_up does not match the stored canonical bound')
    Hcap = chsh_hmin_upper(S)
    if not definitely_gt(Hlo, Hcap):
        raise ValueError('certified H_lo is not strictly above H_cap')
    gain = strictly_below(Hlo / Hcap - 1)
    if not definitely_gt(gain, arb(0)):
        raise ValueError('certified grid-point gain is not positive')
    return dict(eta=e, eta_s=es, Hlo=Hlo, Hcap=Hcap, gain=gain,
                Pup_exact=Pup_exact, S=S, eps=eps, terms=terms, mu=mu, omega=om, g=g)



def _report_eta0900(rec):
    """Exact-rational sample check at eta = 0.900, including the quoted formula."""
    Pup = rec['Pup_exact']
    claimed_beta = Fr(3315405941619, 1099511627776)
    claimed_mu = Fr(12436015, 6959753)
    formula = Fr(1, 2) + claimed_beta + claimed_mu * (rec['eps'] - rec['omega'])
    if rec['mu'] != claimed_mu:
        raise ValueError('eta=0.900 mu is not the quoted 12436015/6959753')
    if formula != Pup:
        raise ValueError('eta=0.900 exact P_up does not match the quoted formula')
    Hlo_thr = Fr('0.198626')
    Hcap_thr = Fr('0.168313')
    # Direct Arb comparison of -log2(P_up_exact) and the CHSH cap, plus the
    # certified outward bounds already used for acceptance.
    log_val = -(aq(Pup).log() / LOG2)
    ok_lo = definitely_gt(log_val, aq(Hlo_thr)) and definitely_gt(rec['Hlo'], aq(Hlo_thr))
    ok_cap = definitely_lt(rec['Hcap'], aq(Hcap_thr))
    print('  [eta=0.900 exact P_up]  formula 1/2 + 3315405941619/1099511627776 + (12436015/6959753)*(eps - omega)  matches  OK')
    print('  [eta=0.900 exact P_up]  P_up = %s / %s' % (Pup.numerator, Pup.denominator))
    print('  [eta=0.900 exact P_up]  -log2(P_up) > 0.198626  (certified H_lo >= %s)  %s' % (
        format_lower(rec['Hlo'], 6), 'OK' if ok_lo else 'FAIL'))
    print('  [eta=0.900 exact P_up]  f(S_up) < 0.168313  (certified H_cap <= %s, S_up = %s)  %s' % (
        format_upper(rec['Hcap'], 6), rec['S'], 'OK' if ok_cap else 'FAIL'))
    if not (ok_lo and ok_cap):
        raise ValueError('eta=0.900 exact-rational rate bounds miss the stated thresholds')


def require_grid(found):
    if len(found) != len(EXPECTED_ETAS):
        raise ValueError('Wrong number of efficiency files')
    if len(set(found)) != len(found):
        raise ValueError('Duplicate rational efficiencies')
    if tuple(found) != EXPECTED_ETAS:
        raise ValueError('Efficiency sequence mismatch')


def main(datadir):
    files = sorted(glob.glob(os.path.join(datadir, '*.json')))
    rows = []
    allok = True
    errors = []
    for f in files:
        es = os.path.basename(f)[:-5]
        try:
            rec = check_one(f)
            pct = format_lower(rec['gain'] * 100, 2)
            print('%s  H_lo >= %s   every CHSH <= %s   gain +%s%%  OK' % (
                rec['eta_s'], format_lower(rec['Hlo'], 6), format_upper(rec['Hcap'], 6), pct,
            ), flush=True)
            if rec['eta'] == Fr('0.9'):
                _report_eta0900(rec)
            rows.append(rec)
        except Exception as ex:
            allok = False
            errors.append((es, repr(ex)))
            print(es, 'FAIL', repr(ex)[:300], flush=True)
    rows.sort(key=lambda r: r['eta'])
    try:
        require_grid([r['eta'] for r in rows])
    except Exception as ex:
        allok = False
        errors.append(('grid', repr(ex)))
        print('GRID FAIL', ex)

    # Every adjacent pair, rational interval membership, claimed band percentages.
    print('interval guarantees (lab at eta in [eta_i, eta_{i+1}] runs the eta_i protocol with clicks discarded; CHSH capped at eta_{i+1}):')
    adj = []
    for r1, r2 in zip(rows, rows[1:]):
        gain_i = strictly_below(r1['Hlo'] / r2['Hcap'] - 1)
        adj.append((r1['eta'], r2['eta'], gain_i))
        if not definitely_gt(gain_i, arb(0)):
            allok = False
            errors.append(('adjacent', 'non-positive certified separation on [%s, %s]' % (r1['eta_s'], r2['eta_s'])))
            print('  adjacent [%s, %s] FAIL: certified gain is not positive' % (r1['eta_s'], r2['eta_s']))
    for a, b, claimed in BANDS:
        band = [g for e1, e2, g in adj if a <= e1 and e2 <= b]
        if not band:
            allok = False
            errors.append(('band', 'empty band [%s, %s]' % (a, b)))
            print('  [%s, %s]: FAIL empty' % (a, b))
            continue
        gmin = certified_min_lower(band)
        ok_band = certifies_at_least(gmin, claimed)
        if not ok_band:
            allok = False
            errors.append(('band', 'claimed %.1f%% not certified on [%s, %s]' % (float(100 * claimed), a, b)))
        print('  [%s, %s]: gain >= %s%%  (%d sub-intervals)  %s' % (
            dec2(a), dec2(b), format_lower(gmin * 100, 2), len(band), 'OK' if ok_band else 'FAIL',
        ))

    if rows:
        grid_min = certified_min_lower([r['gain'] for r in rows])
        # The published 18.0 percent is the displayed upper end of the grid-point
        # range, not a lower guarantee; we only require the lower end 8.8 percent.
        ok_grid = certifies_at_least(grid_min, GRID_GAIN_MIN)
        if not ok_grid:
            allok = False
            errors.append(('grid-gains', 'grid-point gain lower bound 8.8%% not certified'))
        print('grid-point gains: %s%% .. (lower end required 8.8%%)  %s' % (
            format_lower(grid_min * 100, 2), 'OK' if ok_grid else 'FAIL',
        ))
        # Six table efficiencies: lower percentage 9.2 percent (do not replace 8.8).
        table = [r for r in rows if r['eta'] in TABLE_ETAS]
        if len(table) != len(TABLE_ETAS):
            allok = False
            errors.append(('table', 'missing table efficiencies'))
            print('table efficiencies FAIL: expected 6, found %d' % len(table))
        else:
            tmin = certified_min_lower([r['gain'] for r in table])
            ok_tab = certifies_at_least(tmin, TABLE_GAIN_MIN)
            if not ok_tab:
                allok = False
                errors.append(('table', 'table-efficiency gain lower bound 9.2%% not certified'))
            print('table-efficiency gains: >= %s%%  (required 9.2%%)  %s' % (
                format_lower(tmin * 100, 2), 'OK' if ok_tab else 'FAIL',
            ))
    else:
        allok = False
        errors.append(('grid', 'no points'))

    # CHSH-variant coverage at every required endpoint (in fact on the whole grid).
    try:
        from verify_variants import check_variants
        var = check_variants(
            t4_dir=datadir,
            expected=EXPECTED_ETA_STRS,
            quiet=True,
        )
        if not var['ok']:
            allok = False
            errors.append(('variants', var.get('why', 'variant check failed')))
            print('CHSH-variant coverage FAIL:', var.get('why', var))
        else:
            print('CHSH-variant coverage: class-2 NPA upper < certified S_max at all %d grid points (min gap %s at eta=%s)  OK' % (
                var['n'], var['min_gap_str'], var['worst_eta'],
            ))
    except Exception as ex:
        allok = False
        errors.append(('variants', repr(ex)))
        print('CHSH-variant coverage FAIL', repr(ex)[:300])

    print('ALL OK =', allok, '(%d points)' % len(rows))
    if errors and not allok:
        print('failures:', len(errors))
    return 0 if allok else 1


if __name__ == '__main__':
    D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(_H, 'data', 't4')
    sys.exit(main(D))
