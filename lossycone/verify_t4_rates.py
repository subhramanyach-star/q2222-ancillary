#!/usr/bin/env python3
"""Stand-alone re-check of the RATE side of Theorems (certified rates, every-eta interval) of prl_v16, App. sec:lossy.
For every data/t4/<eta>.json it re-derives in exact rationals / ball arithmetic (no SDP solver):
 (1) the honest recorded behaviour p = L_eta(q) from the stored rational device point q (exact loss map);
 (2) the witness g = (sum_i lam_i M_eta^{-T} V_i)/scale with rational convex lam over the stored exact polygon vertices,
     and omega = g.p exactly;
 (3) for each of Eve's guesses e: the exact-rational NPA dual certificate of h_Q(Pi_e + mu*gtilde) (exact LDL^T, via
     verify_all.check_cert), the l1 error |g - gtilde|_1 in ball arithmetic, hence
     P_guess <= max_e [c0_e + beta_e + mu*err] - mu*omega  and  H_lo = -log2 P_up (Lemma 'Two-sided bounds'(a));
 (4) the CHSH cap: exact certificate for the loss-transformed CHSH functional, S_up = c0 + beta, and
     H_cap = -log2[1/2 + 1/2 sqrt(2 - S_up^2/4)] rounded up;
 (5) the gains H_lo/H_cap - 1 and the interval guarantees min_i H_lo(eta_i)/H_cap(eta_{i+1}) - 1.
The cone statement N_Q(q) = C(q) at each point is checked separately by vpt.py.
USAGE: python3 verify_t4_rates.py [datadir]"""
import sys, os, glob, json
from fractions import Fraction as Fr
_H = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [_H, os.path.join(_H, '..', 'rigorous_rates'), os.path.join(_H, '..', 'eta_study')]
from flint import arb, ctx
ctx.prec = 200
from verify_all import check_cert, guess_funcs
import exact as EX
from cert_eta import loss_functional_F, G_CHSH_F
def aq(q): q = Fr(q); return arb(q.numerator)/q.denominator
def up(b): return float(b.mid()) + float(b.rad())
def lo(b): return float(b.mid()) - float(b.rad())
LOG2 = arb(2).log()
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 't4')
rows = []; allok = True
for f in sorted(glob.glob(os.path.join(D, '*.json'))):
    st = json.load(open(f)); es = os.path.basename(f)[:-5]; e = Fr(es)
    try:
        pt = st['point']; q = [Fr(x) for x in pt['z']]; M, m0 = EX.lossmap(e)
        p = [sum((M[i][j]*q[j] for j in range(8)), Fr(0)) + m0[i] for i in range(8)]
        Gs = [EX.push_functional([Fr(x) for x in V['g']], e) for V in pt['verts']]
        lam = [Fr(x) for x in st['lam']]; assert all(l >= 0 for l in lam) and sum(lam) == 1 and len(lam) == len(Gs)
        g0 = [sum((l*G_[j] for l, G_ in zip(lam, Gs)), Fr(0)) for j in range(8)]
        g = [Fr(x) for x in st['g']]; sc = g0[0]/g[0] if g[0] != 0 else None
        assert sc is not None and sc > 0 and all(g0[j] == sc*g[j] for j in range(8)), 'witness not in the pushed cone'
        om = EX.dot(g, p); assert om == Fr(st['omega']), 'score mismatch'
        u = st['upper']; mu = Fr(u['mu']); gt = [Fr(x) for x in u['gtilde']]
        err = sum((abs(aq(x) - aq(y)) for x, y in zip(g, gt)), arb(0))
        terms = []
        for (c0, v, lab), ce in zip(guess_funcs(('A', 1)), u['per_e']):
            c = [v[j] + mu*gt[j] for j in range(8)]
            beta, L = check_cert(ce['cert'], c)
            terms.append(aq(c0) + aq(beta) + aq(mu)*err - aq(mu)*aq(om))
        Pup = max(up(t) for t in terms)
        Hlo = lo(-(aq(Fr(Pup))).log()/LOG2) if Pup < 1 else 0.0
        cj = st['chsh_cert']; c0, ge = loss_functional_F(G_CHSH_F, e)
        beta, L = check_cert(cj, ge); S = c0 + beta; Sa = aq(S)
        Hc = -((arb(1)/2 + (arb(2) - Sa*Sa/4).sqrt()/2).log()/LOG2); Hcap = up(Hc)
        ok = Hlo > Hcap and abs(Hlo - st['H_lo']) < 1e-9
        rows.append((float(e), Hlo, Hcap, ok)); allok &= ok
        print('%.3f  H_lo >= %.6f   every CHSH <= %.6f   gain %+.2f%%  %s' % (float(e), Hlo, Hcap, 100*(Hlo/Hcap - 1), 'OK' if ok else 'FAIL'), flush=True)
    except Exception as ex:
        allok = False; print(es, 'FAIL', repr(ex)[:200], flush=True)
rows.sort()
bands = [(0.75, 0.80), (0.80, 0.85), (0.85, 0.90), (0.90, 0.95), (0.95, 0.98)]
print('interval guarantees (lab at eta in [eta_i, eta_{i+1}] runs the eta_i protocol with clicks discarded; CHSH capped at eta_{i+1}):')
for a, b in bands:
    gs = [r1[1]/r2[2] - 1 for r1, r2 in zip(rows, rows[1:]) if a - 1e-12 <= r1[0] and r2[0] <= b + 1e-12]
    if gs: print('  [%.2f, %.2f]: gain >= %.2f%%  (%d sub-intervals)' % (a, b, 100*min(gs), len(gs)))
print('grid-point gains: %.2f%% .. %.2f%%' % (100*min(r[1]/r[2] - 1 for r in rows), 100*max(r[1]/r[2] - 1 for r in rows)))
print('ALL OK =', allok, '(%d points)' % len(rows))
