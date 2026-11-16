import os
"""Certify that under the source-noise models of Table S4 (white noise, dephasing, and the
experiment-matched white+dephasing rows), the second class of CHSH sign variants (---+) with
missed detections recorded as +1 has certified S_up no larger than the canonical CHSH.
Model (same as noise/ncert.py):
  white: S <= c0 + v * beta
  deph:  S <= c0 + p * beta + (1-p) * L_loc,  L_loc = local max of the loss-transformed functional
  wd:    S <= c0 + v * (p * beta + (1-p) * L_loc)
where beta = exact NPA upper bound on h_Q of the ideal (loss-free) coefficients of the variant.
Canonical certificates for the eight experiment rows are taken from noise/certs/wd_*.json.
Also check a small white/deph grid covering the table ranges."""
import os, sys, json, glob, time, itertools
from fractions import Fraction as F
import os
A = os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
sys.path[:0] = [os.path.join(A,'rigorous_rates'), os.path.join(A,'eta_study'), os.path.join(A,'lossycone')]
import rates_core as R
from verify_all import check_cert
from cert_eta import loss_functional_F, G_CHSH_F

def loss_fun(sv, eta):
    """Return (c0, g_ideal) for sign pattern sv with +1 binning at efficiency eta.
    Recorded score = c0 + g_ideal · z_ideal for a loss-free behaviour z (exact identity)."""
    sa = [sv[0]+sv[1], sv[2]+sv[3]]; sb = [sv[0]+sv[2], sv[1]+sv[3]]
    c0 = (1-eta)**2 * sum(sv)
    # ideal coefficients of the loss-transformed functional: recorded = eta^2 * s·E + eta(1-eta)*s·(A,B) + c0
    # so g_ideal = (eta(1-eta) sa, eta(1-eta) sb, eta^2 s); h_Q of this g_ideal is what we certify
    g = [eta*(1-eta)*sa[0], eta*(1-eta)*sa[1], eta*(1-eta)*sb[0], eta*(1-eta)*sb[1]] + [eta*eta*x for x in sv]
    return c0, g

def local_max(g):
    best = None
    for a0,a1,b0,b1 in itertools.product((1,-1), repeat=4):
        z = [a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]
        val = sum(g[j]*z[j] for j in range(8))
        best = val if best is None or val > best else best
    return best

def S_up(c0, g, model, v, p):
    beta, L, cj, qv = R.certify_hq(g)
    beta2, _ = check_cert(cj, g)  # exact re-check
    assert beta2 == beta
    Lloc = local_max(g)
    if model == 'white':
        S = c0 + v * beta
    elif model == 'deph':
        S = c0 + p * beta + (1-p) * Lloc
    elif model == 'wd':
        S = c0 + v * (p * beta + (1-p) * Lloc)
    else:
        raise ValueError(model)
    return S, beta, Lloc, cj, float(c0) + float(v if model=='white' else (v*p if model=='wd' else p))*qv + (0 if model=='white' else float((1-(v*p if model=='wd' else p))*Lloc if model!='white' else 0))

def compare(eta, model, v, p, tag):
    e = F(eta) if not isinstance(eta, F) else eta
    v = F(v) if not isinstance(v, F) else v
    p = F(p) if not isinstance(p, F) else p
    # canonical
    c0c, gc = loss_fun([1,1,1,-1], e)
    # for canonical we can also use loss_functional_F which matches G_CHSH_F
    c0c2, gc2 = loss_functional_F(G_CHSH_F, e)
    assert c0c == c0c2 and gc == gc2
    Sc, bc, Lc, cjc, _ = S_up(c0c, gc, model, v, p)
    # class 2
    c02, g2 = loss_fun([-1,-1,-1,1], e)
    S2, b2, L2, cj2, q2 = S_up(c02, g2, model, v, p)
    ok = S2 <= Sc
    row = dict(tag=tag, eta=str(e), model=model, v=str(v), p=str(p),
               S_can=float(Sc), S_c2=float(S2), margin=float(Sc-S2), ok=ok,
               cert_can=cjc, cert_c2=cj2, Lloc_can=str(Lc), Lloc_c2=str(L2))
    print('%s eta=%s %s v=%s p=%s  can=%.8f  c2=%.8f  margin=%.3e  %s' % (
        tag, e, model, v, p, float(Sc), float(S2), float(Sc-S2), 'OK' if ok else 'FAIL'), flush=True)
    return row

rows = []
# 1) the eight experiment-matched rows
for f in sorted(glob.glob(os.path.join(A,'noise','certs','wd_*.json'))):
    d = json.load(open(f))
    e = F(d['eta']); base = os.path.basename(f)[:-5]
    # wd_etaX_qV-P
    vp = base.split('_q')[1].split('-'); v, p = F(vp[0]), F(vp[1])
    # re-check stored canonical S_up matches our formula
    ch = d['chsh']; c0, ge = loss_functional_F(G_CHSH_F, e)
    beta, _ = check_cert(ch['cert'], ge)
    Lloc = 4*e*(1-e)+2*e*e
    S_stored = c0 + v*(p*beta + (1-p)*Lloc)
    assert S_stored == F(ch['S_up']), (base, S_stored, ch['S_up'])
    rows.append(compare(e, 'wd', v, p, base))

# 2) white-noise grid covering table ranges
for eta, v in [(F(3,4), F(99,100)), (F(4,5), F(99,100)), (F(17,20), F(99,100)),
               (F(9,10), F(99,100)), (F(19,20), F(99,100)),
               (F(3,4), F(1)), (F(4,5), F(1)), (F(9,10), F(1))]:
    rows.append(compare(eta, 'white', v, F(1), 'white'))

# 3) dephasing grid
for eta, p in [(F(3,4), F(96,100)), (F(4,5), F(96,100)), (F(17,20), F(96,100)),
               (F(9,10), F(96,100)), (F(19,20), F(96,100)),
               (F(3,4), F(99,100)), (F(9,10), F(99,100))]:
    rows.append(compare(eta, 'deph', F(1), p, 'deph'))

# strip bulky certs for summary; keep full in separate file
full = rows
summary = [{k: r[k] for k in r if k not in ('cert_can','cert_c2')} for r in full]
OUT=os.path.dirname(os.path.abspath(__file__)); json.dump(full, open(os.path.join(OUT,'noise_class2_certs.json'),'w'))
json.dump(summary, open(os.path.join(OUT,'noise_class2_summary.json'),'w'), indent=1)
print('ALL OK =', all(r['ok'] for r in rows), 'n=', len(rows),
      'min_margin=', min(r['margin'] for r in rows))
