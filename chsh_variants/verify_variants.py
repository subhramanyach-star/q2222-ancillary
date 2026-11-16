#!/usr/bin/env python3
"""Stand-alone check (exact rationals only, no SDP solver) of Proposition S2(iii) of prl_v18:
at each of the 112 efficiencies of Theorem S6, the largest recorded value of the second class of CHSH variants
(sign pattern (-,-,-,+), missed detections recorded as +1) is below the certified canonical S_up(eta).
Inputs: class2_certs.json (exact NPA certificates, produced by class2_cert.py) and ../lossycone/data/t4/<eta>.json
(canonical CHSH certificate). Expected output: 112 lines 'OK', min margin 3.200e-03 at eta = 0.980, 'ALL OK = True'."""
import os, sys, json
from fractions import Fraction as F
H = os.path.dirname(os.path.abspath(__file__)); A = os.path.join(H, '..')
sys.path[:0] = [os.path.join(A, 'rigorous_rates'), os.path.join(A, 'eta_study')]
from verify_all import check_cert
from cert_eta import loss_functional_F, G_CHSH_F
def class2(eta):
    sv = [-1, -1, -1, 1]; sa = [sv[0]+sv[1], sv[2]+sv[3]]; sb = [sv[0]+sv[2], sv[1]+sv[3]]
    c0 = (1-eta)**2*sum(sv)
    return c0, [eta*(1-eta)*sa[0], eta*(1-eta)*sa[1], eta*(1-eta)*sb[0], eta*(1-eta)*sb[1]] + [eta**2*x for x in sv]
rows = json.load(open(os.path.join(H, 'class2_certs.json'))); ok_all = True; worst = None
assert len(rows) == 112
for r in rows:
    e = F(r['eta']); st = json.load(open(os.path.join(A, 'lossycone', 'data', 't4', r['eta'] + '.json')))
    c0, ge = loss_functional_F(G_CHSH_F, e); b1, _ = check_cert(st['chsh_cert'], ge); S1 = c0 + b1
    c02, c2 = class2(e); assert [F(x) for x in r['cert']['c']] == c2, 'certificate is for a different functional'
    b2, _ = check_cert(r['cert'], c2); S2 = c02 + b2; ok = S2 <= S1; ok_all &= ok
    if worst is None or S1 - S2 < worst[0]: worst = (S1 - S2, r['eta'])
    print('%s  canonical S_up = %.8f  class-2 S_up = %.8f  %s' % (r['eta'], float(S1), float(S2), 'OK' if ok else 'FAIL'))
print('min margin %.3e at eta = %s' % (float(worst[0]), worst[1])); print('ALL OK =', ok_all)
