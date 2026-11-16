#!/usr/bin/env python3
"""Stand-alone re-check of the CONE side of the every-eta theorem (prl_v16, App. sec:lossy) at all certified efficiencies:
runs the unchanged per-point checks of vpt.py (exact re-derivation of q, c0, active cones, exhaustive vertex enumeration of
C(q), forced certificate data, exact first-order/endpoint checks, Arb bisection over c in [-1,1]) on every data/t4/<eta>.json.
USAGE: python3 vpt_t4.py [datadir]"""
import sys, os, glob, json, time
D = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', 't4')
VPT = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'vpt.py')
src = open(VPT).read().split("if __name__ == '__main__':"); head, main = src[0], src[1]
sys.argv = ['vpt.py']; exec(head)
P = {}; CERT = {}
for f in sorted(glob.glob(os.path.join(D, '*.json'))):
    e = os.path.basename(f)[:-5]; st = json.load(open(f)); P[e] = st['point']
    for k, c in st['certs'].items(): CERT['%s/%s' % (e, k)] = c
sys.argv = ['vpt.py'] + sorted(P)
exec(main.replace('\n    ', '\n', ).lstrip() if False else 'if True:' + main)
