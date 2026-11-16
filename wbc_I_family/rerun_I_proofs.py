"""Re-run (not redesign) the Arb branch-and-bound for the saved I-family W1 certificates (Proposition 6, 'Certified
normal cone along I_delta').  Run from this directory:  python3 rerun_I_proofs.py   (about 40 s with 4 processes)."""
import sys, os, json, time
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); os.chdir(HERE)
import polygon_io
from fractions import Fraction
from flint import fmpq
from uniform_cert import Family, prove
if __name__ == '__main__':
    D = polygon_io.load('out/polygon_I.json'); B = D['B']; V = D['R']['W1']['V']; x = B['x']
    pieces = []
    d = json.load(open('out/cert_I_W1_low.json')); pieces += [(p['interval'], p['k'], p['L'], p['XM'], p['XW'], p['R']) for p in d['pieces']]
    d = json.load(open('out/cert_I_W1_3_20.json')); pieces.append((d['t_interval'], d['k'], d['L'], d['XM'], d['XW'], d['R']))
    for iv, k, L, XM, XW, R in pieces:
        Fm = Family(x, V['g'], V['d0'], V['d1'], B['c0'], Fraction(XM), Fraction(XW), k, L, R=[[Fraction(v) for v in r] for r in R])
        a, b = Fraction(iv[0]), Fraction(iv[1]); t0 = time.time()
        res = prove(Fm, fmpq(a.numerator, a.denominator), fmpq(b.numerator, b.denominator), NS=32, procs=4)
        print('t in [%s, %s] k=%d L=%d:' % (iv[0], iv[1], k, L), res, '%.1fs' % (time.time() - t0), flush=True)
