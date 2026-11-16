"""Auto-splitting driver: uniform certificates for one vertex on parameter subintervals, rigorous 2D BnB on each."""
import pickle, sympy as sp, numpy as np, json, sys, time, os
from flint import fmpq
from fractions import Fraction
from uniform_cert import Family, prove
fam, vname, ta, tb, outname = sys.argv[1], sys.argv[2], Fraction(sys.argv[3]), Fraction(sys.argv[4]), sys.argv[5]
marg0 = float(sys.argv[6]) if len(sys.argv) > 6 else 0.002
D = __import__('polygon_io').load('out/polygon_%s.json' % fam); B = D['B']; V = D['R'][vname]['V']; x = B['x']
if fam.startswith('J'): print('regime base', B.get('base'), 'sgn', B.get('sgn'), flush=True)
done = []; todo = [(ta, tb)]
while todo:
    a, b = todo.pop(0); ok = False
    for (k, L) in eval(os.environ.get('KL', '((2, 2), (3, 2))')):
        Fm = Family(x, V['g'], V['d0'], V['d1'], B['c0'], (a + b)/2, (b - a)/2, k, L)
        w = float(b - a); xg = np.linspace(float(a), float(b), 7)
        try:
            eps = Fm.design(xg, marg=marg0, nc=int(os.environ.get('NC', '241')))
        except Exception as e:
            print('design failed', (str(a), str(b)), k, L, e, flush=True); continue
        fc = Fm.float_check(np.linspace(float(a), float(b), 9), n=2001)
        if fc[0] > 1e-12 or fc[2] >= 1: print('float check failed', fc, flush=True); continue
        res = prove(Fm, fmpq(a.numerator, a.denominator), fmpq(b.numerator, b.denominator), NS=32)
        print((str(a), str(b)), k, L, 'eps %.3e' % eps, res, flush=True)
        if res['endpoint_ok'] and res['bnb_ok']:
            done.append(dict(interval=[str(a), str(b)], k=k, L=L, XM=str((a+b)/2), XW=str((b-a)/2), R=[[str(v) for v in r] for r in Fm.R], eps=float(eps), result=res))
            ok = True; break
    if not ok:
        m = (a + b)/2
        if b - a < Fraction(1, 2000): print('giving up at', (str(a), str(b)), flush=True); break
        todo = [(a, m), (m, b)] + todo
    json.dump(dict(family=fam, vertex=vname, target=[str(ta), str(tb)], pieces=done, complete=not todo and ok), open(outname, 'w'), indent=1)
print('covered', [p['interval'] for p in done])
