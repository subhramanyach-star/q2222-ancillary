"""Consistency of tgen3 exact endpoint Taylor coefficients vs Arb series at t=0 and vs direct values at small t."""
import sys
sys.argv=['tgen3.py', sys.argv[1] if len(sys.argv)>1 else 'c_V3t.json']
import tgen3 as T
from flint import arb, fmpq, arb_series
dS=T.hats(arb_series([arb(0),1],prec=T.MT+1)); q=T.ep_quantities(dS,T.pin(dS),lambda e: arb(e)); bad=0
for key,cs in T.EXC.items():
    x=q[key]; ac=T.coeffs(x,T.MT) if isinstance(x,arb_series) else [x]+[arb(0)]*T.MT
    for j in range(T.MT):
        if not ac[j].contains(arb(cs[j])): bad+=1; print('MISMATCH',key,j,cs[j],ac[j])
print('coefficient mismatches:',bad)
for key in [(1,2,('val',-1,1)),(1,3,('val',-1,1)),(-1,2,('val',-1,1)),(-1,0,'uf')]:
    print(key,[str(c) for c in T.EXC[key]])
for a in [fmpq(0),fmpq(1,2**20),fmpq(1,2**12),fmpq(1,1024)]:
    b=a+fmpq(1,2**30); d=T.hats(T.ball(a,b)); print(float(a),T.endpoint_ok(d,1,a,b),T.endpoint_ok(d,-1,a,b))
