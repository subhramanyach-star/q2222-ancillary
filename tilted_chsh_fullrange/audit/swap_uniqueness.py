"""Audit: in-slice uniqueness for the swapped problem (Alice-overlap slice c'=0 at p' = P p_theta).
With all four mu' > 0, complementary slackness forces z_{a|x} := q_b + a q_{E,x.} = (1 + a q_x) v_{a|x} (Euclidean CS, Gram = I at c'=0);
eliminating q_b leaves q_0 d_0 - q_1 d_1 = const with d_x = (v_{+|x} - v_{-|x})/2, so q is unique iff det(d_0, d_1) != 0.
Exact symbolic det in s, plus sign scan."""
import json, sympy as sp
s,C,r=sp.symbols('s C r',positive=True)
SW=json.load(open('../smalltheta/swap_data.json')); p=[sp.sympify(x,locals={'s':s,'C':C,'r':r}) for x in SW['p']]
v={}
for a in (1,-1):
    for x in (0,1):
        t=1+a*p[x]; w=[p[2]+a*p[4+2*x], p[3]+a*p[5+2*x]]; v[(a,x)]=[wi/t for wi in w]
d=[[(v[(1,x)][i]-v[(-1,x)][i])/2 for i in range(2)] for x in (0,1)]
det=sp.simplify((d[0][0]*d[1][1]-d[0][1]*d[1][0]).subs({C:sp.sqrt(1-s**2),r:sp.sqrt(1+s**2)}))
print('p\' =',p); print('det(d0,d1) =',det)
import mpmath as mp
vals=[float(det.subs(s,sv)) for sv in [mp.mpf(k)/1000 for k in range(1,1001)]]
print('min |det| on s in [0.001,1] grid:',min(abs(x) for x in vals),' signs:',set(1 if x>0 else -1 for x in vals))
