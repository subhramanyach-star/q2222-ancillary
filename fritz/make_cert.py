"""Build exact rational NPA-1+AB certificate that the Fritz image of tilted CHSH (a=14/17) is not quantum.
Functional (original sign convention of the image): w.z = -(a0+a1) - b0 - (E00+E01+E10-E11).
Output: cert_q14.json (all rationals as strings)."""
import numpy as np, cvxpy as cp, json
from fractions import Fraction as Fr
from npa import monomials, canon
mons=monomials('1+AB'); n=len(mons)
names=[(('A0',),()),(('A1',),()),((),('B0',)),((),('B1',))]+[((a,),(b,)) for a in ['A0','A1'] for b in ['B0','B1']]
known=[canon(k[0]+k[1]) for k in names]
cls={}
for i in range(n):
    for j in range(n):
        cls.setdefault(canon(tuple(reversed(mons[i]))+mons[j]),[]).append((i,j))
w=[Fr(-1),Fr(-1),Fr(-1),Fr(0),Fr(-1),Fr(-1),Fr(-1),Fr(1)]
T=Fr(61,20)
def target(k):
    if k==((),()): return T
    if k in known: return -w[known.index(k)]
    return Fr(0)
Y=cp.Variable((n,n),symmetric=True); lam=cp.Variable()
cons=[Y-lam*np.eye(n)>>0]+[sum(Y[i,j] for i,j in ent)==float(target(k)) for k,ent in cls.items()]
cp.Problem(cp.Maximize(lam),cons).solve(solver='CLARABEL'); print('min eig margin',lam.value)
Yr=[[Fr(round(Y.value[i,j]*10**6),10**6) for j in range(n)] for i in range(n)]
for i in range(n):
    for j in range(i): Yr[i][j]=Yr[j][i]
for k,ent in cls.items():
    res=target(k)-sum(Yr[i][j] for i,j in ent)
    for i,j in ent: Yr[i][j]+=res/len(ent)
json.dump(dict(monomials=[list(m) for m in mons],w=[str(x) for x in w],T=str(T),
   Y=[[str(x) for x in r] for r in Yr],a='14/17',tilted_bound='52/17'),open('cert_q14.json','w'),indent=0)
print('written')
