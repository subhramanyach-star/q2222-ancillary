#!/usr/bin/env python3
"""Standalone exact verifier (Python stdlib only) for:
  Proposition. For two parties with two binary observables each, marginals included, Fritz's Question 14
  (arXiv:1202.0141, eq. (23)) has a negative answer: F(x)F maps a Tsirelson inequality outside the quantum cone.
Steps (all in exact rational arithmetic):
 1. g = tilted-CHSH Tsirelson inequality  52/17 - (a<A0> + <A0B0>+<A0B1>+<A1B0>-<A1B1>) >= 0, a=14/17.
    Validity on Q: quantum max of a A0 + CHSH is sqrt(8+2a^2) [Acin-Massar-Pironio PRL 108,100402 (2012)];
    here 8+2a^2 = (52/17)^2 exactly.
 2. Image y = F g F^T with Fritz's F^{st} = [[1,0,1],[0,1,0],[1,0,-1]], index order (-1,0,+1) = (A0, none, A1).
 3. Separating functional w; exact identity  T - w.z = <Y, Gamma(z)>  for every NPA level-1+AB moment matrix,
    with Y positive definite (exact LDL^T). Hence w.z <= T on Q, but w.z* = 41/13 > T.
 4. Sanity: a=0 maps CHSH to the Tsirelson point (Fritz Prop. 15)."""
import json, sys
from fractions import Fraction as Fr
D=json.load(open(sys.argv[1] if len(sys.argv)>1 else 'cert_q14.json'))
a=Fr(D['a']); h=Fr(D['tilted_bound'])
assert h>0 and h*h==8+2*a*a, 'tilted bound not exact'
F=[[1,0,1],[0,1,0],[1,0,-1]]
def tensor(beta,h):   # beta=(a0,a1,b0,b1,E00,E01,E10,E11); index 0<->A0/B0, 1<->none, 2<->A1/B1
    a0,a1,b0,b1,E00,E01,E10,E11=beta
    return [[-E00,-a0,-E01],[-b0,h,-b1],[-E10,-a1,-E11]]
def image(g):
    y=[[sum(F[s][t]*g[t][u]*F[r][u] for t in range(3) for u in range(3)) for r in range(3)] for s in range(3)]
    n=y[1][1]; return [y[0][1]/n,y[2][1]/n,y[1][0]/n,y[1][2]/n,y[0][0]/n,y[0][2]/n,y[2][0]/n,y[2][2]/n], n
z,norm=image(tensor([a,0,0,0,1,1,1,-1],h))
print('image normalisation y^{00} =',norm,'>0')
print('image z* =',[str(v) for v in z])
assert z==[Fr(-7,26),Fr(-7,26),0,0,Fr(-17,26),Fr(-17,26),Fr(-17,26),Fr(17,26)]
# no-signalling check of the image (it must be NS by Fritz Thm 4, since g is also a Bell inequality: local max 2+a >= h)
assert h>=2+a   # bound >= local max 2+a, so g is also a valid Bell inequality
P=[(1+s*z[x]+t*z[2+y]+s*t*z[4+2*x+y])/4 for x in (0,1) for y in (0,1) for s in (1,-1) for t in (1,-1)]
assert min(P)>=0; print('image is a valid NS box, min prob =',min(P))
# sanity a=0
# a=0: image correlators are -(2/h)*[[1,1],[1,-1]] with zero marginals; with h=2*sqrt(2) this is the
# Tsirelson point up to flipping Bob's outcomes (a symmetry of Q).  Check the rational part exactly with h=1:
z0,_=image(tensor([0,0,0,0,1,1,1,-1],Fr(1)))
assert z0==[0,0,0,0,-2,-2,-2,2]; print('sanity a=0: image = -(2/h)(1,1,1,-1) correlators, zero marginals -> Tsirelson point at h=2sqrt2')
# NPA certificate
def red(seq):
    out=[]
    for s in seq:
        if out and out[-1]==s: out.pop()
        else: out.append(s)
    return tuple(out)
def canon(word):
    A=red([x for x in word if x[0]=='A']); B=red([x for x in word if x[0]=='B'])
    return min((A,B),(tuple(reversed(A)),tuple(reversed(B))))
mons=[tuple(m) for m in D['monomials']]; n=len(mons)
Y=[[Fr(x) for x in r] for r in D['Y']]; w=[Fr(x) for x in D['w']]; T=Fr(D['T'])
assert all(Y[i][j]==Y[j][i] for i in range(n) for j in range(n))
names=[(('A0',),()),(('A1',),()),((),('B0',)),((),('B1',))]+[((p,),(q,)) for p in ('A0','A1') for q in ('B0','B1')]
known=[canon(k[0]+k[1]) for k in names]
sums={}
for i in range(n):
    for j in range(n):
        k=canon(tuple(reversed(mons[i]))+mons[j]); sums[k]=sums.get(k,0)+Y[i][j]
for k,v in sums.items():
    want = T if k==((),()) else (-w[known.index(k)] if k in known else 0)
    assert v==want, (k,v,want)
print('identity  T - w.z = <Y,Gamma>  holds exactly for all', len(sums),'moment classes')
# exact LDL^T
M=[r[:] for r in Y]; piv=[]
for k in range(n):
    d=M[k][k]; assert d>0, 'not PD'; piv.append(d)
    for i in range(k+1,n):
        f=M[i][k]/d
        for j in range(k+1,n): M[i][j]-=f*M[k][j]
print('Y positive definite: all', n, 'LDL pivots > 0; smallest =', float(min(piv)))
val=sum(wi*zi for wi,zi in zip(w,z))
print('w.z* =',val,'=',float(val),'  >  T =',T,'=',float(T),'  gap =',val-T)
assert val>T
print('VERIFIED: Fritz image of a valid Tsirelson inequality is not a quantum behaviour.')
