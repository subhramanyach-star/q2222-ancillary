import sys, json
from fractions import Fraction as Fr
import numpy as np
base=json.load(open('c_V3t.json'))
sys.argv=['tgen3.py','c_V3t.json']
import tgen3 as T
from flint import arb, fmpq
# Chebyshev coefficients of (1-c)^6 (exact)
from numpy.polynomial import chebyshev as Ch, polynomial as Po
mono=[Fr(1),Fr(-6),Fr(15),Fr(-20),Fr(15),Fr(-6),Fr(1)]
import sympy as sp
c_=sp.symbols('c'); poly=sp.expand((1-c_)**6)
a=[sp.Rational(0)]*7
rem=poly
for j in range(6,-1,-1):
    lc=sp.Poly(rem,c_).coeff_monomial(c_**j); Tj=sp.chebyshevt(j,c_); lcT=sp.Poly(Tj,c_).coeff_monomial(c_**j)
    a[j]=lc/lcT; rem=sp.expand(rem-a[j]*Tj)
assert rem==0
def mk(name,eps,v):
    c=json.loads(json.dumps(base))
    c['Rt']=[[str(Fr(c['Rt'][j][i])+Fr(str(a[j]))*eps*v[i]) for i in range(2)] for j in range(7)]
    json.dump(c,open('c_V3t_negF%s.json'%name,'w')); return c
for name,eps,v in [('a',Fr(1,1024),(1,0)),('b',Fr(-1,1024),(1,0)),('c',Fr(1,1024),(0,1)),('d',Fr(-1,1024),(0,1)),('e',Fr(1,1024),(1,1)),('f',Fr(-1,1024),(1,1)),('g',Fr(1,1024),(1,-1)),('h',Fr(-1,1024),(1,-1))]:
    c=mk(name,eps,v)
    T.RT=[[T.Q(x),T.Q(y)] for x,y in c['Rt']]; T.RMONO=[]
    for i in range(2):
        m=[fmpq(0)]*(len(T.RT)+1)
        for j in range(len(T.RT)):
            for k,vv in enumerate(T._TM[j]): m[k]+=T.RT[j][i]*vv; m[k+1]+=T.RT[j][i]*vv
        m[0]+=T.R0[i]; T.RMONO.append(m)
    T._SH.clear()
    for tv in [fmpq(1,2000),fmpq(1,300)]:
        d=T.hats(arb(tv)); al=T.pin(d); neg=[]
        for i in range(1,401):
            cc=fmpq(-1)+fmpq(i,200)
            if i<20:
                for sub in range(1,10): 
                    c2=fmpq(-1)+fmpq(i-1,200)+fmpq(sub,2000); g=T.G(d,arb(c2),al)
                    neg+=[(float(c2),k) for k in range(4) if g[k]<0]
            g=T.G(d,arb(cc),al)
            neg+=[(float(cc),k) for k in range(4) if g[k]<0]
        print(name,float(tv),len(neg),neg[:2],neg[-1:] )
