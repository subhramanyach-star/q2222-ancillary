import sys, json, numpy as np
sys.argv=['tgen3.py','c_V3s.json']
import tgen3
from flint import arb, fmpq
import mpmath as mp
import g3
mp.mp.dps=50
C=json.load(open('c_V3s.json'))
D=g3.lamv('V3')
def R_full(c):  # for g3: R as function -> emulate via cheb coeff? simpler: direct delta here
    pass
for tv in ['1/100','1/1000','1/64']:
    t=fmpq(*map(int,tv.split('/'))); d=tgen3.hats(arb(t))
    tm=mp.mpf(int(t.p))/int(t.q); dat,c0=g3.data(D,tm); g=dat['g']
    for cv in ['-1','-99/100','0','3/10','98/100','1']:
        c=fmpq(*map(int,cv.split('/'))) if '/' in cv else fmpq(int(cv))
        Gt=tgen3.G(d,arb(c))
        # direct: delta = d0(0)+t*dhat computed by tgen3 vecs
        al=tgen3.pin(d); U,Wh=tgen3.vecs(d,arb(c),al)
        dh=[ (U[1][i]-(d['g'][2+i]/2+d['g'][4+i])) for i in range(2)]
        delta=[mp.mpf(float(tgen3.X0['d0'][i]))+tm*mp.mpf(dh[i].mid().str(40,radius=False)) for i in range(2)]
        cm=mp.mpf(int(c.p))/int(c.q)
        Fd=[f/tm for f in g3.Fs(g,delta,cm)]
        print(tv,cv,[float(x.mid()) for x in Gt],[float(x) for x in Fd])
