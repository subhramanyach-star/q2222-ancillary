"""Independent cross-check of sgen.G against direct F/t^k from the exact swapped data (mpmath, 60 digits)."""
import sys, json, mpmath as mp, sympy as sp
sys.argv=['sgen.py','c_V1s.json']
import sgen
from flint import arb, fmpq
mp.mp.dps=60
SW=json.load(open('swap_data.json'))['vertices']['V1']; s,C,r=sp.symbols('s C r',positive=True)
def val(e,tv):
    sv=mp.sqrt(tv); return mp.mpf(sp.N(sp.sympify(e,locals={'s':s,'C':C,'r':r}).subs({s:sp.sqrt(sp.Rational(tv)),C:sp.sqrt(1-sp.Rational(tv)),r:sp.sqrt(1+sp.Rational(tv))}),70))
R=[[mp.mpf(sp.Rational(x)) for x in sgen.CERT['R'][i]] for i in range(2)]
def nrm(v,c): return mp.sqrt((1+c)/2*(v[0]+v[1])**2+(1-c)/2*(v[0]-v[1])**2)
worst=0
for tv in ['1/1000','1/100','1/63']:
    t=mp.mpf(sp.Rational(tv)); g=[val(e,tv) for e in SW['g']]; d0=[val(e,tv) for e in SW['d0']]; d1=[val(e,tv) for e in SW['d1']]
    dq=sgen.tdata(sgen.Q(tv),sgen.Q(tv))
    for cv in ['-1','-7/10','-1/10','1/5','3/5','99/100']:
        c=mp.mpf(sp.Rational(cv))
        Rc=[sum(R[i][k]*c**k for k in range(len(R[i]))) for i in range(2)]
        dl=[d0[i]+d1[i]*c+t*c*c*Rc[i] for i in range(2)]
        Gs=sgen.Gall(dq,arb(sgen.Q(cv)))
        for k,(s0,s1) in enumerate(sgen.COMBOS):
            u=[g[2+i]/2+dl[i]+s0*g[4+i] for i in range(2)]; w=[g[2+i]/2-dl[i]+s1*g[6+i] for i in range(2)]
            F=1-s0*g[0]-s1*g[1]-nrm(u,c)-nrm(w,c); Fd=F/t**(2 if k==0 else 1)
            gv=mp.mpf(Gs[k].mid().str(50,radius=False)); worst=max(worst,abs(gv-Fd))
            assert Gs[k].contains(arb(str(mp.nstr(Fd,40)))) or abs(gv-Fd)<1e-25, (tv,cv,k,gv,Fd)
print('max |G_sgen - F/t^k| =', mp.nstr(worst,5))
