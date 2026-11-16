import sys, pickle, json, sympy as sp, mpmath as mp
import os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
mp.mp.dps=60
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); import polygon_io; D=polygon_io.load(os.path.join(os.path.dirname(os.path.abspath(__file__)),'out','polygon_I.json')); B=D['B']; V=D['R']['W1']['V']; t=B['x']
print('tight', V.get('tight'))
for tv in [sp.Rational(1,40),sp.Rational(1,10),sp.Rational(1,5),sp.Rational(26,100)]:
    g=[mp.mpf(sp.N(e.subs(t,tv),70)) for e in V['g']]; d0=[mp.mpf(sp.N(e.subs(t,tv),70)) for e in V['d0']]; d1=[mp.mpf(sp.N(e.subs(t,tv),70)) for e in V['d1']]; c0=mp.mpf(sp.N(B['c0'].subs(t,tv),70))
    def Phi(c):
        dl=[d0[i]+d1[i]*(c-c0) for i in range(2)]; tot=0
        for x in (0,1):
            Bx=[g[2+i]/2+(dl[i] if x==0 else -dl[i]) for i in range(2)]; e=[g[4+2*x],g[5+2*x]]
            tot+=max(a*g[x]+mp.sqrt((Bx[0]+a*e[0])**2+(Bx[1]+a*e[1])**2+2*c*(Bx[0]+a*e[0])*(Bx[1]+a*e[1])) for a in (1,-1))
        return tot
    # F_{s0s1} individually
    def Fs(c):
        dl=[d0[i]+d1[i]*(c-c0) for i in range(2)]; r=[]
        for s0 in (1,-1):
            for s1 in (1,-1):
                u=[g[2+i]/2+dl[i]+s0*g[4+i] for i in range(2)]; v=[g[2+i]/2-dl[i]+s1*g[6+i] for i in range(2)]
                A=u[0]**2+u[1]**2+2*c*u[0]*u[1]; Bq=v[0]**2+v[1]**2+2*c*v[0]*v[1]
                r.append(1-s0*g[0]-s1*g[1]-mp.sqrt(A)-mp.sqrt(Bq))
        return r
    p=[0,0,0,0]; 
    print('t=%s: Phi(c0)-1=%s  F(c0)=%s  F\'(c0)=%s'%(tv, mp.nstr(Phi(c0)-1,5),[mp.nstr(f,3) for f in Fs(c0)],[mp.nstr(mp.diff(lambda c: Fs(c)[k],c0),3) for k in range(4)]))
