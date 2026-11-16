import json, mpmath as mp, numpy as np, sympy as sp, sys
sys.path.insert(0,'../fullrange')
from midpin import cands
mp.mp.dps=60
J=json.load(open('../fullrange/thm3_data.json'))
S_,C_,R_=sp.symbols('s C r',positive=True)
def lamv(name):
    V=J['vertices'][name]; return {k:[sp.lambdify((S_,C_,R_),sp.sympify(x,locals={'s':S_,'C':C_,'r':R_}),'mpmath') for x in V[k]] for k in ('g','d0','d1')}
def data(D,t):
    t=mp.mpf(t); s=mp.sqrt(t); v=(s,mp.sqrt(1-t),mp.sqrt(1+t))
    return {k:[f(*v) for f in fs] for k,fs in D.items()}, (1-t)/(1+t)
def cheb(c,n):
    T=[mp.mpf(1),c]
    while len(T)<n+1: T.append(2*c*T[-1]-T[-2])
    return T[:n+1]
def nrm(u,c): return mp.sqrt(max(u[0]**2+u[1]**2+2*c*u[0]*u[1],0))
def make_delta(dat,c0,t,R,mpin=None):
    d0,d1,g=dat['d0'],dat['d1'],dat['g']; n=len(R)-1
    def Q(c): T=cheb(c,n); return [sum(R[j][i]*T[j] for j in range(n+1)) for i in range(2)]
    al=[0,0]
    if mpin is not None:
        e,i,j=mpin; K=cands(g,e); tgt=(K[i]+K[j])/2; Qe=Q(mp.mpf(e))
        Dal=(tgt-(d0[0]+d1[0]*(e-c0))-e*(d0[1]+d1[1]*(e-c0)))/((e-c0)**2)-t*(Qe[0]+e*Qe[1])
        al=[Dal/2,e*Dal/2]
    def delta(c):
        q=Q(c); return [d0[i]+d1[i]*(c-c0)+(c-c0)**2*(t*q[i]+al[i]) for i in range(2)]
    return delta
def Fs(g,dl,c):
    out=[]
    for s0 in (1,-1):
        for s1 in (1,-1):
            u=[g[2+i]/2+dl[i]+s0*g[4+i] for i in range(2)]; w=[g[2+i]/2-dl[i]+s1*g[6+i] for i in range(2)]
            out.append(1-s0*g[0]-s1*g[1]-nrm(u,c)-nrm(w,c))
    return out
if __name__=='__main__':
    name='V3'; D=lamv(name); R=np.load(sys.argv[1]); R=[[mp.mpf(float(a)),mp.mpf(float(b))] for a,b in R]
    mpin=eval(sys.argv[2]) if len(sys.argv)>2 else None
    for t in ['1e-2','1e-3','1e-4','1e-5']:
        t=mp.mpf(t); dat,c0=data(D,t)
        g=dat['g']; K=cands(g,-1); print('t',t,'cands(-1)',[mp.nstr(k,12) for k in K])
        delta=make_delta(dat,c0,t,R,mpin)
        worst=mp.inf
        cgrid=[-1+mp.mpf(10)**(-e) for e in np.arange(12,0,-0.5)]+[mp.mpf(x) for x in np.linspace(-0.9,0.9,37)]+[1-mp.mpf(10)**(-e) for e in np.arange(1,12,0.5)]+[mp.mpf(-1),mp.mpf(1),c0]
        rows=[]
        for c in sorted(cgrid):
            G=[f/t for f in Fs(g,delta(c),c)]; m=min(G); rows.append((c,m,G))
            worst=min(worst,m)
        print('  worst G=%s'%mp.nstr(worst,6))
        for c,m,G in rows:
            if c<-0.99 or m<0: print('   1+c=%.3e  G=%s'%(float(1+c),[mp.nstr(x,5) for x in G]))
