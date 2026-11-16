import numpy as np, json, sys, warnings
warnings.filterwarnings('ignore')
from npa import pguess
from families import *
LEV='3'
VS=[0.999,0.99,0.98,0.95,0.90]
def H(pg): return -np.log2(pg) if pg==pg else float('nan')
def rate(guess,g,v):
    pg,st=pguess(LEV,guess,g=g,omega=v*(g@P)); return H(pg),st
def full(guess,pv):
    pg,st=pguess(LEV,guess,pobs=pv); return H(pg),st
out={}
rng=np.random.default_rng(1)
which=sys.argv[1]
if which=='PT':
    P=np.array([0,0,0,0,1,1,1,-1])/s2
    gc=g_chsh/(2*s2)
    rho=1-1/s2
    def g_r(r0,r1):
        return gc+np.array([r0/s2+r1/s2, r0/s2-r1/s2, -r0, -r1,0,0,0,0])
    # vertices of octagon: g0 at (rho,0) direction; g1 at angle pi/4 radius? check hL
    fun={'CHSH (centre)':gc,'half-way to g0':g_r(rho/2,0),'vertex g0':g_r(rho,0)}
    for nm,g in fun.items(): print(nm,'g.P',g@P,'hL',hL(g))
    for guess,lab in [(('A',0),'H_min(A0|E)'),(('AB',0,0),'H_min(AB|00,E)')]:
        res={}
        for nm,g in fun.items(): res[nm]=[rate(guess,g,v)[0] for v in VS]
        res['full statistics']=[full(guess,v*P)[0] for v in VS]
        out[lab]=res
elif which in ('th1','th2','th3'):
    th={'th1':np.pi/8,'th2':3*np.pi/16,'th3':27.76/180*np.pi}[which]
    P,order,V=hexagon(th)
    gt=gtilt(th); gch=g_chsh/(g_chsh@P); cen=sum(order)/6
    print('theta/pi',th/np.pi,'CHSH(p)',g_chsh@P,'hL(gch)',hL(gch))
    fun={'tilted CHSH g_theta':gt,'CHSH':gch,'hexagon centroid':cen}
    for i,nm in enumerate(['V1','V2','V1p','V3','V4','V3p']): fun['vertex '+nm]=order[i]
    guess=('A',1)
    res={nm:[rate(guess,g,v)[0] for v in VS] for nm,g in fun.items()}
    res['full statistics']=[full(guess,v*P)[0] for v in VS]
    # search over hexagon: random convex combinations, then local refine at each v
    best={}
    samples=[rng.dirichlet(np.ones(6)*0.7) for _ in range(120)]
    for j,v in enumerate(VS):
        bv=-1;bl=None
        for lam in samples:
            g=sum(l*o for l,o in zip(lam,order)); h,_=rate(guess,g,v)
            if h>bv: bv,bl=h,lam
        # coordinate refinement
        step=0.1
        for it in range(40):
            imp=False
            for k in range(6):
                for sgn in (1,-1):
                    lam=bl.copy(); lam[k]=max(0,lam[k]+sgn*step); lam/=lam.sum()
                    g=sum(l*o for l,o in zip(lam,order)); h,_=rate(guess,g,v)
                    if h>bv+1e-9: bv,bl,imp=h,lam,True
            if not imp: step/=2
            if step<1e-3: break
        gb=sum(l*o for l,o in zip(bl,order))
        best[str(v)]={'H':bv,'lam':bl.tolist(),'g':gb.tolist(),'dist_to_gtilt':float(np.linalg.norm(gb-gt))}
        print(v,'best',bv,'gtilt',res['tilted CHSH g_theta'][j],'full',res['full statistics'][j],flush=True)
    res['best in hexagon (search)']=[best[str(v)]['H'] for v in VS]
    out['H_min(A1|E)']=res; out['best']=best; out['theta']=th
elif which in ('Id1','Id2'):
    d={'Id1':np.pi/12,'Id2':np.pi/6}[which]
    P,g,Y,W1,W2=Idelta(d)
    gch=g_chsh/(g_chsh@P)
    fun={'I_delta (centre)':g,'CHSH':gch,'half-way to W1':g+0.5*(W1[0]*Y[0]+W1[1]*Y[1]),'vertex W1':g+W1[0]*Y[0]+W1[1]*Y[1]}
    print('delta',d,'CHSH(p)',g_chsh@P)
    for guess,lab in [(('AB',0,0),'H_min(AB|00,E)'),(('A',0),'H_min(A0|E)')]:
        res={nm:[rate(guess,gg,v)[0] for v in VS] for nm,gg in fun.items()}
        res['full statistics']=[full(guess,v*P)[0] for v in VS]
        out[lab]=res
    out['delta']=d
out['v']=VS; out['level']=LEV
json.dump(out,open(f'rates_{which}.json','w'),indent=1)
for k,res in out.items():
    if isinstance(res,dict) and k!='best':
        print(k)
        for nm,vals in res.items(): print('  %-28s'%nm,' '.join('%.4f'%x for x in vals))
