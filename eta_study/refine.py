"""Refinement pass: re-optimise each (eta, family) from several warm starts (CHSH-optimal device, tilted optimum,
neighbouring-eta optima).  Starting from the CHSH-optimal device guarantees H(family) >= H(CHSH) (every family's
constraint set contains the CHSH constraint: alpha=0, resp. the octagon centre at theta'=pi/4).  Writes
out/opt{suf}_v{V}_zrefine.jsonl; load_opt keeps the best record per key."""
import sys, json, numpy as np, warnings; warnings.simplefilter('ignore')
from scipy.optimize import minimize
from core import *
from cone import polygon
V=float(sys.argv[1]); GTAG=sys.argv[2] if len(sys.argv)>2 else 'A1'; FAMS=sys.argv[3].split(',') if len(sys.argv)>3 else ['tilted','cone','full']
GUESS=('AB',int(GTAG[2]),int(GTAG[3])) if GTAG.startswith('AB') else ('A',int(GTAG[1])); SUF='' if GTAG=='A1' else '_'+GTAG
OUT='out/opt%s_v%s_zrefine.jsonl'%(SUF,V)
pg_t=PGm(1,1,guess=GUESS); pg_f=PGm(0,8,guess=GUESS); PC={}
_poly={}
def poly(thp):
    k=round(thp,9)
    if k not in _poly:
        p,G=polygon(thp); G=list(G)
        while len(G)<8: G.append(G[0])
        _poly[k]=np.array(G)
    return _poly[k]
pg_c=PGm(8,0,guess=GUESS)
def obj(fam,par,eta):
    z=z_honest(par[:5],eta,V)
    if fam=='tilted': return pg_t.H(z,Rg=G_CHSH[None,:],Re=np.eye(8)[[0]])
    if fam=='full': return pg_f.H(z,Re=np.eye(8))
    thp=float(np.clip(par[5],0.01*np.pi,np.pi/4)); return pg_c.H(z,Rg=poly(thp))
R=load_opt(V,SUF); etas=sorted(set(e for e,f in R if f=='chsh_opt'))
for fam in FAMS:
    for i,eta in enumerate(etas):
        if (eta,'chsh_opt') not in R or R[(eta,'chsh_opt')]['H']<=1e-9: continue
        cur=R.get((eta,fam)); starts=[]
        c=list(R[(eta,'chsh_opt')]['par'][:5])
        if fam=='cone':
            starts.append(c+[np.pi/4]); starts.append(c+[0.24*np.pi])
            if (eta,'tilted') in R: starts.append(list(R[(eta,'tilted')]['par'][:5])+[0.24*np.pi])
        else:
            starts.append(c)
            if fam=='full' and (eta,'tilted') in R: starts.append(list(R[(eta,'tilted')]['par'][:5]))
        if cur: starts.append(list(cur['par']))
        for j in (i-1,i+1):
            if 0<=j<len(etas) and (etas[j],fam) in R: starts.append(list(R[(etas[j],fam)]['par']))
        best=(cur['H'] if cur else -1,None)
        for x0 in starts:
            if fam=='cone' and len(x0)<6: continue
            f=lambda x:-obj(fam,x,eta)
            r=minimize(f,np.array(x0,float),method='Nelder-Mead',options={'maxiter':500,'xatol':1e-8,'fatol':1e-12,'adaptive':True})
            if -r.fun>best[0]+1e-10: best=(-r.fun,r.x)
        if best[1] is None: print('eta %.3f %s no improvement (%.6g)'%(eta,fam,best[0]),flush=True); continue
        x=best[1]; rec=dict(eta=eta,fam=fam,H=float(best[0]),S=float(G_CHSH@z_honest(x[:5],eta,V)),par=list(map(float,x)),refined=True)
        if fam=='cone': rec['thetap']=float(np.clip(x[5],0.01*np.pi,np.pi/4))
        if fam=='tilted':
            obj(fam,x,eta); dg,de=pg_t.duals(); rec['alpha']=float(-de[0]/dg[0]) if abs(dg[0])>1e-12 else None
        with open(OUT,'a') as fo: fo.write(json.dumps(rec)+'\n')
        R[(eta,fam)]=rec
        print('eta %.3f %-7s H %.6g -> %.6g'%(eta,fam,cur['H'] if cur else float('nan'),rec['H']),flush=True)
