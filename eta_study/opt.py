"""Optimise honest devices (state tilt theta, four x-z-plane angles) at each detector efficiency eta for:
  chsh_maxent : CHSH score only, |Phi+>                  -> H_min = f(S), S = 2sqrt2 eta^2 + 2(1-eta)^2 (exact)
  chsh_opt    : CHSH score only, theta + settings optimised (maximise S; H_min(S) is increasing)
  tilted      : best member of the tilted-CHSH family alpha<A0> + CHSH (alpha optimised) = constraints {CHSH >=, <A0> =}
  cone        : best witness in the Thm-2 cone C_theta' (theta' optimised) = constraints {V_i >= } for all vertices
  full        : full statistics (upper reference; = best arbitrary single witness for H_min)
Entropy: H_min(A1|X=1,E) from NPA level 1+AB (valid lower bound; certified later).  Resumable: out/opt_v{v}.jsonl."""
import warnings; warnings.simplefilter('ignore')
import numpy as np, json, os, sys, time, warnings
warnings.filterwarnings('ignore')
from scipy.optimize import minimize
from core import *
from cone import polygon
V=float(sys.argv[1]) if len(sys.argv)>1 else 1.0
FAMS=sys.argv[2].split(',') if len(sys.argv)>2 else ['chsh','tilted','full','cone']
GTAG=sys.argv[3] if len(sys.argv)>3 else 'A1'
GUESS=('AB',int(GTAG[2]),int(GTAG[3])) if GTAG.startswith('AB') else ('A',int(GTAG[1]))
SUF='' if GTAG=='A1' else '_'+GTAG
OUT='out/opt%s_v%s_%s.jsonl'%(SUF,V,'+'.join(FAMS))
import glob
ETAS=[float(x) for x in sys.argv[4].split(',')] if len(sys.argv)>4 else sorted(set([0.995,0.67,0.68,0.69,0.70,0.71,0.72,0.73,0.74,0.75,0.76,0.77,0.78,0.79,0.80,0.81,0.82,0.825,0.828,0.83,0.835,0.84,0.85,0.86,0.87,0.88,0.89,0.90,0.92,0.94,0.96,0.98,0.99]))
done={}
for fn in glob.glob('out/opt%s_v%s*.jsonl'%(SUF,V)):
    for l in open(fn): d=json.loads(l); done[(d['eta'],d['fam'])]=d
pg_t=PGm(1,1,guess=GUESS); pg_c=PGm(8,0,guess=GUESS); pg_f=PGm(0,8,guess=GUESS); pg_s=PGm(1,0,guess=GUESS)
E_A0=np.eye(8)[[0]]; I8=np.eye(8)
def S_of(par,eta): return G_CHSH@z_honest(par,eta,V)
_poly={}
def poly(thp):
    k=round(thp,9)
    if k not in _poly:
        p,G=polygon(thp); G=list(G)
        while len(G)<8: G.append(G[0])
        _poly[k]=np.array(G)
    return _poly[k]
def obj(fam,par,eta):
    z=z_honest(par,eta,V)
    if fam=='tilted': return pg_t.H(z,Rg=G_CHSH[None,:],Re=E_A0)
    if fam=='full': return pg_f.H(z,Re=I8)
    if fam=='cone':
        thp=float(np.clip(par[5],0.01*np.pi,np.pi/4)); return pg_c.H(z,Rg=poly(thp))
def run_nm(f,x0s,maxiter):
    best=(-np.inf,None)
    for x0 in x0s:
        r=minimize(lambda x:-f(x),x0,method='Nelder-Mead',options={'maxiter':maxiter,'xatol':1e-7,'fatol':1e-10,'adaptive':True})
        r=minimize(lambda x:-f(x),r.x,method='Nelder-Mead',options={'maxiter':maxiter,'xatol':1e-8,'fatol':1e-11,'adaptive':True})
        if -r.fun>best[0]: best=(-r.fun,r.x)
    return best
rng=np.random.default_rng(7)
prev={}
def std_starts(eta):
    th=[np.pi/4,np.pi/8,0.15*np.pi,0.08*np.pi]
    out=[]
    for t in th:
        s=np.sin(2*t); r=np.sqrt(1+s*s); m=np.arccos(1/r)
        out.append(np.array([t,0,np.pi/2,m,-m]))
    return out
def save(rec):
    with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
    done[(rec['eta'],rec['fam'])]=rec
    print('v=%g eta=%.3f %-11s H=%.6f S=%.5f th=%.4fpi %s'%(V,rec['eta'],rec['fam'],rec['H'],rec['S'],rec['par'][0]/np.pi,rec.get('extra','')),flush=True)
for eta in sorted(ETAS,reverse=True):
    t0=time.time()
    # CHSH families
    if 'chsh' in FAMS and (eta,'chsh_maxent') not in done:
        S=V*2*np.sqrt(2)*eta**2+2*(1-eta)**2 if V==1.0 else S_of([np.pi/4,0,np.pi/2,np.pi/4,-np.pi/4],eta)
        HH=Hmin_chsh(S) if GTAG=='A1' else pg_s.H(z_honest([np.pi/4,0,np.pi/2,np.pi/4,-np.pi/4],eta,V),Rg=G_CHSH[None,:])
        save(dict(eta=eta,fam='chsh_maxent',H=HH,Hvn=Hvn_chsh(S),S=float(S),par=[np.pi/4,0,np.pi/2,np.pi/4,-np.pi/4]))
    if 'chsh' in FAMS and (eta,'chsh_opt') not in done:
        x0s=std_starts(eta)+([np.array(prev['chsh_opt'])] if 'chsh_opt' in prev else [])
        Sb,x=run_nm(lambda p:S_of(p,eta),x0s,3000)
        HH=Hmin_chsh(Sb) if GTAG=='A1' else pg_s.H(z_honest(x,eta,V),Rg=G_CHSH[None,:])
        save(dict(eta=eta,fam='chsh_opt',H=HH,Hvn=Hvn_chsh(Sb),S=float(Sb),par=list(map(float,x))))
    if (eta,'chsh_opt') not in done:
        Sb,x=run_nm(lambda p:S_of(p,eta),std_starts(eta)+([np.array(prev['chsh_opt'])] if 'chsh_opt' in prev else []),3000)
        done[(eta,'chsh_opt')]=dict(par=list(map(float,x)))
    prev['chsh_opt']=done[(eta,'chsh_opt')]['par']
    for fam in [f for f in ('tilted','full','cone') if f in FAMS]:
        if (eta,fam) in done: prev[fam]=done[(eta,fam)]['par']; continue
        x0s=[np.array(prev[fam])] if fam in prev else []
        base=std_starts(eta)+[np.array(done[(eta,'chsh_opt')]['par'])]
        if fam=='cone':
            base=[np.concatenate([b,[t]]) for b in base for t in (b[0] if b[0]>0.03 else 0.1,)]
            if (eta,'tilted') in done: tp=done[(eta,'tilted')]['par']; base.append(np.array(tp+[max(tp[0],0.03)]))
        else:
            if fam=='full' and (eta,'tilted') in done: base.append(np.array(done[(eta,'tilted')]['par']))
        x0s=x0s+base
        Hb,x=run_nm(lambda p:obj(fam,p,eta),x0s[:3],450)
        rec=dict(eta=eta,fam=fam,H=float(Hb),S=float(S_of(x[:5],eta)),par=list(map(float,x)))
        if fam=='tilted':
            obj(fam,x,eta); dg,de=pg_t.duals(); rec['alpha']=float(-de[0]/dg[0]) if abs(dg[0])>1e-12 else None; rec['extra']='alpha=%s'%rec['alpha']
        if fam=='cone': rec['thetap']=float(np.clip(x[5],0.01*np.pi,np.pi/4)); rec['extra']="theta'=%.4fpi"%(rec['thetap']/np.pi)
        save(rec); prev[fam]=rec['par']
    print('  eta %.3f done in %.0fs'%(eta,time.time()-t0),flush=True)
