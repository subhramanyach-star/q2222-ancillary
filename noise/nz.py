"""Noise study: honest devices with source noise, then symmetric loss eta (no-click -> +1).
Noise models on |psi_th> = cos th|00> + sin th|11>:
  white : rho = v |psi><psi| + (1-v) I/4            -> z_ideal -> v z_ideal (marginals and correlators)
  deph  : rho = p |psi><psi| + (1-p) Delta(|psi><psi|) (coherence |00><11| scaled by p; = mixture with the
          classically correlated state cos^2|00><00| + sin^2|11><11|)  -> marginals unchanged, sin(2th) -> p sin(2th)
Fidelity with the target pure state: white F = v + (1-v)/4 ; deph F = 1 - (1-p) sin^2(2th)/2.
Families (each re-optimised over its own state tilt th and 4 x-z angles, at each (eta, noise)):
  chsh   : max CHSH S over th, angles  -> H_min = -log2(1/2+1/2 sqrt(2-S^2/4)) (tight, Pironio et al 2010)
  tilted : constraints {CHSH >=, <A0> =}  (= best alpha<A0>+CHSH, SDP duality), NPA 1+AB float
  cone   : constraints {V_i(theta') >=} for the vertices of C_theta' (= best witness in C_theta'), theta' optimised
Float only (NPA 1+AB upper bound on P_guess up to solver tolerance). Certification: ncert.py.
Usage: nz.py MODEL ETA q1,q2,...  (MODEL wd: white v AND dephasing p, q = v:p)
   (q = v or p, processed in the given order, warm-started along the list)
Output: out/opt_MODEL.jsonl (append, resumable)."""
import warnings; warnings.simplefilter('ignore')
import sys, os, json, time, glob
import os; _H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(_H,'..','eta_study')); sys.path.insert(0,os.path.join(_H,'..','rigorous_rates'))
import numpy as np
from fractions import Fraction as F
from scipy.optimize import minimize
from core import PGm, G_CHSH, Hmin_chsh, lossy
from cone import polygon
def z_ideal(par,model,q):
    th=par[0]; a0,a1,b0,b1=par[1:5]; C=np.cos(2*th); S=np.sin(2*th)
    if model=='white': v,pc=q,1.0
    elif model=='wd': v,pc=q
    else: v,pc=1.0,q
    A=[v*np.cos(a0)*C,v*np.cos(a1)*C]; B=[v*np.cos(b0)*C,v*np.cos(b1)*C]
    E=[v*(np.cos(a)*np.cos(b)+pc*np.sin(a)*np.sin(b)*S) for a in (a0,a1) for b in (b0,b1)]
    return np.array(A+B+E)
def zN(par,eta,model,q): return lossy(z_ideal(par,model,q),eta,1.0)
def fidelity(model,q,th):
    if model=='wd': v,p=q; return v*(1-(1-p)*np.sin(2*th)**2/2)+(1-v)/4
    return q+(1-q)/4 if model=='white' else 1-(1-q)*np.sin(2*th)**2/2
pg_t=PGm(1,1); pg_c=PGm(8,0); E_A0=np.eye(8)[[0]]
_poly={}
def poly(thp):
    k=round(thp,9)
    if k not in _poly:
        if len(_poly)>5000: _poly.clear()
        p,G=polygon(thp); G=list(G)
        while len(G)<8: G.append(G[0])
        _poly[k]=np.array(G)
    return _poly[k]
def obj(fam,x,eta,model,q):
    z=zN(x[:5],eta,model,q)
    if fam=='chsh': return float(G_CHSH@z)
    if fam=='tilted': return pg_t.H(z,Rg=G_CHSH[None,:],Re=E_A0)
    thp=float(np.clip(x[5],0.01*np.pi,np.pi/4-1e-9)); return pg_c.H(z,Rg=poly(thp))
def nm(f,x0s,it):
    best=(-np.inf,None)
    for x0 in x0s:
        x0=np.asarray(x0,float)
        r=minimize(lambda x:-f(x),x0,method='Nelder-Mead',options=dict(maxiter=it,xatol=1e-8,fatol=1e-11,adaptive=True))
        r=minimize(lambda x:-f(x),r.x,method='Nelder-Mead',options=dict(maxiter=it,xatol=1e-9,fatol=1e-12,adaptive=True))
        if -r.fun>best[0]: best=(-r.fun,r.x)
    return best
def prior_certs(eta):
    """devices from the earlier v=1 / v=0.99 certificates at this eta (warm starts)."""
    out={'tilted':[],'cone':[]}
    for fn in glob.glob(os.path.join(_H,'..','eta_study','certs','eta%s_v*.json'%eta)):
        if 'AB' in fn: continue
        d=json.load(open(fn))
        for fam in out:
            if fam in d:
                par=[float(F(s)) for s in d[fam]['par_q']]
                if fam=='cone': par=par+[float(F(d['cone']['thetap']))]
                out[fam].append(par)
    return out
def chsh_starts():
    out=[]
    for t in np.linspace(0.06,0.25,9)*np.pi:
        s=np.sin(2*t); m=np.arccos(1/np.sqrt(1+s*s)); out.append([t,0,np.pi/2,m,-m])
    return out
if __name__=='__main__':
    model=sys.argv[1]; eta=float(sys.argv[2])
    qs=[tuple(float(y) for y in x.split(':')) if model=='wd' else float(x) for x in sys.argv[3].split(',')]
    fams=sys.argv[4].split(',') if len(sys.argv)>4 else ['chsh','tilted','cone']
    OUT='out/opt_%s.jsonl'%model; done={}
    if os.path.exists(OUT):
        for l in open(OUT): d=json.loads(l); done[(d['eta'],tuple(d['q']) if isinstance(d['q'],list) else d['q'],d['fam'])]=d
    pr=prior_certs(eta); prev={}
    for q in qs:
        t0=time.time()
        for fam in fams:
            k=(eta,q,fam)
            if k in done: prev[fam]=done[k]['par']; continue
            if fam=='chsh':
                x0s=chsh_starts()+([prev['chsh']] if 'chsh' in prev else [])
                # cheap screening then polish best 3
                sc=sorted(x0s,key=lambda x:-obj('chsh',np.array(x),eta,model,q))
                Hb,x=nm(lambda x:obj('chsh',x,eta,model,q),sc[:4],3000)
                rec=dict(eta=eta,q=q,model=model,fam=fam,S=float(Hb),H=Hmin_chsh(Hb) if Hb>2 else 0.0,par=list(map(float,x)))
            else:
                x0s=([prev[fam]] if fam in prev else [])+pr[fam][:2]
                if fam=='cone' and 'tilted' in prev: x0s.append(list(prev['tilted'][:5])+[0.24*np.pi])
                if fam=='tilted' and 'chsh' in prev: x0s.append(prev['chsh'])
                if fam=='cone' and 'chsh' in prev: x0s.append(list(prev['chsh'][:5])+[0.245*np.pi])
                Hb,x=nm(lambda x:obj(fam,x,eta,model,q),x0s[:4],600)
                rec=dict(eta=eta,q=q,model=model,fam=fam,H=float(Hb),S=float(G_CHSH@zN(x[:5],eta,model,q)),par=list(map(float,x)))
                if fam=='cone': rec['thetap']=float(np.clip(x[5],0.01*np.pi,np.pi/4-1e-9))
            rec['F']=float(fidelity(model,q,rec['par'][0])); rec['secs']=round(time.time()-t0,1)
            with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
            done[k]=rec; prev[fam]=rec['par']
            print('%s eta=%.4f q=%s %-6s H=%.6f S=%.5f th=%.4fpi F=%.4f'%(model,eta,str(q),fam,rec['H'],rec['S'],rec['par'][0]/np.pi,rec['F']),flush=True)
