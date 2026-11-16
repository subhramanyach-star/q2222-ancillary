"""von Neumann H(A1|X=1,E) under detector efficiency eta, devices optimised for vN.
Search: BFF (m=4, one solver trial; float-certified) as objective, Nelder-Mead over (theta, 4 angles[, alpha]).
Final: exact-PSD BFF certificate (bffexact.RateExact, m=8) at the best device -> certified lower bound and affine
min-tradeoff function -> EAT rates (Liu et al. 2021 Thm 3 + spot checking, ../rates_v11/eat.py) for n=1e8, 1e10.
CHSH: tight closed form at the honest S_max (achievable) and at the certified S_up (upper bound for every CHSH protocol).
Resumable: out/vn_v{V}.jsonl."""
import sys, os, json, time, warnings, io, contextlib
warnings.simplefilter('ignore')
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rates_v11'))
import numpy as np
from scipy.optimize import minimize
from bffcert import Rate
from bffexact import RateExact
import eat
from core import *
V=float(sys.argv[1]); ETAS=[float(x) for x in sys.argv[2].split(',')]; FAMS=sys.argv[3].split(',') if len(sys.argv)>3 else ['chsh','tilted','full']
OUT='out/vn_v%s.jsonl'%V
done=set()
if os.path.exists(OUT):
    for l in open(OUT): d=json.loads(l); done.add((d['eta'],d['fam']))
R=load_opt(V)
certs={}
if os.path.exists('out/cert_results.jsonl'):
    for l in open('out/cert_results.jsonl'): d=json.loads(l); certs[(d['eta'],d['v'],tuple(d.get('guess',['A',1])))]=d
r6=Rate(m=4); r6.TRIALS=r6.TRIALS[:1]; r6.HARD_NODES=0
RX=RateExact(m=8)
def quiet(f,*a,**k):
    return f(*a,**k)
def eats(K0,kap,z): return [float(eat.eat_rate(K0,np.asarray(kap),z,n)['rate']) for n in (1e8,1e10)]
def save(rec):
    with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
    print('eta=%.3f %-6s H_cert=%.6f EAT(1e8,1e10)=(%.4f,%.4f) %s'%(rec['eta'],rec['fam'],rec['H_cert'],rec['eat'][0],rec['eat'][1],rec.get('note','')),flush=True)
for eta in ETAS:
    if (eta,'chsh_opt') not in R: continue
    zc=z_honest(R[(eta,'chsh_opt')]['par'],eta,V); Sc=float(G_CHSH@zc)
    if 'chsh' in FAMS and (eta,'chsh') not in done:
        Sup=certs.get((eta,V,('A',1)),{}).get('S_up')
        e=eats(*eat.chsh_tangent(Sc),zc) if Sc>2+1e-9 else [0.0,0.0]
        save(dict(eta=eta,fam='chsh',H_cert=Hvn_chsh(Sc),H_upper_every_chsh=(Hvn_chsh(Sup) if Sup else None),S=Sc,S_up=Sup,eat=e,par=R[(eta,'chsh_opt')]['par']))
    if Sc<=2+1e-9: continue
    for fam in [f for f in ('tilted','full') if f in FAMS]:
        if (eta,fam) in done: continue
        t0=time.time()
        if fam=='tilted':
            def f(x):
                z=z_honest(x[:5],eta,V); g=np.array([x[5],0,0,0,1,1,1,-1.]); b=r6.bound('score',g=g,omega=g@z-1e-10)
                return -b['H_cert'] if np.isfinite(b['H_cert']) else 1.0
            starts=[list(R[(eta,'chsh_opt')]['par'])+[0.0]]
            if (eta,'tilted') in R: starts.append(list(R[(eta,'tilted')]['par'][:5])+[R[(eta,'tilted')].get('alpha') or 0.0])
        else:
            def f(x):
                z=z_honest(x[:5],eta,V); b=r6.bound('full',z=z); return -b['H_cert'] if np.isfinite(b['H_cert']) else 1.0
            starts=[list(R[(eta,'chsh_opt')]['par'])]
            for k in ((eta,'tilted'),(eta,'full')):
                if k in R: starts.append(list(R[k]['par'][:5]))
        best=(np.inf,None)
        for x0 in starts:
            r=minimize(f,np.array(x0,float),method='Nelder-Mead',options={'maxiter':150,'xatol':1e-6,'fatol':1e-7,'adaptive':True})
            if r.fun<best[0]: best=(r.fun,r.x)
        x=best[1]; z=z_honest(x[:5],eta,V)
        if fam=='tilted':
            g=np.array([x[5],0,0,0,1,1,1,-1.]); b=RX.bound('score',g=g,omega=g@z-1e-10)
        else: b=RX.bound('full',z=z)
        rec=dict(eta=eta,fam=fam,H_cert=float(b['H_cert']),H_num=float(b['H_num']),H_search_m4=float(-best[0]),K0=float(b['K0']),
                 kappa=list(map(float,b['kappa'])),par=list(map(float,x[:5])),z=list(map(float,z)),S=float(G_CHSH@z),secs=round(time.time()-t0,1))
        if fam=='tilted': rec['alpha']=float(x[5]); rec['g']=list(map(float,g)); rec['omega']=float(g@z-1e-10); rec['note']='alpha=%.4f'%x[5]
        rec['eat']=eats(rec['K0'],rec['kappa'],z) if rec['H_cert']>0 else [0.0,0.0]
        save(rec)
