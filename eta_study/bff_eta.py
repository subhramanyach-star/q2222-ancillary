"""von Neumann H(A1|X=1,E) under detector efficiency eta: dual-certified BFF lower bounds (../rates_v11/bffcert.py,
exact-rational residual certification) for the tilted / cone witnesses and full statistics at the H_min-optimised
honest devices (alpha re-scanned for vN), CHSH by the tight closed form 1-h(1/2+1/2 sqrt(S^2/4-1)) evaluated at the
certified S_up(eta) (an upper bound for every CHSH protocol) and at the honest S_max (achievable).
EAT finite-size rates (Liu et al. 2021 Thm 3 + spot checking, ../rates_v11/eat.py) from the certified affine
min-tradeoff functions.  Resumable: out/bff_v{v}.jsonl."""
import sys, os, json, time, glob, warnings
warnings.simplefilter('ignore')
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rates_v11'))
import numpy as np
from bffexact import RateExact as Rate
import eat
from core import z_honest, G_CHSH, Hvn_chsh
V=float(sys.argv[1]); ETAS=[float(x) for x in sys.argv[2].split(',')]; M=int(sys.argv[3]) if len(sys.argv)>3 else 8
OUT='out/bff_v%s.jsonl'%V
done=set()
if os.path.exists(OUT):
    for l in open(OUT): d=json.loads(l); done.add((d['eta'],d['fam'],d['m']))
from core import load_opt
allr=load_opt(V)
certs={}
if os.path.exists('out/cert_results.jsonl'):
    for l in open('out/cert_results.jsonl'): d=json.loads(l); certs[(d['eta'],d['v'])]=d
RT={M:Rate(m=M,xstar=1)}; R6=Rate(m=6,xstar=1)
def save(rec):
    with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
    print('eta=%.3f %-7s m=%d H_cert=%.5f (num %.5f) EAT1e8=%.4f EAT1e10=%.4f %s'%(rec['eta'],rec['fam'],rec['m'],rec['H_cert'],rec['H_num'],rec['eat'][0],rec['eat'][1],rec.get('alpha','')),flush=True)
def eats(K0,kap,z):
    return [eat.eat_rate(K0,np.asarray(kap),z,n)['rate'] for n in (1e8,1e10)]
for eta in ETAS:
    # CHSH
    if (eta,'chsh',0) not in done and (eta,'chsh_opt') in allr:
        c=allr[(eta,'chsh_opt')]; z=z_honest(c['par'],eta,V); S=float(G_CHSH@z)
        Sup=certs.get((eta,V),{}).get('S_up')
        if S>2.0+1e-9:
            K0,kap=eat.chsh_tangent(S); e=eats(K0,kap,z)
        else: e=[0.0,0.0]
        save(dict(eta=eta,fam='chsh',m=0,H_cert=Hvn_chsh(S),H_num=Hvn_chsh(S),H_upper_all_chsh=(Hvn_chsh(Sup) if Sup else None),S=S,eat=e))
    for fam in ('tilted','cone','full'):
        if (eta,fam,M) in done or (eta,fam) not in allr: continue
        r=allr[(eta,fam)]; z=z_honest(r['par'],eta,V); t0=time.time(); extra={}
        if fam=='full':
            b=RT[M].bound('full',z=z)
        else:
            if fam=='tilted':
                a0=r.get('alpha') or 0.0; best=None
                for a in sorted(set([0.0,a0*0.5,a0,a0*1.5,a0*2.0,a0*3.0])):
                    g=np.array([a,0,0,0,1,1,1,-1.]); bb=R6.bound('score',g=g,omega=g@z-1e-10)
                    if best is None or bb['H_cert']>best[1]: best=(a,bb['H_cert'])
                g=np.array([best[0],0,0,0,1,1,1,-1.]); extra['alpha']=best[0]
            else:
                cr=json.load(open('certs/eta%s_v%s.json'%(eta,V))) if os.path.exists('certs/eta%s_v%s.json'%(eta,V)) else None
                if cr is None or 'cone' not in cr: continue
                from fractions import Fraction; g=np.array([float(Fraction(x)) for x in cr['cone']['upper']['gtilde']])
            b=RT[M].bound('score',g=g,omega=g@z-1e-10); extra['g']=list(map(float,g))
        if b['status']!='ok': print('fail',eta,fam,b['status']); continue
        rec=dict(eta=eta,fam=fam,m=M,H_cert=float(b['H_cert']),H_num=float(b['H_num']),K0=float(b['K0']),kappa=list(map(float,b['kappa'])),
                 z=list(map(float,z)),secs=round(time.time()-t0,1),**extra)
        rec['eat']=eats(rec['K0'],rec['kappa'],z) if rec['H_cert']>0 else [0.0,0.0]
        save(rec)
