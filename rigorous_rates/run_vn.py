"""Rigorous UPPER bounds on the von Neumann rate H(A1|X=1,E) (explicit mixed-state qubit attacks, exact rational data,
ball-arithmetic entropies), to pair with the certified BFF LOWER bounds of ../rates_v11.  Checkpointed:
out/vn_results.jsonl, certs_vn/<id>.json."""
import json, os, sys, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
from fractions import Fraction as Fr
from flint import arb
import families as F
import vn_attack as V
os.makedirs('certs_vn',exist_ok=True)
OUT='out/vn_results.jsonl'; done=set()
if os.path.exists(OUT): done={json.loads(l)['id'] for l in open(OUT)}
TH={'pi/8':np.pi/8,'27.76deg':27.76/180*np.pi}
def seeds(th,v,rng,n=12):
    s=[V.noisy_target_x(th,v),V.noisy_target_x(th,min(1,v+0.002)),V.noisy_target_x(th,1.0)]
    for _ in range(n): s.append(V.noisy_target_x(th,v)+np.concatenate([rng.normal(0,.2,4),rng.normal(0,.2,6),rng.normal(0,1.0,4)]))
    return s
def verify(att,gx,rows_exact=None,p_eq=None,t=None):
    """rows_exact: list of (g Fractions..., omega float-free) -> we check with Fractions vs arb omega via families_arb."""
    pis=[Fr(float(p)).limit_denominator(10**14) for p in att['pi']]; s=sum(pis); pis=[p/s for p in pis]
    comps=[V.freeze(c) for c in att['cols']]; data=[V.exact_component(fz,gx) for fz in comps]
    Hb=sum((arb(p.numerator)/p.denominator*H for p,(z,H) in zip(pis,data)),arb(0))
    zbar=[sum(p*z[j] for p,(z,H) in zip(pis,data)) for j in range(8)]
    return pis,comps,Hb,zbar
def main():
    from families_arb import build, dot, fr
    import rates_core as R
    rng=np.random.default_rng(7)
    cases=[]
    for th in ['pi/8','27.76deg']:
        for v in ['0.999','0.99','0.98','0.95']:
            for kind in ['tilted','cone','full']: cases.append((th,v,kind))
    for th,v,kind in cases:
        cid=f"vn_{th.replace('/','_')}_{kind}_v{v}"
        if cid in done: continue
        t0=time.time(); thf=TH[th]; vf=float(v); P,order,_=F.hexagon(thf); gt=F.gtilt(thf)
        if kind=='tilted': rows=[(gt,vf*(gt@P)+1e-10)]; rhs=None
        elif kind=='cone': rows=[(o,vf*(o@P)+1e-10) for o in order]; rhs=None
        else: rows=[]; T=Fr(1,10**8); rhs=list(vf*P/(1-float(T)))
        att=V.colgen_vn(1,rows,seeds(thf,vf,rng),iters=40,starts=12,rng=rng,rhs_eq=rhs)
        pis,comps,Hb,zbar=verify(att,1)
        # exact constraint check in ball arithmetic against the exact target
        pA,gA,orderA=build({'fam':'theta','theta':th,'g':'tilted'}); Vv=fr(v)
        zb=[arb(x.numerator)/x.denominator for x in zbar]
        if kind=='tilted': ok=bool((dot(gA,zb)-Vv*dot(gA,pA))>0)
        elif kind=='cone': ok=all(bool((dot(o,zb)-Vv*dot(o,pA))>0) for o in orderA)
        else:
            TT=arb(1)/10**8; w=[(Vv*pA[j]-(1-TT)*zb[j])/TT for j in range(8)]; ok=R.in_local(w); Hb=(1-TT)*Hb
        Hhi=float(Hb.mid()+Hb.rad()) if ok else None
        rec={'id':cid,'theta':th,'v':v,'kind':kind,'H_hi':Hhi,'lp':att['value'],'ncols':len(comps),'secs':round(time.time()-t0,1),'last':att['hist'][-1]}
        json.dump({'rec':rec,'weights':[str(p) for p in pis],'components':comps},open(f'certs_vn/{cid}.json','w'))
        with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
        print(rec,flush=True)
if __name__=='__main__': main()
