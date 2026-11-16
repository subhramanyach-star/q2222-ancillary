"""Runner: two-sided rigorous H_min bounds for every number in prl_v13 (checkpointed, resumable).
Writes out/results.jsonl (one record per case) and certs/<id>.json (all data needed by verify_all.py)."""
import sys, json, os, time, warnings
warnings.filterwarnings('ignore')
import numpy as np
from fractions import Fraction
from flint import arb
import rates_core as R, qubit
from families_arb import build, dot, fr
OUT='out/results.jsonl'
done=set()
if os.path.exists(OUT):
    for l in open(OUT): done.add(json.loads(l)['id'])
def Fq(x): return Fraction(x)
def H(P): return float(-np.log2(P))
def seeds(spec):
    if spec['fam']=='theta':
        from families_arb import theta_of
        return [qubit.target_strategy(float(theta_of(spec['theta']).mid()))]
    if spec['fam']=='Id':
        ph={'pi/6':np.pi/6,'pi/12':np.pi/12}[spec['delta']]
        return [qubit.phiplus_strategy([0,np.pi/2+ph,np.pi/2,-ph])]
    if spec['fam']=='PT':
        return [qubit.phiplus_strategy([0,np.pi/2,np.pi/4,-np.pi/4])]
    return []
def run_single(cid,spec,guess,v):
    p,g,_=build(spec); V=fr(v); om=V*dot(g,p)
    att=R.colgen(guess,rows_ge=[(R.mid(g),float(om.mid())+1e-11)],init=seeds(spec),iters=120)
    lo,ainfo=R.verify_attack(guess,att,rows_ge=[(g,om)])
    mu=Fraction(max(0.0,float(-att['dual_ub'][0]))).limit_denominator(10**8)
    up,uinfo=R.upper_single(guess,g,om,mu)
    up=min(up,1.0)
    return lo,up,{'attack':ainfo,'upper':uinfo}
def run_full(cid,spec,guess,v,T=Fraction(1,10**8)):
    p,_,_=build({k:x for k,x in spec.items() if k!='g'}); V=fr(v); pv=[V*x for x in p]
    tgt=[float(x.mid())/(1-float(T)) for x in pv]
    att=R.colgen(guess,rows_eq=True,rhs_eq=tgt,init=seeds(spec),iters=200,tol=1e-11)
    lo,ainfo=R.verify_attack(guess,att,p_eq=pv,t=T)
    up,uinfo=R.upper_full(guess,att['dual_lin'],pv)
    return lo,min(up,1.0),{'attack':ainfo,'upper':uinfo}
def run_cone(cid,spec,guess,v):
    p,_,order=build(spec); V=fr(v)
    rows=[(g,V*dot(g,p)) for g in order]
    att=R.colgen(guess,rows_ge=[(R.mid(g),float(om.mid())+1e-11) for g,om in rows],init=seeds(spec),iters=150)
    lo,ainfo=R.verify_attack(guess,att,rows_ge=rows)
    m=np.maximum(-np.asarray(att['dual_ub']),0); mu=float(m.sum())
    if mu<=1e-12:
        return lo,1.0,{'attack':ainfo,'upper':{'trivial':True},'gstar':None}
    lam=[Fraction(x/mu).limit_denominator(10**9) for x in m]; s=sum(lam); lam=[x/s for x in lam]
    gs='lam:'+','.join(str(x) for x in lam); sp=dict(spec); sp['g']=gs
    p,g,_=build(sp); om=V*dot(g,p)
    up,uinfo=R.upper_single(guess,g,om,Fraction(mu).limit_denominator(10**8))
    return lo,min(up,1.0),{'attack':ainfo,'upper':uinfo,'gstar':gs}
def main(caselist):
    for c in caselist:
        if c['id'] in done: continue
        t=time.time(); guess=tuple(c['guess'])
        if c['kind']=='single': lo,up,info=run_single(c['id'],c['spec'],guess,c['v'])
        elif c['kind']=='full': lo,up,info=run_full(c['id'],c['spec'],guess,c['v'])
        elif c['kind']=='cone': lo,up,info=run_cone(c['id'],c['spec'],guess,c['v'])
        rec={'id':c['id'],'case':c,'P_lo':lo,'P_up':up,'H_lo':H(up),'H_up':(H(lo) if lo else None),
             'levels':[e.get('level') for e in info['upper'].get('per_e',[])],'secs':round(time.time()-t,1)}
        if c['kind']=='cone': rec['gstar']=info.get('gstar')
        json.dump({'case':c,'P_lo':lo,'P_up':up,**info},open('certs/%s.json'%c['id'],'w'))
        with open(OUT,'a') as f: f.write(json.dumps(rec)+'\n')
        print(c['id'],'H in [%.6f, %s]'%(rec['H_lo'],'%.6f'%rec['H_up'] if rec['H_up'] is not None else 'None'),rec['levels'],rec['secs'],flush=True)
if __name__=='__main__':
    main(json.load(open(sys.argv[1])))
