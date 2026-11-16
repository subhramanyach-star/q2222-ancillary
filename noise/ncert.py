"""Rigorous two-sided H_min(A1|X=1,E) at (eta, noise) for the devices found by nz.py (pipeline of ../cert_eta.py:
exact-rational NPA duals + 200-bit ball arithmetic).
WINNERS (tilted alpha<A0>+CHSH with dyadic alpha; best witness in C_theta' as an exact convex combination of the
  exact polygon vertices): honest dyadic qubit device under the noise model + loss reaches g.z >= omega (arb);
  Lemma 1 => H_min >= H_lo for EVERY quantum device reaching that score.
CHSH SIDE:
  white : every device whose source is v*(any state) + (1-v)*I/4 (any dimension, any measurements):
          S <= c0 + v*beta (beta = exact NPA certificate of the loss-transformed CHSH)  -- tight (qubit attack).
  deph  : every device whose source is p*(any state) + (1-p)*(any separable state) (contains the dephased model):
          S <= c0 + p*beta + (1-p)*L,  L = local max of the loss-transformed CHSH (= 4eta(1-eta)+2eta^2, exact).
          Valid but NOT tight for the dephased qubit family; the tight comparator is the numerically optimised
          dephased CHSH device (nz.py, float), reported alongside.
  H_up(S_up) via an explicit qubit attack (ball arithmetic) as in cert_eta.chsh_side.
Usage: ncert.py MODEL ETA q1,q2,..   -> certs/MODEL_eta_q.json, out/cert_MODEL.jsonl"""
import warnings; warnings.simplefilter('ignore')
import sys, os, json, time
import os; _H=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(_H,'..','eta_study')); sys.path.insert(0,os.path.join(_H,'..','rigorous_rates'))
import numpy as np
from fractions import Fraction as F
from flint import arb
import cert_eta as CE
from cert_eta import aq, dy, tilted_g, best_alpha, cone_witness, polygon_arb, H_of, ball_up
import rates_core as R
from core import G_CHSH, Hmin_chsh
import nz
QP=None
def honest_model(model):
    def f(par_q,eta,q):
        th=aq(par_q[0]); a=[aq(x) for x in par_q[1:5]]; C=(2*th).cos(); S=(2*th).sin(); eta=aq(eta); q=aq(q)
        if model=='wd': v,pc=aq(QP[0]),aq(QP[1])
        else: v,pc=(q,arb(1)) if model=='white' else (arb(1),q)
        A=[v*a[0].cos()*C,v*a[1].cos()*C]; B=[v*a[2].cos()*C,v*a[3].cos()*C]
        E=[v*(a[x].cos()*a[2+y].cos()+pc*a[x].sin()*a[2+y].sin()*S) for x in (0,1) for y in (0,1)]
        A2=[eta*t+(1-eta) for t in A]; B2=[eta*t+(1-eta) for t in B]
        E2=[eta*eta*E[2*x+y]+eta*(1-eta)*(A[x]+B[y])+(1-eta)**2 for x in (0,1) for y in (0,1)]
        return A2+B2+E2
    return f
def chsh_side(model,eta,q):
    c0,ge=CE.loss_functional_F(CE.G_CHSH_F,eta)
    beta,L,cj,qv=R.certify_hq(ge); assert beta>=0
    Lloc=4*eta*(1-eta)+2*eta*eta; best=None
    for a0 in (1,-1):
        for a1 in (1,-1):
            for b0 in (1,-1):
                for b1 in (1,-1):
                    z=[a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]; val=sum(ge[j]*z[j] for j in range(8))
                    best=val if best is None or val>best else best
    assert best==Lloc, (best,Lloc)
    if model=='wd': v,p=QP; S_up=c0+v*(p*beta+(1-p)*Lloc)
    else: S_up=c0+q*beta if model=='white' else c0+q*beta+(1-q)*Lloc
    att=R.colgen(CE.GUESS,rows_ge=[(G_CHSH.copy(),float(S_up)+1e-11)],iters=120)
    lo,ainfo=R.verify_attack(CE.GUESS,att,rows_ge=[([aq(x) for x in CE.G_CHSH_F],aq(S_up))])
    Hup=H_of(aq(F(lo).limit_denominator(10**15))-arb(10)**-15) if lo else None
    return dict(S_up=str(S_up),S_up_f=float(S_up),level=L,cert=cj,c0=str(c0),beta=str(beta),L_local=str(Lloc),
                H_up=(ball_up(Hup) if lo else None),P_lo=lo,attack=ainfo,H_closed=Hmin_chsh(float(S_up)))
def load(model):
    R_={}
    for l in open('out/opt_%s.jsonl'%model):
        d=json.loads(l); k=(d['eta'],tuple(d['q']) if isinstance(d['q'],list) else d['q'],d['fam'])
        if k not in R_ or d['H']>R_[k]['H']: R_[k]=d
    return R_
def run(model,eta,q,recs):
    global QP
    etaq=F(eta).limit_denominator(10000); t0=time.time()
    if model=='wd': QP=tuple(F(x).limit_denominator(10000) for x in q); qq=QP[0]
    else: qq=F(q).limit_denominator(10000)
    CE.honest_arb=honest_model(model)          # winner() evaluates the honest score with this model
    out=dict(model=model,eta=str(etaq),q=str(qq)); summ=dict(model=model,eta=eta,q=q)
    ch=chsh_side(model,etaq,qq); out['chsh']=ch; summ['S_up']=ch['S_up_f']; summ['H_chsh_up_cert']=ch['H_up']
    if (eta,q,'chsh') in recs: summ['S_chsh_best']=recs[(eta,q,'chsh')]['S']; summ['H_chsh_best']=recs[(eta,q,'chsh')]['H']
    for fam in ('tilted','cone'):
        if (eta,q,fam) not in recs: continue
        par=recs[(eta,q,fam)]['par']; z=nz.zN(par[:5],float(etaq),model,tuple(float(x) for x in QP) if model=='wd' else float(qq))
        try:
            if fam=='tilted':
                a,_=best_alpha(z); a_=dy(a,20); w=CE.winner(par,etaq,qq,tilted_g(a_),'tilted'); w['alpha']=str(a_)
            else:
                thq=dy(min(recs[(eta,q,fam)]['thetap'],np.pi/4-1e-9),30); lamq,_=cone_witness(z,thq)
                p,verts=polygon_arb(thq); g=[sum((aq(l)*V[j] for l,V in zip(lamq,verts)),arb(0)) for j in range(8)]
                w=CE.winner(par,etaq,qq,g,'cone'); w['thetap']=str(thq); w['lam']=[str(x) for x in lamq]
            out[fam]=w; summ[fam]=[w['H_lo'],w['H_up']]; summ[fam+'_float']=recs[(eta,q,fam)]['H']
        except Exception as e:
            summ[fam+'_error']=repr(e)[:200]; summ[fam+'_float']=recs[(eta,q,fam)]['H']
    out['secs']=round(time.time()-t0,1)
    json.dump(out,open('certs/%s_eta%s_q%s.json'%(model,eta,'-'.join(map(str,q)) if model=='wd' else q),'w'))
    with open('out/cert_%s.jsonl'%model,'a') as f: f.write(json.dumps(summ)+'\n')
    print(json.dumps(summ),flush=True)
if __name__=='__main__':
    model=sys.argv[1]; eta=float(sys.argv[2])
    qs=[tuple(float(y) for y in x.split(':')) if model=='wd' else float(x) for x in sys.argv[3].split(',')]
    recs=load(model); done=set()
    if os.path.exists('out/cert_%s.jsonl'%model):
        for l in open('out/cert_%s.jsonl'%model): d=json.loads(l); done.add((d['eta'],tuple(d['q']) if isinstance(d['q'],list) else d['q']))
    for q in qs:
        if (eta,q) in done: continue
        if not all((eta,q,f) in recs for f in ('chsh','tilted','cone')): print('skip (no opt)',eta,q); continue
        run(model,eta,q,recs)
