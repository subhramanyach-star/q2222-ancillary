#!/usr/bin/env python3
"""Stand-alone check of Proposition S2(iv) of prl_v19: under white noise / dephasing / experiment-matched
(white+dephasing) source noise, class-2 CHSH variants have certified S_up below the canonical one.
Inputs: noise_class2_certs.json. Expected: ALL OK = True, min margin >= 1.96e-2."""
import os, sys, json, itertools
from fractions import Fraction as F
H=os.path.dirname(os.path.abspath(__file__)); A=os.path.join(H,'..')
sys.path[:0]=[os.path.join(A,'rigorous_rates'), os.path.join(A,'eta_study')]
from verify_all import check_cert
def loss_fun(sv, eta):
    sa=[sv[0]+sv[1],sv[2]+sv[3]]; sb=[sv[0]+sv[2],sv[1]+sv[3]]
    c0=(1-eta)**2*sum(sv)
    g=[eta*(1-eta)*sa[0],eta*(1-eta)*sa[1],eta*(1-eta)*sb[0],eta*(1-eta)*sb[1]]+[eta*eta*x for x in sv]
    return c0,g
def local_max(g):
    best=None
    for a0,a1,b0,b1 in itertools.product((1,-1),repeat=4):
        z=[a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1]; val=sum(g[j]*z[j] for j in range(8))
        best=val if best is None or val>best else best
    return best
rows=json.load(open(os.path.join(H,'noise_class2_certs.json'))); ok_all=True; worst=None
for r in rows:
    e,v,p=F(r['eta']),F(r['v']),F(r['p']); model=r['model']
    for label,sv,cjkey in [('can',[1,1,1,-1],'cert_can'),('c2',[-1,-1,-1,1],'cert_c2')]:
        c0,g=loss_fun(sv,e); beta,_=check_cert(r[cjkey],g); L=local_max(g)
        if model=='white': S=c0+v*beta
        elif model=='deph': S=c0+p*beta+(1-p)*L
        else: S=c0+v*(p*beta+(1-p)*L)
        if label=='can': Sc=S
        else: S2=S
    ok=S2<=Sc; ok_all&=ok
    if worst is None or Sc-S2<worst[0]: worst=(Sc-S2,r['tag'])
    print('%s %s  can=%.8f c2=%.8f %s'%(r['tag'],model,float(Sc),float(S2),'OK' if ok else 'FAIL'))
print('min margin %.3e at %s'% (float(worst[0]),worst[1])); print('ALL OK =',ok_all, 'n=',len(rows))
