#!/usr/bin/env python3
"""Standalone check of the von Neumann attack upper bounds (certs_vn/*.json): exact rational states/observables,
constraint check in ball arithmetic, entropy in ball arithmetic. Needs only python-flint."""
import json, glob, os
from fractions import Fraction as Fr
from flint import arb
from vn_exact import exact_component
from families_arb import build, dot, fr
from verify_all import in_local
res=[]
for fn in sorted(glob.glob('certs_vn/*.json')):
    d=json.load(open(fn)); r=d['rec']; pis=[Fr(x) for x in d['weights']]
    assert all(p>=0 for p in pis) and sum(pis)==1
    data=[exact_component(c,1) for c in d['components']]
    Hb=sum((arb(p.numerator)/p.denominator*H for p,(z,H) in zip(pis,data)),arb(0))
    zb=[arb(sum(p*z[j] for p,(z,H) in zip(pis,data)).numerator)/sum(p*z[j] for p,(z,H) in zip(pis,data)).denominator for j in range(8)]
    pA,gA,orderA=build({'fam':'theta','theta':r['theta'],'g':'tilted'}); V=fr(r['v'])
    if r['kind']=='tilted': ok=bool((dot(gA,zb)-V*dot(gA,pA))>0)
    elif r['kind']=='cone': ok=all(bool((dot(o,zb)-V*dot(o,pA))>0) for o in orderA)
    else:
        T=arb(1)/10**8; w=[(V*pA[j]-(1-T)*zb[j])/T for j in range(8)]; ok=in_local(w); Hb=(1-T)*Hb
    Hhi=float(Hb.mid()+Hb.rad()) if ok else None
    res.append({'id':r['id'],'H_hi':Hhi,'ok':ok}); print('%-28s H(A1|E) <= %s'%(r['id'],'%.6f'%Hhi if ok else 'FAIL'))
json.dump(res,open('out/vn_verified.json','w'),indent=1)
