# Negative controls for the pinned_term path (Wa strip 22, near pi/4). Perturbation dP(c) = k (c-1)(1+c)^6/64 added to the
# Chebyshev part (vanishes at the pinned end c=1, so the pin correction alpha is unchanged); report where float F < 0.
import sys, json, numpy as np
sys.argv=['vgen.py','c_Wa.json']
exec(open('vgen.py').read().split("def lambda_run")[0])
from numpy.polynomial import chebyshev as Ch, polynomial as Po
pc = Po.polymul([-1,1], Po.polypow([1,1],6))/64.0
cc = Ch.poly2cheb(pc)                      # 8 Chebyshev coefficients (degree 7)
from fractions import Fraction
cq = [Fraction(x).limit_denominator(10**9) for x in cc]
J=json.load(open('c_Wa.json')); S=J['strips'][22]; xa0,xb0=q(S['xa']),q(S['xb'])
grid=[1-fmpq(1,2**e) for e in range(1,40)]+[fmpq(i,200)-1 for i in range(400)]
def scan(R):
    out=[]
    for xi in range(5):
        x=xa0+(xb0-xa0)*fmpq(xi,4); dat=theta_data(arb(x)); P=build(dat,R)
        for c in grid:
            v=[float(z.mid()) for z in F(dat,P,arb(c))]
            if min(v)<0: out.append((float(c),float(x),min(v),v.index(min(v))))
    return out
for side in ():
  for k in ['-3','-1','-3/10','-1/10','1/10','3/10','1','3']:
    R=[[q(a),q(b)] for a,b in S['R']]
    for j,c_ in enumerate(cq): R[j][side]+=fmpq(c_.numerator,c_.denominator)*q(k)
    o=scan(R)
    print(side,k,len(o), (min(t[0] for t in o),max(t[0] for t in o), min(t[2] for t in o), set(t[3] for t in o)) if o else '', flush=True)
print('--- small k; writing c_Wa_negpin_*.json')
for side,k,tag in [(0,'1/100','a'),(0,'-1/1000','b'),(1,'1/100000','c')]:
    R=[[q(a),q(b)] for a,b in S['R']]
    for j,c_ in enumerate(cq): R[j][side]+=fmpq(c_.numerator,c_.denominator)*q(k)
    o=scan(R)
    print(side,k,len(o), (min(t[0] for t in o),max(t[0] for t in o), min(t[2] for t in o), set(t[3] for t in o)) if o else '', flush=True)
    S2=dict(S); S2['R']=[[str(a),str(b)] for a,b in R]; S2['negative_control']='R_1.. perturbed by %s*(c-1)(1+c)^6/64 on side %d'%(k,side)
    json.dump(dict(vertex='Wa',data='oct',pins=J['pins'],strips=[S2]),open('c_Wa_negpin_%s.json'%tag,'w'))
