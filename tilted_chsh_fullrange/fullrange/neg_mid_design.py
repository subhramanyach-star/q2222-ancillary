# Negative controls for the mid-pin certificates (V3l strip 8): perturbation k (c+1)(1-c)^6/64 (vanishes at the mid-pinned
# end c=-1, so alpha is unchanged); find k with float F<0 near c=-1, write c_V3l_negmid_*.json.
import sys, json
sys.argv=['vgen.py','c_V3l.json']
exec(open('vgen.py').read().split("def lambda_run")[0])
from numpy.polynomial import chebyshev as Ch, polynomial as Po
from fractions import Fraction
cq=[Fraction(x).limit_denominator(10**9) for x in Ch.poly2cheb(Po.polymul([1,1],Po.polypow([1,-1],6))/64.0)]
J=json.load(open('c_V3l.json')); S=J['strips'][8]; xa0,xb0=q(S['xa']),q(S['xb'])
grid=[-1+fmpq(1,2**e) for e in range(1,40)]+[fmpq(i,2000)-1 for i in range(0,4000,5)]
def scan(R):
    out=[]
    for xi in range(5):
        x=xa0+(xb0-xa0)*fmpq(xi,4); dat=theta_data(arb(x)); P=build(dat,R)
        for c in grid:
            v=[float(z.mid()) for z in F(dat,P,arb(c))]
            if min(v)<0: out.append((float(c),float(x),min(v)))
    return out
made=0
for side in (0,1):
  for k in ['1/1000','1/300','1/100','1/30','1/10','3/10','1','-1/1000','-1/300','-1/100','-1/30','-1/10','-3/10','-1']:
    R=[[q(a),q(b)] for a,b in S['R']]
    for j,c_ in enumerate(cq): R[j][side]+=fmpq(c_.numerator,c_.denominator)*q(k)
    o=scan(R)
    if o:
        lo,hi,mn=min(t[0] for t in o),max(t[0] for t in o),min(t[2] for t in o); print(side,k,len(o),lo,hi,mn,flush=True)
        S2=dict(S); S2['R']=[[str(a),str(b)] for a,b in R]; S2['negative_control']='k=%s comp %d: float F<0 (min %.3g) on c in [%.4f,%.4f]'%(k,side,mn,lo,hi)
        json.dump(dict(vertex=J['vertex'],data=J['data'],pins=J['pins'],strips=[S2]),open('c_V3l_negmid_%d.json'%made,'w')); made+=1
        break
print('made',made)
