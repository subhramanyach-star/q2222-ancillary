import sys, json
from fractions import Fraction as Fr
base=json.load(open('c_V3t.json'))
sys.argv=['tgen3.py','c_V3t.json']
import tgen3 as T
from flint import arb, fmpq
for name,dv in [('a',(Fr(1,256),0)),('b',(-Fr(1,256),0)),('c',(0,Fr(1,256))),('d',(0,-Fr(1,256)))]:
    c=json.loads(json.dumps(base)); c['Rt'][0]=[str(Fr(c['Rt'][0][i])+dv[i]) for i in range(2)]
    json.dump(c,open('c_V3t_negE%s.json'%name,'w'))
    T.RT=[[T.Q(a),T.Q(b)] for a,b in c['Rt']]
    # rebuild monomials
    import importlib
    T.RMONO=[]
    for i in range(2):
        m=[fmpq(0)]*(len(T.RT)+1)
        for j in range(len(T.RT)):
            for k,v in enumerate(T._TM[j]): m[k]+=T.RT[j][i]*v; m[k+1]+=T.RT[j][i]*v
        m[0]+=T.R0[i]; T.RMONO.append(m)
    T._SH.clear()
    for tv in [fmpq(1,2000),fmpq(1,300)]:
        d=T.hats(arb(tv)); al=T.pin(d); neg=[]
        for i in range(1,401):
            cc=fmpq(-1)+fmpq(i,200); g=T.G(d,arb(cc),al)
            for k in range(4):
                if g[k]<0: neg.append((float(cc),k))
        print(name,float(tv),len(neg),neg[:2],neg[-1:] )
