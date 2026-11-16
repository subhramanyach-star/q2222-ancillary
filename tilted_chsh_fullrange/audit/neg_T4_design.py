"""Audit negative control for tgen3's T4/T2 path (c in (c0,1), t small): R -> R + k(1-c)v, i.e. R0 += 2kv, Rt_0 -= kv
(delta(+1) and the endpoint check unchanged).  Float probe of G on a c-grid at a few t; writes c_V3t_negT4_*.json."""
import sys, json, os
os.chdir('../smalltheta'); sys.argv=['tgen3.py','c_V3t.json']; sys.path.insert(0,'.')
import tgen3 as T
from flint import arb, fmpq
base=json.load(open('c_V3t.json'))
def load(C):
    T.R0=[T.Q(a) for a in C['R0']]; T.RT=[[T.Q(a),T.Q(b)] for a,b in C['Rt']]
    TM=T._cheb_mono(len(T.RT)-1); RM=[]
    for i in range(2):
        m=[fmpq(0)]*(len(T.RT)+1)
        for j in range(len(T.RT)):
            for k,v in enumerate(TM[j]): m[k]+=T.RT[j][i]*v; m[k+1]+=T.RT[j][i]*v
        m[0]+=T.R0[i]; RM.append(m)
    T.RMONO=RM; T._SH.clear()
found=[]
for k in ['16']:
  for v in [(1,0),(-1,0),(0,1),(0,-1)]:
    C=json.loads(json.dumps(base)); kq=fmpq(*map(int,k.split('/')))
    C['R0']=[str(T.Q(C['R0'][i])+2*kq*v[i]) for i in range(2)]
    C['Rt'][0]=[str(T.Q(C['Rt'][0][i])-kq*v[i]) for i in range(2)]
    load(C); worst=(1e9,None)
    for tv in [fmpq(1,4096),fmpq(1,2048),fmpq(1,1100)]:
        d=T.hats(arb(tv)); c0=float((1-tv)/(1+tv))
        for j in range(1,200):
            c=fmpq(int(c0*2**30),2**30)+(1-fmpq(int(c0*2**30),2**30))*fmpq(j,200)
            G=T.G(d,arb(c)); m=min(float(x.mid()) for x in G)
            if m<worst[0]: worst=(m,(float(tv),float(c)))
        Gc1=min(float(x.mid()) for x in T.G(d,arb(1)))
    print(k,v,'min G on (c0,1):',worst,' G(1) min:',Gc1)
    if worst[0]<-1e-9 and Gc1>-1e-12: found.append((k,v,C,worst))
for n,(k,v,C,w) in enumerate(found[:2]):
    C['tstrips']=[base['tstrips'][0]]; C['note']='AUDIT negative control (T4/T2 region): R += %s (1-c) %s; float min G %.3g at (t,c)=%s'%(k,v,w[0],w[1])
    json.dump(C,open('../audit/c_V3t_negT4_%d.json'%n,'w'),indent=0); print('wrote',n,C['note'])
