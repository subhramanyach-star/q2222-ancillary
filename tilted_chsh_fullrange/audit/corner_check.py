"""Audit: exact check that u^_{-}(c=-1, t=0) = 0 (P = M = 0 at the corner), required by tgen3.edge_dc's Hadamard branch
(ta = 0), which is asserted only in its docstring.  Also checks Taylor coefficients of P, M at the corner."""
import sys; sys.path.insert(0,'../smalltheta'); sys.argv=['tgen3.py','c_V3t.json']
import os; os.chdir('../smalltheta')
import tgen3 as T, flint
from flint import fmpq, fmpq_series
for cert in ['c_V3t.json','c_V3m.json']:
    import json
    C=json.load(open(cert)); T.R0=[T.Q(a) for a in C['R0']]; T.RT=[[T.Q(a),T.Q(b)] for a,b in C['Rt']]
    # rebuild RMONO for this cert
    TM=T._cheb_mono(len(T.RT)-1); RM=[]
    for i in range(2):
        m=[fmpq(0)]*(len(T.RT)+1)
        for j in range(len(T.RT)):
            for k,v in enumerate(TM[j]): m[k]+=T.RT[j][i]*v; m[k+1]+=T.RT[j][i]*v
        m[0]+=T.R0[i]; RM.append(m)
    T.RMONO=RM; T._SH.clear()
    flint.ctx.cap=8
    d=T.hats_exact(); Re=T.Rpoly(fmpq(-1)); al=T.pin_generic(d,Re)
    U,_=T.vecs(d,fmpq(-1),al)
    v=U[-1]; P=v[0]+v[1]; M=v[0]-v[1]
    cs=lambda x: [str(c) for c in (x.coeffs() if isinstance(x,fmpq_series) else [x])][:4]
    print(cert,'P(-1,t) coeffs',cs(P),' M(-1,t) coeffs',cs(M))
    flint.ctx.cap=64
