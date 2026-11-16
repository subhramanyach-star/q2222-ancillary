"""Rigorous upper bounds on the quantum maximum h_Q(c) of a 2222 Bell expression c.z (+constant) by
exact-rational NPA (sum-of-squares) dual certificates.
Certificate: symmetric rational Z >= 0 (verified exactly) with  sum_{(i,j) in class k} Z_ij = -c_k (k = z-moments),
= 0 for every other non-identity moment class; then for EVERY quantum realisation (any dimension)
  beta - c.z = <Z, Gamma> >= 0,  beta = trace Z  (identity class = diagonal),
where Gamma is the real part of the moment matrix (PSD, entries = moments).  So h_Q(c) <= beta."""
import numpy as np, cvxpy as cp, json
from fractions import Fraction
from flint import fmpq
from npa_struct import Struct
from exact_psd import is_psd
_S={}
def S(level):
    if level not in _S: _S[level]=Struct(level)
    return _S[level]
def _fq(x): x=Fraction(x); return fmpq(x.numerator,x.denominator)
def solve_dual(c,level):
    s=S(level); n=s.n
    Z=cp.Variable((n,n),symmetric=True)
    tgt=[0.0]*s.m
    for t,k in enumerate(s.kz): tgt[k]=-float(c[t])
    cons=[Z>>0]
    for k in range(s.m):
        if k==s.k1: continue
        cons.append(sum(Z[i,j] for (i,j) in s.cls[k])==tgt[k])
    prob=cp.Problem(cp.Minimize(cp.trace(Z)),cons)
    try:
        prob.solve(solver='CLARABEL',tol_gap_abs=1e-10,tol_gap_rel=1e-10,tol_feas=1e-10,max_iter=500)
    except Exception:
        prob.solve(solver='CLARABEL')
    return Z.value, prob.value
def round_cert(Zn,c,level,eps,den=2**40):
    s=S(level); n=s.n
    w,V=np.linalg.eigh((Zn+Zn.T)/2); w=np.maximum(w,0)+eps
    Zp=(V*w)@V.T
    Q=[[None]*n for _ in range(n)]
    for i in range(n):
        for j in range(i,n):
            q=Fraction(round(Zp[i,j]*den),den); Q[i][j]=q; Q[j][i]=q
    tgt={k:Fraction(0) for k in range(s.m)}
    for t,k in enumerate(s.kz): tgt[k]=-Fraction(c[t])
    for k in range(s.m):
        if k==s.k1: continue
        cl=s.cls[k]; cur=sum(Q[i][j] for (i,j) in cl); r=(tgt[k]-cur)/len(cl)
        for (i,j) in cl: Q[i][j]+=r   # classes are closed under transpose, so symmetry is kept
    beta=sum(Q[i][i] for i in range(n))
    return Q,beta
def certify(c,level='3',eps_list=(1e-9,1e-8,1e-7,1e-6,1e-5)):
    """c: 8 rationals (Fraction).  Returns dict with exact beta (Fraction) and Z, or None."""
    Zn,val=solve_dual(c,level)
    if Zn is None: return None
    for eps in eps_list:
        Q,beta=round_cert(Zn,c,level,eps)
        ok,info=is_psd([[_fq(x) for x in row] for row in Q])
        if ok: return {'level':level,'c':[str(Fraction(x)) for x in c],'beta':beta,'Z':Q,'solver':val,'eps':eps}
    return None
def verify(cd):
    """Standalone exact check of a certificate dict (c, level, Z as strings or Fractions). Returns beta or raises."""
    level=cd['level']; s=Struct(level); n=s.n
    c=[Fraction(x) for x in cd['c']]; Z=[[Fraction(x) for x in row] for row in cd['Z']]
    assert len(Z)==n and all(len(r)==n for r in Z)
    for i in range(n):
        for j in range(n): assert Z[i][j]==Z[j][i]
    tgt={k:Fraction(0) for k in range(s.m)}
    for t,k in enumerate(s.kz): tgt[k]=-c[t]
    for k in range(s.m):
        if k==s.k1: continue
        assert sum(Z[i][j] for (i,j) in s.cls[k])==tgt[k], 'linear constraint %d'%k
    ok,info=is_psd([[_fq(x) for x in row] for row in Z]); assert ok, info
    return sum(Z[i][i] for i in range(n))
def to_json(cd):
    return {'level':cd['level'],'c':cd['c'],'beta':str(cd['beta']),'Z':[[str(x) for x in row] for row in cd['Z']]}
