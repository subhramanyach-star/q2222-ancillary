"""NPA (real moment matrix) guessing-probability SDPs for the CHSH scenario.
Behaviour coords z=(a0,a1,b0,b1,E00,E01,E10,E11); P(ab|xy)=(1+a alpha_x+b beta_y+ab E_xy)/4.
Observables +-1. Real relaxation (Re of a complex moment matrix is feasible) -> valid upper bound on P_guess.
"""
import itertools, numpy as np, cvxpy as cp

OPS_A=['A0','A1']; OPS_B=['B0','B1']
def reduce(word):
    a=[o for o in word if o[0]=='A']; b=[o for o in word if o[0]=='B']
    def red(s):
        out=[]
        for o in s:
            if out and out[-1]==o: out.pop()
            else: out.append(o)
        return out
    return tuple(red(a)), tuple(red(b))
def canon(word):
    a,b=reduce(word)
    ar,br=tuple(reversed(a)),tuple(reversed(b))
    return min((a,b),(ar,br))
def words(level):
    if level=='1+AB':
        W=[()]+[(o,) for o in OPS_A+OPS_B]+[(x,y) for x in OPS_A for y in OPS_B]
    else:
        L=int(level); W=set()
        for n in range(L+1):
            for w in itertools.product(OPS_A+OPS_B,repeat=n):
                a,b=reduce(w); W.add(a+b)
        W=sorted(W,key=lambda w:(len(w),w))
    return W
def adj(w): return tuple(reversed(w))
class NPA:
    def __init__(s,level):
        s.W=words(level); n=len(s.W)
        keys={}; s.idx=np.zeros((n,n),dtype=int)
        for i,u in enumerate(s.W):
            for j,v in enumerate(s.W):
                k=canon(adj(u)+v)
                if k not in keys: keys[k]=len(keys)
                s.idx[i,j]=keys[k]
        s.keys=keys; s.n=n; s.m=len(keys)
        s.k1=keys[((),())]
        s.kA=[keys[((x,),())] for x in OPS_A]; s.kB=[keys[((),(y,))] for y in OPS_B]
        s.kE=[keys[((x,),(y,))] for x in OPS_A for y in OPS_B]
        # sparse basis matrices
        s.F=[np.zeros((n,n)) for _ in range(s.m)]
        for i in range(n):
            for j in range(n): s.F[s.idx[i,j]][i,j]=1.0
    def block(s):
        y=cp.Variable(s.m)
        G=sum(y[k]*s.F[k] for k in range(s.m))
        return y,G
    def z(s,y):
        return cp.hstack([y[s.kA[0]],y[s.kA[1]],y[s.kB[0]],y[s.kB[1]]]+[y[k] for k in s.kE])

def pguess(level, guess, g=None, omega=None, pobs=None, solver='CLARABEL', verbose=False):
    """guess: ('A',x) local guess of A_x; ('AB',x,y) global guess.
    constraints: g.z >= omega (single functional) and/or z = pobs (full statistics)."""
    N=NPA(level)
    if guess[0]=='A': labels=[(a,) for a in (1,-1)]
    else: labels=[(a,b) for a in (1,-1) for b in (1,-1)]
    ys=[];cons=[];obj=0;ztot=0;norm=0
    for lab in labels:
        y,G=N.block(); ys.append(y)
        cons.append(G>>0)
        zz=N.z(y); ztot=ztot+zz; norm=norm+y[N.k1]
        if guess[0]=='A':
            x=guess[1]; a=lab[0]
            obj=obj+(y[N.k1]+a*y[N.kA[x]])/2
        else:
            x,yy=guess[1],guess[2]; a,b=lab
            obj=obj+(y[N.k1]+a*y[N.kA[x]]+b*y[N.kB[yy]]+a*b*y[N.kE[2*x+yy]])/4
    cons.append(norm==1)
    if g is not None: cons.append(np.asarray(g)@ztot>=omega)
    if pobs is not None: cons.append(ztot==np.asarray(pobs))
    prob=cp.Problem(cp.Maximize(obj),cons)
    return _solve(prob,verbose)

def hmax(level,g):
    """NPA upper bound on h_Q(g)."""
    N=NPA(level); y,G=N.block()
    prob=cp.Problem(cp.Maximize(np.asarray(g)@N.z(y)),[G>>0,y[N.k1]==1])
    prob.solve(solver='CLARABEL',tol_gap_abs=1e-8,tol_gap_rel=1e-8,tol_feas=1e-8)
    return prob.value

def _solve(prob,verbose=False):
    for opts in [dict(solver='CLARABEL',tol_gap_abs=1e-8,tol_gap_rel=1e-8,tol_feas=1e-8,max_iter=400),
                 dict(solver='CLARABEL',tol_gap_abs=1e-7,tol_gap_rel=1e-7,tol_feas=1e-7,max_iter=400),
                 dict(solver='COPT')]:
        try:
            prob.solve(verbose=verbose,**opts)
            if prob.status in ('optimal','optimal_inaccurate'):
                return prob.value, prob.status+'/'+opts['solver']
        except Exception as e:
            pass
    return np.nan,'failed'

def pguess_full_dual(level, guess, pobs):
    """Full-statistics guessing probability; returns (value, status, g*) where g* is the dual functional
    of the constraint sum_e z_e = pobs (so that the single constraint g*.q >= g*.pobs gives the same bound)."""
    N=NPA(level)
    labels=[(a,) for a in (1,-1)] if guess[0]=='A' else [(a,b) for a in (1,-1) for b in (1,-1)]
    cons=[];obj=0;ztot=0;norm=0
    for lab in labels:
        y,G=N.block(); cons.append(G>>0); zz=N.z(y); ztot=ztot+zz; norm=norm+y[N.k1]
        if guess[0]=='A':
            x=guess[1]; obj=obj+(y[N.k1]+lab[0]*y[N.kA[x]])/2
        else:
            x,yy=guess[1],guess[2]; a,b=lab
            obj=obj+(y[N.k1]+a*y[N.kA[x]]+b*y[N.kB[yy]]+a*b*y[N.kE[2*x+yy]])/4
    eq=ztot==np.asarray(pobs)
    cons+= [norm==1, eq]
    prob=cp.Problem(cp.Maximize(obj),cons)
    val,st=_solve(prob)
    return val,st,np.array(eq.dual_value).ravel()
