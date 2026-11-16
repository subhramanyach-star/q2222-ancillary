#!/usr/bin/env python3
"""STANDALONE VERIFIER (no SDP/LP solver, no floating-point optimisation).
Needs only python3 + python-flint (pip install python-flint).  For every certs/<id>.json it re-derives, in exact
rational and ball (arb, 200-bit) arithmetic:
  lower bound on H_min : Lemma 1 / Lemma 2 of PROOFS.md with exact-rational NPA certificates (PSD by exact LDL)
  upper bound on H_min : the stored explicit qubit attack (mixture weights, dyadic angles and state) reproduces
                         the constraint (score >= omega, all six cone vertices, or p_v exactly via a local remainder)
Prints the verified interval [H_lo, H_hi] for each case and writes out/verified.json.
Usage: python3 verify_all.py [id-substring ...]"""
import sys, os, json, glob, itertools
from fractions import Fraction
from flint import arb, arb_mat, fmpq, ctx
ctx.prec=200
from npa_struct import Struct
from exact_psd import is_psd
from families_arb import build, dot, fr
def aq(q): q=Fraction(q); return arb(q.numerator)/q.denominator
def up(x): return x.mid()+x.rad()
def lo(x): return x.mid()-x.rad()
LOG2=arb(2).log()
def guess_funcs(guess):
    out=[]
    if guess[0]=='A':
        x=guess[1]
        for e in (1,-1):
            v=[Fraction(0)]*8; v[x]=Fraction(e,2); out.append((Fraction(1,2),v,(e,)))
    else:
        x,y=guess[1],guess[2]
        for a in (1,-1):
            for b in (1,-1):
                v=[Fraction(0)]*8; v[x]=Fraction(a,4); v[2+y]=Fraction(b,4); v[4+2*x+y]=Fraction(a*b,4)
                out.append((Fraction(1,4),v,(a,b)))
    return out
_S={}
def check_cert(cj,c_expected):
    """exact check of an NPA dual certificate; returns beta (Fraction) with h_Q(c.z) <= beta."""
    L=cj['level']
    if L not in _S: _S[L]=Struct(L)
    s=_S[L]; n=s.n
    c=[Fraction(x) for x in cj['c']]
    assert c==list(c_expected), 'certificate objective mismatch'
    Z=[[Fraction(x) for x in row] for row in cj['Z']]
    assert len(Z)==n and all(len(r)==n for r in Z)
    for i in range(n):
        for j in range(i+1,n): assert Z[i][j]==Z[j][i]
    tgt={k:Fraction(0) for k in range(s.m)}
    for t,k in enumerate(s.kz): tgt[k]=-c[t]
    for k in range(s.m):
        if k==s.k1: continue
        assert sum(Z[i][j] for (i,j) in s.cls[k])==tgt[k]
    ok,_=is_psd([[fmpq(x.numerator,x.denominator) for x in row] for row in Z]); assert ok,'Z not PSD'
    beta=sum(Z[i][i] for i in range(n)); assert beta==Fraction(cj['beta'])
    return beta,L
def zvec_arb(st):
    ang=[aq(x) for x in st['ang_q']]; psi=[aq(x) for x in st['psi_q']]
    nrm=sum((p*p for p in psi),arb(0)).sqrt(); psi=[p/nrm for p in psi]
    def ob(spec,a):
        if spec is not None: return arb_mat([[spec,0],[0,spec]])
        return arb_mat([[a.cos(),a.sin()],[a.sin(),-a.cos()]])
    I=arb_mat([[1,0],[0,1]])
    kron=lambda M,N: arb_mat([[M[i//2,k//2]*N[i%2,k%2] for k in range(4)] for i in range(4)])
    Ax=[kron(ob(st['sa'][x],ang[x]),I) for x in (0,1)]; By=[kron(I,ob(st['sb'][y],ang[2+y])) for y in (0,1)]
    v=arb_mat([[p] for p in psi]); vt=v.transpose(); e=lambda M:(vt*M*v)[0,0]
    return [e(Ax[0]),e(Ax[1]),e(By[0]),e(By[1])]+[e(Ax[x]*By[y]) for x in (0,1) for y in (0,1)]
def in_local(w):
    a=w[:2]; b=w[2:4]; E=w[4:]
    for x,y,s,t in itertools.product((0,1),(0,1),(1,-1),(1,-1)):
        if not (1+s*a[x]+t*b[y]+s*t*E[2*x+y])>=0: return False
    tot=E[0]+E[1]+E[2]+E[3]
    return all((tot-2*E[k])<=2 and (tot-2*E[k])>=-2 for k in range(4))
def verify(d):
    c=d['case']; guess=tuple(c['guess']); G=guess_funcs(guess); V=fr(c['v']); kind=c['kind']; spec=c['spec']
    p,g,order=build(spec)
    U=d['upper']; levels=[]
    # ---------- upper bound on P_guess ----------
    if kind=='full':
        pv=[V*x for x in p]; ft=[Fraction(x) for x in U['f']]; terms=[]
        for (c0,v,lab),pe in zip(G,U['per_e']):
            beta,L=check_cert(pe['cert'],[v[j]+ft[j] for j in range(8)]); levels.append(L); terms.append(aq(c0)+aq(beta))
        fp=dot([aq(x) for x in ft],pv); Pup=max(up(t-fp) for t in terms)
    elif U.get('trivial'):
        Pup=1.0
    else:
        if kind=='cone':
            gs=d['gstar']; lam=[Fraction(x) for x in gs[4:].split(',')]
            assert all(x>=0 for x in lam) and sum(lam)==1, 'g* not a convex combination of the hexagon vertices'
            sp=dict(spec); sp['g']=gs; p,g,order=build(sp)
        om=V*dot(g,p); mu=Fraction(U['mu']); assert mu>=0
        gt=[Fraction(x) for x in U['gtilde']]
        err=sum((abs(x-aq(q)) for x,q in zip(g,gt)),arb(0)); terms=[]
        for (c0,v,lab),pe in zip(G,U['per_e']):
            beta,L=check_cert(pe['cert'],[v[j]+mu*gt[j] for j in range(8)]); levels.append(L)
            terms.append(aq(c0)+aq(beta)+aq(mu)*err-aq(mu)*om)
        Pup=max(up(t) for t in terms)
    Pup=min(float(Pup),1.0) if not isinstance(Pup,float) else min(Pup,1.0)
    # ---------- lower bound on P_guess (explicit attack) ----------
    A=d['attack']; pis=[Fraction(x) for x in A['weights']]
    assert all(x>=0 for x in pis) and sum(pis)==1
    zs=[zvec_arb(st) for st in A['strategies']]
    gi={lab:(c0,v) for c0,v,lab in G}
    Pm=sum((aq(pi)*(aq(gi[tuple(e)][0])+dot([aq(x) for x in gi[tuple(e)][1]],z)) for pi,e,z in zip(pis,A['guesses'],zs)),arb(0))
    zbar=[sum((aq(pi)*z[j] for pi,z in zip(pis,zs)),arb(0)) for j in range(8)]
    ok=True
    if kind=='single':
        p0,g0,_=build(spec); ok=bool((dot(g0,zbar)-V*dot(g0,p0))>0)
    elif kind=='cone':
        p0,_,order0=build({k:x for k,x in spec.items() if k!='g'})
        ok=all(bool((dot(gv,zbar)-V*dot(gv,p0))>0) for gv in order0)
    else:
        T=aq(A['t']); pv=[V*x for x in p]; w=[(pv[j]-(1-T)*zbar[j])/T for j in range(8)]
        ok=in_local(w); Pm=(1-T)*Pm+T
    Plo=float(lo(Pm)) if ok else None
    Hlo=-(arb(Pup).log()/LOG2) if Pup<1 else arb(0)
    Hhi=(-(arb(Plo).log()/LOG2)) if Plo else None
    return {'id':c['id'],'H_lo':float(lo(Hlo)) if Pup<1 else 0.0,'H_hi':float(up(Hhi)) if Hhi is not None else None,
            'P_up':Pup,'P_lo':Plo,'attack_ok':ok,'levels':levels,'npa3_float':c.get('npa3')}
if __name__=='__main__':
    pats=sys.argv[1:]; res=[]
    for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'certs','*.json'))):
        cid=os.path.basename(fn)[:-5]
        if pats and not any(s in cid for s in pats): continue
        r=verify(json.load(open(fn))); res.append(r)
        print('%-34s H in [%.6f, %s]  levels=%s  npa3=%s'%(r['id'],r['H_lo'],'%.6f'%r['H_hi'] if r['H_hi'] is not None else 'FAIL',
              ','.join(sorted(set(r['levels']))) or '-', '%.4f'%r['npa3_float'] if r['npa3_float'] is not None else '-'),flush=True)
    os.makedirs('out',exist_ok=True)
    if not pats: json.dump(res,open('out/verified.json','w'),indent=1)
