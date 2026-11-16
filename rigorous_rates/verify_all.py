#!/usr/bin/env python3
"""STANDALONE VERIFIER (no SDP/LP solver, no floating-point optimisation).
Needs only python3 + python-flint (pip install python-flint).  For every certs/<id>.json it re-derives, in exact
rational and ball (arb, 200-bit) arithmetic:
  lower bound on H_min : Lemma 1 / Lemma 2 of PROOFS.md with exact-rational NPA certificates (PSD by exact LDL)
  upper bound on H_min : the stored explicit qubit attack (mixture weights, dyadic angles and state) reproduces
                         the constraint (score >= omega, all six cone vertices, or p_v exactly via a local remainder)
Prints the verified interval [H_lo, H_hi] for each case and writes out/verified.json.
Usage: python3 verify_all.py [--lower-bounds-only] [id-substring ...]
Default (two-sided): exit 0 iff a nonempty collection is verified and every case has
a present, valid, successful attack (attack_ok and a finite H_hi).  The explicitly named
--lower-bounds-only mode skips attack evaluation entirely (it does not read d['attack']),
may report success without an attack, and never claims two-sided verification."""
import sys, os, json, glob, itertools
from fractions import Fraction
from flint import arb, arb_mat, fmpq, ctx
ctx.prec=200

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

from npa_struct import Struct
from exact_psd import is_psd
from families_arb import build, dot, fr
from arb_bounds import (
    aq, LOG2, certified_max_upper, strictly_below, strictly_above,
    hmin_lower_from_P_upper, hmin_upper_from_P_lower,
    format_lower, format_upper, certified_float_lower, certified_float_upper,
    definitely_gt,
)
# Certified replacements for the former mid/rad float conversions.
# Kept under the old names so that importers can be updated gradually.
def up(x):
    return strictly_above(x)
def lo(x):
    return strictly_below(x)
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
    if c!=list(c_expected):
        raise ValueError('certificate objective mismatch')
    Z=[[Fraction(x) for x in row] for row in cj['Z']]
    if not (len(Z)==n and all(len(r)==n for r in Z)):
        raise ValueError('certificate matrix has the wrong dimension')
    for i in range(n):
        for j in range(i+1,n):
            if Z[i][j]!=Z[j][i]:
                raise ValueError('certificate matrix is not symmetric')
    tgt={k:Fraction(0) for k in range(s.m)}
    for t,k in enumerate(s.kz): tgt[k]=-c[t]
    for k in range(s.m):
        if k==s.k1: continue
        if sum(Z[i][j] for (i,j) in s.cls[k])!=tgt[k]:
            raise ValueError('certificate linear constraint %d failed' % k)
    ok,_=is_psd([[fmpq(x.numerator,x.denominator) for x in row] for row in Z])
    if not ok:
        raise ValueError('Z not PSD')
    beta=sum(Z[i][i] for i in range(n))
    if beta!=Fraction(cj['beta']):
        raise ValueError('certificate beta does not match the trace of Z')
    return beta,L

def _require_obs_specs(specs):
    if not isinstance(specs, (list, tuple)) or len(specs)!=2:
        raise ValueError('need two observable specifications per party')
    for s in specs:
        if s is not None and s not in (-1, +1):
            raise ValueError('observable spec must be None, -1, or +1')

def check_attack_input(A, G, kind):
    """Physical-data checks on an explicit attack, before it is evaluated."""
    if not all(k in A for k in ('weights', 'strategies', 'guesses')):
        raise ValueError('attack is missing weights, strategies, or guesses')
    weights, strategies, guesses = A['weights'], A['strategies'], A['guesses']
    if not (len(weights)==len(strategies)==len(guesses) and len(weights)>0):
        raise ValueError('weights, strategies, and guesses must have equal nonzero length')
    pis=[Fraction(x) for x in weights]
    if not all(x>=0 for x in pis) or sum(pis)!=1:
        raise ValueError('attack weights must be a nonnegative partition of unity')
    valid={lab for _,_,lab in G}
    for st,e in zip(strategies, guesses):
        if st.get('ang_q') is None or len(st['ang_q'])!=4:
            raise ValueError('strategy must have four angles')
        if st.get('psi_q') is None or len(st['psi_q'])!=4:
            raise ValueError('strategy must have four state components')
        if all(Fraction(x)==0 for x in st['psi_q']):
            raise ValueError('state vector is zero')
        _require_obs_specs(st.get('sa'))
        _require_obs_specs(st.get('sb'))
        lab=tuple(e)
        if lab not in valid:
            raise ValueError('invalid guessing label')
    if kind=='full' or 't' in A:
        if 't' not in A:
            raise ValueError('local remainder is missing t')
        t=Fraction(A['t'])
        if not (t>0 and t<=1):
            raise ValueError('local remainder weight t must satisfy 0 < t <= 1')
    return pis

def zvec_arb(st):
    ang=[aq(x) for x in st['ang_q']]; psi=[aq(x) for x in st['psi_q']]
    if len(st['ang_q'])!=4 or len(st['psi_q'])!=4:
        raise ValueError('strategy must have four angles and four state components')
    if all(Fraction(x)==0 for x in st['psi_q']):
        raise ValueError('state vector is zero')
    nrm=sum((p*p for p in psi),arb(0)).sqrt(); psi=[p/nrm for p in psi]
    def ob(spec,a):
        if spec is not None:
            if spec not in (-1, +1):
                raise ValueError('observable spec must be None, -1, or +1')
            return arb_mat([[spec,0],[0,spec]])
        return arb_mat([[a.cos(),a.sin()],[a.sin(),-a.cos()]])
    _require_obs_specs(st.get('sa')); _require_obs_specs(st.get('sb'))
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
def verify(d, lower_bounds_only=False):
    c=d['case']; guess=tuple(c['guess']); G=guess_funcs(guess); V=fr(c['v']); kind=c['kind']; spec=c['spec']
    p,g,order=build(spec)
    U=d['upper']; levels=[]
    # ---------- upper bound on P_guess ----------
    if kind=='full':
        if len(U.get('per_e', [])) != len(G):
            raise ValueError('Missing or extra guessing certificates')
        pv=[V*x for x in p]; ft=[Fraction(x) for x in U['f']]
        if len(ft) != 8 or len(pv) != 8:
            raise ValueError('Invalid Bell-vector dimension')
        terms=[]
        for (c0,v,lab),pe in zip(G,U['per_e']):
            beta,L=check_cert(pe['cert'],[v[j]+ft[j] for j in range(8)]); levels.append(L); terms.append(aq(c0)+aq(beta))
        fp=dot([aq(x) for x in ft],pv); Pup=certified_max_upper([t-fp for t in terms])
    elif U.get('trivial'):
        Pup=arb(1)
    else:
        if kind=='cone':
            gs=d['gstar']; lam=[Fraction(x) for x in gs[4:].split(',')]
            if not (all(x>=0 for x in lam) and sum(lam)==1):
                raise ValueError('g* not a convex combination of the hexagon vertices')
            sp=dict(spec); sp['g']=gs; p,g,order=build(sp)
        om=V*dot(g,p); mu=Fraction(U['mu'])
        if mu<0:
            raise ValueError('Negative threshold multiplier')
        gt=[Fraction(x) for x in U['gtilde']]
        if len(gt) != 8 or len(g) != 8:
            raise ValueError('Invalid Bell-vector dimension')
        if len(U.get('per_e', [])) != len(G):
            raise ValueError('Missing or extra guessing certificates')
        err=sum((abs(x-aq(q)) for x,q in zip(g,gt)),arb(0)); terms=[]
        for (c0,v,lab),pe in zip(G,U['per_e']):
            beta,L=check_cert(pe['cert'],[v[j]+mu*gt[j] for j in range(8)]); levels.append(L)
            terms.append(aq(c0)+aq(beta)+aq(mu)*err-aq(mu)*om)
        Pup=certified_max_upper(terms)
    if Pup >= 1 or not (Pup < 1):
        # If the enclosure is not certified < 1, cap the guessing probability at 1.
        if not (Pup < 1):
            Pup=arb(1)
    # ---------- lower bound on P_guess (explicit attack) ----------
    Hlo_ball=hmin_lower_from_P_upper(Pup)
    if lower_bounds_only:
        # Do not read or validate d['attack'].  This mode never reports a
        # two-sided interval.
        return {'id':c['id'],
                'H_lo':certified_float_lower(Hlo_ball),
                'H_hi':None,
                'H_lo_dec':format_lower(Hlo_ball, 6),
                'H_hi_dec':None,
                'P_up':certified_float_upper(Pup),'P_lo':None,'attack_ok':False,
                'levels':levels,'npa3_float':c.get('npa3'),
                'lower_bounds_only':True}
    if 'attack' not in d:
        raise ValueError('two-sided verification requires a present attack')
    A=d['attack']
    pis=check_attack_input(A, G, kind)
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
        t=Fraction(A['t']); T=aq(t); pv=[V*x for x in p]; w=[(pv[j]-(1-T)*zbar[j])/T for j in range(8)]
        ok=in_local(w); Pm=(1-T)*Pm+T
    if not ok:
        Plo=None
        Hhi_ball=None
    else:
        Plo_ball=strictly_below(Pm)
        Plo=certified_float_lower(Plo_ball)
        Hhi_ball=hmin_upper_from_P_lower(Plo_ball)
        if Hhi_ball is not None and not Hhi_ball.is_finite():
            Hhi_ball=None
            ok=False
    return {'id':c['id'],
            'H_lo':certified_float_lower(Hlo_ball),
            'H_hi':certified_float_upper(Hhi_ball) if Hhi_ball is not None else None,
            'H_lo_dec':format_lower(Hlo_ball, 6),
            'H_hi_dec':format_upper(Hhi_ball, 6) if Hhi_ball is not None else None,
            'P_up':certified_float_upper(Pup),'P_lo':Plo,'attack_ok':ok,'levels':levels,'npa3_float':c.get('npa3'),
            'lower_bounds_only':False}
if __name__=='__main__':
    lower_only=False; pats=[]
    for a in sys.argv[1:]:
        if a in ('--lower-bounds-only','--lower_bounds_only'):
            lower_only=True
        else:
            pats.append(a)
    res=[]
    for fn in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'certs','*.json'))):
        cid=os.path.basename(fn)[:-5]
        if pats and not any(s in cid for s in pats): continue
        r=verify(json.load(open(fn)), lower_bounds_only=lower_only); res.append(r)
        if lower_only:
            print('%-34s H_lo >= %s  (lower bounds only)  levels=%s  npa3=%s'%(
                r['id'], r['H_lo_dec'],
                ','.join(sorted(set(r['levels']))) or '-',
                '%.4f'%r['npa3_float'] if r['npa3_float'] is not None else '-'),flush=True)
        else:
            print('%-34s H in [%s, %s]  levels=%s  npa3=%s'%(r['id'],r['H_lo_dec'], r['H_hi_dec'] if r['H_hi_dec'] is not None else 'FAIL',
                  ','.join(sorted(set(r['levels']))) or '-', '%.4f'%r['npa3_float'] if r['npa3_float'] is not None else '-'),flush=True)
    os.makedirs('out',exist_ok=True)
    if not pats and not lower_only: json.dump(res,open('out/verified.json','w'),indent=1)
    def _finite_hi(r):
        h=r['H_hi']
        if h is None:
            return False
        try:
            import math
            return math.isfinite(h)
        except (TypeError, ValueError):
            return False
    if lower_only:
        allok=bool(res) and all(r['H_lo'] is not None for r in res)
        print('RESULT LOWER BOUNDS ONLY =', allok)
    else:
        allok=bool(res) and all(r['attack_ok'] and _finite_hi(r) for r in res)
        print('RESULT TWO-SIDED ALL OK =', allok)
    sys.exit(0 if allok else 1)
