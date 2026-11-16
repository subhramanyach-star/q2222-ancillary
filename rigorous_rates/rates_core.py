"""Two-sided rigorous bounds on Eve's guessing probability (hence H_min) in the 2222 scenario.
UPPER bound on P_guess (=> LOWER bound on H_min), valid for all quantum systems:
  Lemma 1 (single score g.z >= omega):  P_g <= max_e h_Q(Pi_e + mu g) - mu*omega   for every mu >= 0,
  Lemma 2 (full statistics z = p):       P_g <= max_e h_Q(Pi_e + f)   - f.p        for every f in R^8,
  with h_Q bounded by exact-rational NPA certificates (cert.py) and irrational data handled in ball arithmetic.
LOWER bound on P_guess (=> UPPER bound on H_min): an explicit attack = finite mixture of qubit strategies with
  Eve knowing the label, reproducing the constraint (verified in ball arithmetic)."""
import numpy as np, json, time
from fractions import Fraction
from scipy.optimize import minimize_scalar, linprog
from flint import arb
import cert, qubit
from families_arb import dot, DET
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
def fl(v): return np.array([float(x) for x in v])
def mid(g): return np.array([float(x.mid()) for x in g])
def dyad(g): return [Fraction(float(x.mid())) for x in g]
def l1err(g,gt): return sum((abs(x-arb(q.numerator)/q.denominator) for x,q in zip(g,gt)),arb(0))
def aq(q): q=Fraction(q); return arb(q.numerator)/q.denominator
def hq_float(c,const=0.0): return qubit.qmax(c,const,starts=6)[0]
# ---------------- single score ----------------
def dual_curve(guess,gm,omega_f,mu_hi=20.0):
    G=guess_funcs(guess)
    def F(mu): return max(hq_float(fl(v)+mu*gm,float(c0)) for c0,v,_ in G)-mu*omega_f
    r=minimize_scalar(F,bounds=(0,mu_hi),method='bounded',options={'xatol':1e-7})
    return r.x,r.fun
def certify_hq(c,levels=('1+AB','2','3'),tol=2e-8):
    """returns (beta Fraction, level, cert json); lowest level whose bound is within tol of the qubit value."""
    qv=hq_float(fl(c)); best=None
    for L in levels:
        cd=cert.certify(c,L)
        if cd is None: continue
        b=cd['beta']
        if best is None or b<best[0]: best=(b,L,cert.to_json(cd))
        if float(b)-qv<tol: break
    return best[0],best[1],best[2],qv
def upper_single(guess,g,omega,mu):
    """rigorous upper bound on P_g for score g.z>=omega (g, omega arb), at rational mu."""
    mu=Fraction(mu); gt=dyad(g); err=l1err(g,gt); G=guess_funcs(guess); terms=[]; certs=[]
    for c0,v,lab in G:
        c=[v[j]+mu*gt[j] for j in range(8)]
        beta,L,cj,qv=certify_hq(c)
        terms.append(aq(c0)+aq(beta)+aq(mu)*err); certs.append({'label':lab,'level':L,'beta_float':float(beta),'qubit':qv,'cert':cj})
    P=max(terms,key=lambda x:x.mid())
    for t in terms: P=P if (P-t)>0 or P is t else (P if P.mid()>=t.mid() else t)
    Pup=max(float((t).mid())+float(t.rad()) for t in terms)
    Pball=[t-aq(mu)*omega for t in terms]
    upper=max(float(b.mid())+float(b.rad()) for b in Pball)
    return upper,{'mu':str(mu),'gtilde':[str(x) for x in gt],'l1err_upper':float(err.mid())+float(err.rad()),'per_e':certs}
# ---------------- attacks (column generation) ----------------
def colgen(guess,rows_ge=(),rows_eq=None,rhs_eq=None,iters=60,tol=1e-10,init=None):
    """maximise sum_k pi_k Pi_{e_k}(z_k) over qubit columns s.t. sum pi_k r.z_k >= w for (r,w) in rows_ge,
    sum pi_k z_k = rhs_eq if rows_eq.  Returns columns, weights, value (float)."""
    G=guess_funcs(guess); cols=[]
    def add(st):
        z=qubit.zvec(st); vals=[float(c0)+fl(v)@z for c0,v,_ in G]; k=int(np.argmax(vals))
        cols.append({'st':st,'z':z,'e':k,'val':vals[k]})
    # seeds: deterministic points (as trivial blocks) and init strategies
    for sa in [(1,1),(1,-1),(-1,1),(-1,-1)]:
        for sb in [(1,1),(1,-1),(-1,1),(-1,-1)]:
            add({'sa':sa,'sb':sb,'ang':np.zeros(4),'psi':np.array([1.,0,0,0])})
    for st in (init or []): add(st)
    for r,w in rows_ge: add(qubit.qmax(np.asarray(r),starts=6)[1])
    hist=[]
    for it in range(iters):
        Z=np.array([c['z'] for c in cols]).T; val=np.array([c['val'] for c in cols]); K=len(cols)
        Aub=[];bub=[]
        for r,w in rows_ge: Aub.append(-(np.asarray(r)@Z)); bub.append(-w)
        Aeq=[np.ones(K)]; beq=[1.0]
        if rows_eq is not None:
            for j in range(8): Aeq.append(Z[j]); beq.append(rhs_eq[j])
        res=linprog(-val,A_ub=np.array(Aub) if Aub else None,b_ub=bub if bub else None,A_eq=np.array(Aeq),b_eq=beq,bounds=(0,None),method='highs',options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
        if res.status!=0:
            if it==0 and rows_eq is None: return None
            return None
        mub=res.ineqlin.marginals if Aub else []; meq=res.eqlin.marginals
        # pricing functional: Pi_e(z) + sum mub_i*(-r_i.z) + meq0 + meq[1:].z
        lin=np.zeros(8)
        for (r,w),m in zip(rows_ge,mub): lin-=m*np.asarray(r)
        if rows_eq is not None: lin+=meq[1:]
        cand=[]
        for c0,v,_ in G:
            cand+=qubit.qmax(fl(v)+lin,float(c0),starts=8,return_all=True,rng=np.random.default_rng(it))
        cand.sort(key=lambda t:-t[0])
        rc=cand[0][0]+meq[0]
        hist.append((-res.fun,rc))
        if rc<tol: break
        if len(hist)>8 and hist[-8][0]>=hist[-1][0]-1e-15: break
        added=[]
        for hv,st in cand:
            if hv+meq[0]<max(tol,0.2*rc): break
            z=qubit.zvec(st)
            if any(np.abs(z-a).max()<1e-7 for a in added): continue
            added.append(z); add(st)
            if len(added)>=6: break
    pi=res.x; keep=[k for k in range(len(pi)) if pi[k]>1e-13]
    return {'cols':[cols[k] for k in keep],'pi':pi[keep],'value':-res.fun,'hist':hist,'dual_lin':lin,'dual_nu':meq[0],'dual_ub':list(mub)}
def verify_attack_single(guess,att,g,omega,margin=0):
    """ball-arithmetic check that the mixture satisfies g.z>=omega; returns rigorous lower bound on P_g or None."""
    G=guess_funcs(guess)
    pis=[Fraction(float(p)).limit_denominator(10**12) for p in att['pi']]; s=sum(pis); pis=[p/s for p in pis]
    zs=[qubit.zvec_arb(qubit.freeze(c['st'])) for c in att['cols']]
    score=sum((aq(p)*dot(g,z) for p,z in zip(pis,zs)),arb(0))
    ok=(score-omega)>0
    P=sum((aq(p)*(aq(G[c['e']][0])+dot([aq(x) for x in G[c['e']][1]],z)) for p,c,z in zip(pis,att['cols'],zs)),arb(0))
    return (float(P.mid())-float(P.rad())) if ok else None, {'score_minus_omega':float((score-omega).mid()),'weights':[str(p) for p in pis],'strategies':[qubit.freeze(c['st']) for c in att['cols']],'guesses':[list(G[c['e']][2]) for c in att['cols']]}
# ---------------- full statistics and multi-row attacks ----------------
def upper_full(guess,f,p):
    """Lemma 2 with dyadic f~ (exact): P_g <= max_e h_Q(Pi_e+f~) - f~.p   (p arb)."""
    ft=[Fraction(float(x)) for x in f]; G=guess_funcs(guess); terms=[]; certs=[]
    for c0,v,lab in G:
        c=[v[j]+ft[j] for j in range(8)]
        beta,L,cj,qv=certify_hq(c)
        terms.append(aq(c0)+aq(beta)); certs.append({'label':lab,'level':L,'beta_float':float(beta),'qubit':qv,'cert':cj})
    fp=dot([aq(x) for x in ft],p)
    upper=max(float((t-fp).mid())+float((t-fp).rad()) for t in terms)
    return upper,{'f':[str(x) for x in ft],'per_e':certs}
def in_local(w):
    """Fine's theorem: z is local iff all 16 probabilities >=0 and all 8 CHSH inequalities hold (ball check)."""
    a=w[:2]; b=w[2:4]; E=w[4:]
    for x in (0,1):
        for y in (0,1):
            for s in (1,-1):
                for t in (1,-1):
                    if not (1+s*a[x]+t*b[y]+s*t*E[2*x+y])>=0: return False
    tot=E[0]+E[1]+E[2]+E[3]
    for k in range(4):
        q=tot-2*E[k]
        if not (q<=2 and q>=-2): return False
    return True
def verify_attack(guess,att,rows_ge=(),p_eq=None,t=None):
    """rows_ge: list of (g arb, omega arb); p_eq: arb behaviour, mixture approximates p_eq/(1-t) and the
    remainder is absorbed by a local (hence perfectly guessable) component of weight t."""
    G=guess_funcs(guess)
    pis=[Fraction(float(p)).limit_denominator(10**14) for p in att['pi']]; s=sum(pis); pis=[p/s for p in pis]
    zs=[qubit.zvec_arb(qubit.freeze(c['st'])) for c in att['cols']]
    Pm=sum((aq(p)*(aq(G[c['e']][0])+dot([aq(x) for x in G[c['e']][1]],z)) for p,c,z in zip(pis,att['cols'],zs)),arb(0))
    info={'weights':[str(p) for p in pis],'strategies':[qubit.freeze(c['st']) for c in att['cols']],'guesses':[list(G[c['e']][2]) for c in att['cols']]}
    ok=True; slack=[]
    for g,om in rows_ge:
        sc=sum((aq(p)*dot(g,z) for p,z in zip(pis,zs)),arb(0)); ok=ok and bool((sc-om)>0); slack.append(float((sc-om).mid()))
    info['score_slacks']=slack
    if p_eq is not None:
        T=aq(t); y=[sum((aq(p)*z[j] for p,z in zip(pis,zs)),arb(0)) for j in range(8)]
        w=[(p_eq[j]-(1-T)*y[j])/T for j in range(8)]
        loc=in_local(w); ok=ok and loc; info['t']=str(t); info['absorbed_w_maxabs']=max(abs(float(x.mid())) for x in w); info['w_local']=loc
        Pm=(1-T)*Pm+T
    info['ok']=ok
    return (float(Pm.mid())-float(Pm.rad())) if ok else None, info
