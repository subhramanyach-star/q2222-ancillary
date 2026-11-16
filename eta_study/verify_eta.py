#!/usr/bin/env python3
"""STANDALONE VERIFIER for the detector-efficiency study (no SDP/LP solver; python3 + python-flint + numpy).
For every certs/eta*.json it re-derives in exact rational / 200-bit ball arithmetic:
 (1) CHSH side: the loss-transformed CHSH functional c0 + g_eta.z (exact in eta), the exact-rational NPA certificate
     h_Q(g_eta) <= beta (exact LDL^T), S_up = c0 + beta; the stored explicit qubit attack has CHSH >= S_up; its
     guessing probability P_lo  =>  every CHSH-based protocol at this eta has H_min(A1|X=1,E) <= -log2 P_lo.
 (2) winners (tilted alpha<A0>+CHSH; witness in the Thm-2 cone C_theta'): the honest qubit device (dyadic theta and
     angles, exact eta, v) reaches g.z >= omega; Lemma 1 with exact NPA certificates gives P_guess <= P_up for every
     quantum device reaching omega  =>  H_min >= -log2 P_up;  plus the stored attack (H_min <= -log2 P_lo).
     For the cone witness: g = sum lam_i V_i with rational convex lam, V_i the exact vertices of C_theta' (each checked:
     V_i.p_theta' = 1, V_i.L <= 1 for all 16 deterministic L, V_i in the slice plane).
Prints a table and writes out/verified_eta.json.  A separation at eta is PROVED when winner H_lo > CHSH H_up."""
import sys, os, glob, json
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rigorous_rates')); sys.path.insert(0,HERE)
from fractions import Fraction as F
from flint import arb, ctx
ctx.prec=200
from verify_all import check_cert, zvec_arb, guess_funcs, aq, up, lo, LOG2
from families_arb import dot
import families_arb as FA
GUESS=('A',1); G=guess_funcs(GUESS); GI={lab:(c0,v) for c0,v,lab in G}
def set_guess(gu):
    global GUESS,G,GI
    GUESS=tuple(gu); G=guess_funcs(GUESS); GI={lab:(c0,v) for c0,v,lab in G}
CHSH=[F(x) for x in (0,0,0,0,1,1,1,-1)]
def lossF(g,eta):
    E=g[4:]; c0=(1-eta)*(g[0]+g[1]+g[2]+g[3])+(1-eta)**2*sum(E)
    ga=[eta*g[x]+eta*(1-eta)*(E[2*x]+E[2*x+1]) for x in (0,1)]; gb=[eta*g[2+y]+eta*(1-eta)*(E[y]+E[2+y]) for y in (0,1)]
    return c0,ga+gb+[eta**2*e for e in E]
def honest(par_q,eta,v):
    th=aq(par_q[0]); a=[aq(x) for x in par_q[1:5]]; C=(2*th).cos(); S=(2*th).sin(); eta=aq(eta); v=aq(v)
    A=[v*a[0].cos()*C,v*a[1].cos()*C]; B=[v*a[2].cos()*C,v*a[3].cos()*C]
    E=[v*(a[x].cos()*a[2+y].cos()+a[x].sin()*a[2+y].sin()*S) for x in (0,1) for y in (0,1)]
    return [eta*t+(1-eta) for t in A]+[eta*t+(1-eta) for t in B]+[eta*eta*E[2*x+y]+eta*(1-eta)*(A[x]+B[y])+(1-eta)**2 for x in (0,1) for y in (0,1)]
def attack(A,rows):
    pis=[F(x) for x in A['weights']]; assert all(x>=0 for x in pis) and sum(pis)==1
    zs=[zvec_arb(st) for st in A['strategies']]
    Pm=sum((aq(pi)*(aq(GI[tuple(e)][0])+dot([aq(x) for x in GI[tuple(e)][1]],z)) for pi,e,z in zip(pis,A['guesses'],zs)),arb(0))
    zbar=[sum((aq(pi)*z[j] for pi,z in zip(pis,zs)),arb(0)) for j in range(8)]
    ok=all(bool(dot(g,zbar)-om>0) for g,om in rows)
    return (float(lo(Pm)) if ok else None)
def H(P): return -(arb(P).log()/LOG2)
def upper_P(U,g,omega):
    mu=F(U['mu']); assert mu>=0; gt=[F(x) for x in U['gtilde']]
    err=sum((abs(x-aq(q)) for x,q in zip(g,gt)),arb(0)); terms=[]; lev=[]
    for (c0,v,lab),pe in zip(G,U['per_e']):
        beta,L=check_cert(pe['cert'],[v[j]+mu*gt[j] for j in range(8)]); lev.append(L)
        terms.append(aq(c0)+aq(beta)+aq(mu)*err-aq(mu)*omega)
    return min(max(float(up(t)) for t in terms),1.0),lev
def cone_vertices(thq):
    import cone as CN, numpy as np
    from flint import arb_mat
    th=aq(thq); p,c0,(s,C,r)=FA.ptheta(th); Gm=FA.gammas(p,c0); m=[1+C,-s*s,1-C,-s*s]
    pf,Gf=CN.polygon(float(F(thq))); out=[]
    for gf in Gf:
        act=[i for i,L in enumerate(CN.DET) if abs(gf@L-1)<1e-7]; got=None
        for i in range(len(act)):
            for j in range(i+1,len(act)):
                M=arb_mat([m,[dot(g,p) for g in Gm],[dot(g,FA.DET[act[i]]) for g in Gm],[dot(g,FA.DET[act[j]]) for g in Gm]])
                try: mu=M.solve(arb_mat([[0],[1],[1],[1]]))
                except Exception: continue
                g=[sum((mu[k,0]*Gm[k][jj] for k in range(4)),arb(0)) for jj in range(8)]
                if max(abs(float(g[jj].mid())-gf[jj]) for jj in range(8))<1e-6: got=g; break
            if got: break
        assert got is not None
        assert (dot(got,p)-1).contains(0), 'V.p != 1'
        assert all(bool(dot(got,L)-1<=1e-50) or (dot(got,L)-1).contains(0) for L in FA.DET), 'deterministic cut violated'
        out.append(got)
    return p,out
def verify(d):
    set_guess(d.get('guess',['A',1]))
    eta=F(d['eta']); v=F(d['v']); res={'eta':float(eta),'v':float(v),'guess':list(GUESS)}
    c=d['chsh']; c0,ge=lossF(CHSH,eta)
    beta,L=check_cert(c['cert'],ge); assert beta>=0; Sup=c0+v*beta; assert Sup==F(c['S_up'])
    Plo=attack(c['attack'],[([aq(x) for x in CHSH],aq(Sup))]); assert Plo is not None, 'CHSH attack fails'
    res['S_up']=float(Sup); res['chsh_H_up']=float(up(H(Plo))); res['chsh_level']=L
    for fam in ('tilted','cone'):
        if fam not in d: continue
        w=d[fam]
        if fam=='tilted': g=[aq(F(w['alpha']))]+[arb(0)]*3+[arb(1),arb(1),arb(1),arb(-1)]
        else:
            lam=[F(x) for x in w['lam']]; assert all(x>=0 for x in lam) and sum(lam)==1
            p,verts=cone_vertices(F(w['thetap'])); assert len(verts)==len(lam)
            g=[sum((aq(l)*V[j] for l,V in zip(lam,verts)),arb(0)) for j in range(8)]
        om=F(w['omega']); z=honest(w['par_q'],eta,v); assert bool(dot(g,z)-aq(om)>0), 'honest device misses threshold'
        Pup,lev=upper_P(w['upper'],g,aq(om)); Pl=attack(w['attack'],[(g,aq(om))])
        res[fam]=[float(lo(H(Pup))) if Pup<1 else 0.0, float(up(H(Pl))) if Pl else None]; res[fam+'_levels']=sorted(set(lev))
        res[fam+'_proved_gap']=res[fam][0]-res['chsh_H_up']
    return res
if __name__=='__main__':
    pats=sys.argv[1:]; out=[]
    for fn in sorted(glob.glob(os.path.join(HERE,'certs','eta*.json'))):
        if pats and not any(s in fn for s in pats): continue
        r=verify(json.load(open(fn))); out.append(r)
        s='%-10s eta=%.4f v=%.3f  S_up=%.6f  CHSH H<=%.6f'%(''.join(map(str,r['guess'])),r['eta'],r['v'],r['S_up'],r['chsh_H_up'])
        for fam in ('tilted','cone'):
            if fam in r: s+='  %s H in [%.6f,%.6f] gap %+.6f %s'%(fam,r[fam][0],r[fam][1] or -1,r[fam+'_proved_gap'],'PROVED' if r[fam+'_proved_gap']>0 else '')
        print(s,flush=True)
    if not pats: json.dump(out,open(os.path.join(HERE,'out','verified_eta.json'),'w'),indent=1)
