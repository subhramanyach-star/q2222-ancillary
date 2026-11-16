"""Rigorous two-sided H_min(A1|X=1,E) statements under detector efficiency eta (exact rational eta, v).
LOSER (every CHSH-based protocol):  S_obs <= S_max(eta) <= S_up := c0 + beta, beta = exact-rational NPA certificate
  for the loss-transformed CHSH functional (Lemma 3 of ../rigorous_rates/PROOFS.md; valid in every dimension), and an
  explicit qubit attack with CHSH >= S_up has guessing probability >= P_lo (ball arithmetic)  =>  H_CHSH <= -log2 P_lo.
WINNER (a stated protocol): honest device (dyadic theta, angles) at efficiency eta and visibility v reaches score
  g.z >= omega (checked in arb); Lemma 1 with exact-rational NPA certificates gives P_guess <= P_up for EVERY quantum
  device reaching that score  =>  H_min >= -log2 P_up.  Witnesses: tilted CHSH alpha<A0>+CHSH (dyadic alpha), and a
  witness in the Thm-2 cone C_theta' (exact convex combination of the exact polygon vertices, arb).
Also an explicit attack for the winner (two-sided interval).  Writes certs/<id>.json + out/cert_results.jsonl."""
import sys, json, os, time, warnings
warnings.simplefilter('ignore')
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rigorous_rates'))
import numpy as np
from fractions import Fraction as F
from flint import arb, arb_mat, ctx
ctx.prec=200
import rates_core as R, qubit
from families_arb import dot, DET as DETA
import families_arb as FA
from scipy.optimize import minimize_scalar
from core import PGm, z_honest, G_CHSH
import cone as CN
LOG2=arb(2).log()
def aq(q): q=F(q); return arb(q.numerator)/q.denominator
def dy(x,bits=40): return F(int(round(float(x)*2**bits)),2**bits)
def ball_up(x): return float(x.mid())+float(x.rad())
def ball_lo(x): return float(x.mid())-float(x.rad())
GUESS=('A',1)
G_CHSH_F=[F(x) for x in (0,0,0,0,1,1,1,-1)]
def loss_functional_F(g,eta):
    """exact: g.z_eta = c0 + g_eta.z_ideal for the binned behaviour (eta Fraction)."""
    E=g[4:]; c0=(1-eta)*(g[0]+g[1]+g[2]+g[3])+(1-eta)**2*sum(E)
    ga=[eta*g[x]+eta*(1-eta)*(E[2*x]+E[2*x+1]) for x in (0,1)]
    gb=[eta*g[2+y]+eta*(1-eta)*(E[y]+E[2+y]) for y in (0,1)]
    return c0,ga+gb+[eta**2*e for e in E]
def honest_arb(par_q,eta,v):
    th=aq(par_q[0]); a=[aq(x) for x in par_q[1:5]]; C=(2*th).cos(); S=(2*th).sin(); eta=aq(eta); v=aq(v)
    A=[v*a[0].cos()*C,v*a[1].cos()*C]; B=[v*a[2].cos()*C,v*a[3].cos()*C]
    E=[v*(a[x].cos()*a[2+y].cos()+a[x].sin()*a[2+y].sin()*S) for x in (0,1) for y in (0,1)]
    A2=[eta*t+(1-eta) for t in A]; B2=[eta*t+(1-eta) for t in B]
    E2=[eta*eta*E[2*x+y]+eta*(1-eta)*(A[x]+B[y])+(1-eta)**2 for x in (0,1) for y in (0,1)]
    return A2+B2+E2
def H_of(P): return -(P.log()/LOG2)
# ------------------------------------------------------------------ loser
def chsh_side(eta,v=F(1)):
    """visibility v: z_eta = loss(v z_ideal) => S = c0 + v g_eta.z_ideal <= c0 + v beta (beta >= 0)."""
    c0,ge=loss_functional_F(G_CHSH_F,eta)
    beta,L,cj,qv=R.certify_hq(ge)
    assert beta>=0
    S_up=c0+v*beta
    # explicit attack at CHSH >= S_up (S_up is rational)
    att=R.colgen(GUESS,rows_ge=[(G_CHSH.copy(),float(S_up)+1e-11)],iters=120)
    lo,ainfo=R.verify_attack(GUESS,att,rows_ge=[([aq(x) for x in G_CHSH_F],aq(S_up))])
    Hup=H_of(aq(F(lo).limit_denominator(10**15)) - arb(10)**-15) if lo else None
    return dict(S_up=str(S_up),S_up_f=float(S_up),S_max_qubit=float(c0)+float(v)*qv,level=L,cert=cj,c0=str(c0),
                H_up=(ball_up(Hup) if lo else None),P_lo=lo,attack=ainfo,
                H_closed=float(-np.log2(0.5+0.5*np.sqrt(max(0,2-float(S_up)**2/4)))))
# ------------------------------------------------------------------ winner
def best_mu(g,om):
    gm=R.mid(g); r=R.dual_curve(GUESS,gm,float(om.mid()) if isinstance(om,arb) else float(om)); return r[0]
def winner(par,eta,v,g_arb,label,mu_hint=None):
    par_q=[dy(x) for x in par[:5]]
    z=honest_arb(par_q,eta,v); sc=dot(g_arb,z)
    omega=F(int(np.floor((float(sc.mid())-1e-12)*2**44)),2**44)
    assert bool(sc-aq(omega)>0), 'honest score below threshold'
    mu=best_mu(g_arb,aq(omega))
    up,uinfo=R.upper_single(GUESS,g_arb,aq(omega),F(mu).limit_denominator(10**8))
    att=R.colgen(GUESS,rows_ge=[(R.mid(g_arb),float(omega)+1e-11)],init=[dict(sa=(None,None),sb=(None,None),ang=np.array(par[1:5]),psi=np.array([np.cos(par[0]),0,0,np.sin(par[0])]))],iters=120)
    lo,ainfo=R.verify_attack(GUESS,att,rows_ge=[(g_arb,aq(omega))])
    return dict(label=label,par_q=[str(x) for x in par_q],omega=str(omega),honest_score_minus_omega=float((sc-aq(omega)).mid()),
                P_up=up,H_lo=float(-np.log2(up)),P_lo=lo,H_up=(float(-np.log2(lo)) if lo else None),upper=uinfo,attack=ainfo)
def tilted_g(alpha): return [aq(alpha)]+[arb(0)]*3+[arb(1),arb(1),arb(1),arb(-1)]
def best_alpha(z):
    pg=PGm(1,0,guess=GUESS)
    f=lambda a: pg(z,Rg=np.array([[a,0,0,0,1,1,1,-1.]]))
    grid=np.linspace(-0.5,2.5,121); vals=[f(a) for a in grid]; k=int(np.argmin(vals))
    r=minimize_scalar(f,bounds=(grid[max(k-1,0)],grid[min(k+1,120)]),method='bounded',options={'xatol':1e-8})
    return (r.x,r.fun) if r.fun<=vals[k] else (grid[k],vals[k])
# ---- exact polygon vertices of C_theta' in arb (active deterministic pairs from the float computation)
def polygon_arb(thq):
    th=aq(thq); p,c0,(s,C,r)=FA.ptheta(th); Gm=FA.gammas(p,c0); m=[1+C,-s*s,1-C,-s*s]
    pf,Gf=CN.polygon(float(thq)); DETf=CN.DET
    verts=[]
    for gf in Gf:
        act=[i for i,L in enumerate(DETf) if abs(gf@L-1)<1e-7]
        # g = sum mu_i Gm_i with mu.m=0, g.p=1, g.L=1 for two active L  (4 linear equations)
        best=None
        for i in range(len(act)):
            for j in range(i+1,len(act)):
                L1,L2=FA.DET[act[i]],FA.DET[act[j]]
                M=arb_mat([m,[dot(g,p) for g in Gm],[dot(g,L1) for g in Gm],[dot(g,L2) for g in Gm]])
                try: mu=M.solve(arb_mat([[0],[1],[1],[1]]))
                except Exception: continue
                g=[sum((mu[k,0]*Gm[k][jj] for k in range(4)),arb(0)) for jj in range(8)]
                if max(abs(float(g[jj].mid())-gf[jj]) for jj in range(8))<1e-6: best=g; break
            if best: break
        assert best is not None, 'vertex not reconstructed'
        verts.append(best)
    return p,verts
def cone_witness(z,thq):
    """best witness in C_theta' for data z: dual weights of the multi-row SDP -> rational convex weights."""
    p,Gf=CN.polygon(float(thq)); pc=PGm(len(Gf),0,guess=GUESS)
    P=pc(z,Rg=Gf); lam=np.maximum(pc.duals()[0],0); lam=lam/lam.sum()
    lamq=[F(x).limit_denominator(10**9) for x in lam]; s=sum(lamq); lamq=[x/s for x in lamq]
    return lamq,P
def run_point(eta,v,recs,tag):
    eta=F(eta).limit_denominator(1000); v=F(v).limit_denominator(1000); out={}
    out['eta']=str(eta); out['v']=str(v); out['guess']=list(GUESS); t0=time.time()
    out['chsh']=chsh_side(eta,v)
    for fam in ('tilted','cone'):
        if fam not in recs: continue
        par=recs[fam]['par']; z=z_honest(par,float(eta),float(v))
        if fam=='tilted':
            a,_=best_alpha(z); aq_=dy(a,20); g=tilted_g(aq_)
            out['tilted']=winner(par,eta,v,g,'tilted'); out['tilted']['alpha']=str(aq_)
        else:
            thq=dy(min(recs['cone']['thetap'],np.pi/4-1e-9),30); lamq,_=cone_witness(z,thq)
            p,verts=polygon_arb(thq)
            g=[sum((aq(l)*V[j] for l,V in zip(lamq,verts)),arb(0)) for j in range(8)]
            out['cone']=winner(par,eta,v,g,'cone'); out['cone']['thetap']=str(thq); out['cone']['lam']=[str(x) for x in lamq]
            out['cone']['g_dot_p_minus_1']=float((dot(g,p)-1).mid())
    out['secs']=round(time.time()-t0,1)
    json.dump(out,open('certs/%s.json'%tag,'w'))
    summ=dict(id=tag,guess=list(GUESS),eta=float(eta),v=float(v),S_up=out['chsh']['S_up_f'],H_chsh_up=out['chsh']['H_up'],H_chsh_closed=out['chsh']['H_closed'])
    for fam in ('tilted','cone'):
        if fam in out: summ[fam]=[out[fam]['H_lo'],out[fam]['H_up']]
    with open('out/cert_results.jsonl','a') as f: f.write(json.dumps(summ)+'\n')
    print(json.dumps(summ),flush=True)
    return out
if __name__=='__main__':
    import glob
    V=float(sys.argv[1]); etas=[float(x) for x in sys.argv[2].split(',')]
    if len(sys.argv)>3:
        gs=sys.argv[3]; GUESS=('AB',int(gs[2]),int(gs[3])) if gs.startswith('AB') else ('A',int(gs[1]))
    gtag='' if GUESS==('A',1) else '_'+''.join(map(str,GUESS))
    done=set()
    if os.path.exists('out/cert_results.jsonl'):
        for l in open('out/cert_results.jsonl'): d=json.loads(l); done.add(d['id'])
    from core import load_opt
    allr=load_opt(V,'' if GUESS==('A',1) else '_'+''.join(map(str,GUESS)))
    for eta in etas:
        tag='eta%s_v%s%s'%(eta,V,gtag)
        if tag in done: continue
        recs={f:allr[(eta,f)] for f in ('tilted','cone') if (eta,f) in allr}
        run_point(eta,V,recs,tag)
