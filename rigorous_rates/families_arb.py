"""Ball-arithmetic (python-flint arb, 200 bits) versions of the targets and Bell functionals of
anc_submission_v7/rates_numerics/families.py.  Every quantity is an arb ball that provably contains the exact value.
Coordinates z=(a0,a1,b0,b1,E00,E01,E10,E11)."""
from flint import arb, arb_mat, ctx
from fractions import Fraction
import itertools
ctx.prec=200
PI=arb.pi(); S2=arb(2).sqrt()
def A(x): return x if isinstance(x,arb) else arb(Fraction(x).numerator)/Fraction(x).denominator
def vec(lst): return [A(x) for x in lst]
def dot(a,b): return sum((x*y for x,y in zip(a,b)),arb(0))
def detpt(A0,A1,B0,B1): return [arb(v) for v in (A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1)]
DET=[detpt(*L) for L in itertools.product([1,-1],repeat=4)]
def theta_of(name):
    return {'pi/8':PI/8,'3pi/16':3*PI/16,'27.76deg':arb(2776)/18000*PI,'pi/4':PI/4}[name]
def ptheta(th):
    s=(2*th).sin(); C=(2*th).cos(); r=(1+s*s).sqrt()
    p=[C,arb(0),C/r,C/r,1/r,1/r,s*s/r,-s*s/r]; c0=(1-s*s)/(1+s*s)
    return p,c0,(s,C,r)
def gammas(p,c0):
    det=1-c0*c0; Gi=[[1/det,-c0/det],[-c0/det,1/det]]
    al=p[:2]; be=p[2:4]; E=[p[4:6],p[6:8]]; out={}
    for x in (0,1):
        for a in (1,-1):
            t=1+a*al[x]; w=[be[k]+a*E[x][k] for k in (0,1)]
            v=[(Gi[k][0]*w[0]+Gi[k][1]*w[1])/t for k in (0,1)]
            g=[arb(0)]*8; g[x]=arb(-a); g[2]=v[0]; g[3]=v[1]; g[4+2*x]=a*v[0]; g[5+2*x]=a*v[1]
            out[(a,x)]=g
    return [out[k] for k in [(1,0),(1,1),(-1,0),(-1,1)]]
def sig(g): return [g[0],-g[1],g[3],g[2],g[5],g[4],-g[7],-g[6]]
def hexagon(th):
    p,c0,(s,C,r)=ptheta(th); Gm=gammas(p,c0)
    m=[1+C,-s*s,1-C,-s*s]
    facets={'V1':[(1,1,1,1),(1,-1,-1,1)],'V2':[(1,-1,-1,1),(1,1,1,-1)],'V3':[(1,-1,1,1),(-1,1,1,-1)],'V4':[(-1,1,1,-1),(-1,-1,-1,1)]}
    V={}
    for k,(L1,L2) in facets.items():
        M=arb_mat([m,[dot(g,p) for g in Gm],[dot(g,detpt(*L1)) for g in Gm],[dot(g,detpt(*L2)) for g in Gm]])
        mu=M.solve(arb_mat([[0],[1],[1],[1]]))
        mu=[mu[i,0] for i in range(4)]
        V[k]=[sum((mu[i]*Gm[i][j] for i in range(4)),arb(0)) for j in range(8)]
    order=[V['V1'],V['V2'],sig(V['V1']),V['V3'],V['V4'],sig(V['V3'])]
    return p,order
def gtilt(th):
    p,c0,(s,C,r)=ptheta(th); al=2*C/r; n=(8+2*al*al).sqrt()
    return [al/n,arb(0),arb(0),arb(0),1/n,1/n,1/n,-1/n]
G_CHSH=[arb(v) for v in (0,0,0,0,1,1,1,-1)]
def chsh_norm(P):
    k=dot(G_CHSH,P); return [x/k for x in G_CHSH]
def Idelta(d):
    cd,sd,c2,s2d=d.cos(),d.sin(),(2*d).cos(),(2*d).sin()
    E=[arb(0),cd,cd,-s2d]; p=[arb(0)]*4+E
    IQ=2*cd**3/(c2*sd)
    g=[x/IQ for x in [arb(0)]*4+[arb(1),1/sd,1/sd,-1/c2]]
    c0=-sd; det=1-c0*c0; Gi=[[1/det,-c0/det],[-c0/det,1/det]]
    Y=[]
    for x in (0,1):
        y=[arb(0)]*8; e=E[2*x:2*x+2]
        y[2]=Gi[0][0]*e[0]+Gi[0][1]*e[1]; y[3]=Gi[1][0]*e[0]+Gi[1][1]*e[1]; y[x]=arb(-1); Y.append(y)
    W1=[-(1-cd)*c2/(2*cd**3), sd/(2*cd**3)]; W2=[-sd*c2/(2*cd**3),(cd-c2)/(2*cd**3)]
    return p,g,Y,W1,W2
def comb(g,d,Y,t):
    return [g[j]+t*(d[0]*Y[0][j]+d[1]*Y[1][j]) for j in range(8)]
def PT():
    P=[arb(0)]*4+[x/S2 for x in (arb(1),arb(1),arb(1),arb(-1))]
    gc=[x/(2*S2) for x in G_CHSH]; rho=1-1/S2
    def g_r(r0,r1): return [gc[0]+r0/S2+r1/S2, gc[1]+r0/S2-r1/S2, gc[2]-r0, gc[3]-r1]+gc[4:]
    return P,{'CHSH (centre)':gc,'half-way to g0':g_r(rho/2,0),'vertex g0':g_r(rho,0)}
def hL(g): return max((dot(g,L) for L in DET),key=lambda x:x.mid())
# ---------------- case builder (shared by runner and standalone verifier) ----------------
def fr(s): q=Fraction(s); return arb(q.numerator)/q.denominator
def build(spec):
    """spec: dict(fam=..., ...). Returns (P_target (arb list), functional g (arb list) or None)."""
    fam=spec['fam']
    if fam=='theta':
        th=theta_of(spec['theta']); p,order=hexagon(th); w=spec.get('g')
        if w is None: return p,None,order
        if w=='tilted': g=gtilt(th)
        elif w=='chsh': g=chsh_norm(p)
        elif w.startswith('vertex:'): g=order[int(w[7:])]
        elif w.startswith('lam:'):
            lam=[fr(x) for x in w[4:].split(',')]
            g=[sum((lam[i]*order[i][j] for i in range(6)),arb(0)) for j in range(8)]
        return p,g,order
    if fam=='Id':
        d={'pi/6':PI/6,'pi/12':PI/12}[spec['delta']]; p,g0,Y,W1,W2=Idelta(d); w=spec.get('g')
        if w is None: return p,None,None
        g={'centre':g0,'chsh':chsh_norm(p),'half':comb(g0,W1,Y,arb(1)/2),'W1':comb(g0,W1,Y,arb(1))}[w]
        return p,g,None
    if fam=='PT':
        p,fun=PT(); w=spec.get('g')
        return p,(fun[w] if w else None),None
