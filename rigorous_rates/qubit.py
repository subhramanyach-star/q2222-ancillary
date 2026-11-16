"""Explicit two-qubit strategies (Jordan blocks) for the 2222 scenario, in float (search) and arb (verification).
A strategy: Alice side either a qubit with A_x = cos(a_x) Z + sin(a_x) X, or a trivial block A_x = s_x*Id
(s_x = +-1); same for Bob; a real unit vector psi in C^2 (x) C^2.  Any 2222 quantum maximum is attained by such a
strategy (Jordan's lemma)."""
import numpy as np, itertools
from scipy.optimize import minimize
from flint import arb, arb_mat
from fractions import Fraction
I2=np.eye(2); Zm=np.array([[1.,0],[0,-1]]); Xm=np.array([[0.,1],[1,0]])
def obs(spec,ang):
    # spec None -> qubit observable with angle ang; spec=+-1 -> trivial
    return spec*I2 if spec is not None else np.cos(ang)*Zm+np.sin(ang)*Xm
def ops(strat):
    sa,sb,ang=strat['sa'],strat['sb'],strat['ang']
    Ax=[np.kron(obs(sa[x],ang[x]),I2) for x in (0,1)]; By=[np.kron(I2,obs(sb[y],ang[2+y])) for y in (0,1)]
    return Ax,By
def zvec_ops(Ax,By,psi):
    e=lambda M: psi@M@psi
    return np.array([e(Ax[0]),e(Ax[1]),e(By[0]),e(By[1])]+[e(Ax[x]@By[y]) for x in (0,1) for y in (0,1)])
def bell_op(c,Ax,By):
    W=c[0]*Ax[0]+c[1]*Ax[1]+c[2]*By[0]+c[3]*By[1]
    for x in (0,1):
        for y in (0,1): W=W+c[4+2*x+y]*Ax[x]@By[y]
    return W
SPECS=[None,1,-1]
def qmax(c,const=0.0,starts=12,rng=None,x0=None,return_all=False):
    """max over qubit strategies of const + c.z.  Only the qubit-qubit Jordan block and the 16 deterministic
    points need to be searched: a block where one party is trivial gives a local behaviour, whose value is at most
    the deterministic maximum.  Returns (value, strategy dict with psi)."""
    rng=rng or np.random.default_rng(0); c=np.asarray(c,float)
    best=(-np.inf,None)
    for sa in itertools.product((1,-1),repeat=2):
        for sb in itertools.product((1,-1),repeat=2):
            st={'sa':sa,'sb':sb,'ang':np.zeros(4)}
            val=-np.inf
            Ax,By=ops(st); val=np.linalg.eigvalsh(bell_op(c,Ax,By))[-1]
            if val>best[0]: best=(val,st)
    sa=sb=(None,None)
    def f(t):
        Ax,By=ops({'sa':sa,'sb':sb,'ang':t})
        return -np.linalg.eigvalsh(bell_op(c,Ax,By))[-1]
    inits=[rng.uniform(-np.pi,np.pi,4) for _ in range(starts)]
    if x0 is not None: inits=[np.asarray(x0)]+inits
    locs=[]
    for t0 in inits:
        r=minimize(f,t0,method='BFGS',options={'gtol':1e-12})
        locs.append((-r.fun,r.x))
        if -r.fun>best[0]: best=(-r.fun,{'sa':sa,'sb':sb,'ang':r.x})
    if return_all:
        out=[]
        for val,ang in sorted(locs,key=lambda t:-t[0]):
            st={'sa':sa,'sb':sb,'ang':ang}; Ax,By=ops(st); w,V=np.linalg.eigh(bell_op(c,Ax,By)); st['psi']=V[:,-1]
            out.append((val+const,st))
        return out
    val,st=best
    Ax,By=ops(st); w,V=np.linalg.eigh(bell_op(c,Ax,By)); st['psi']=V[:,-1]
    return val+const,st
def target_strategy(th):
    """|psi_theta>=cos th|00>+sin th|11>, A0=Z, A1=X, B_y=cos(m)Z+(-1)^y sin(m)X, cos m=1/r (tilted-CHSH optimal)."""
    s=np.sin(2*th); r=np.sqrt(1+s*s); m=np.arccos(1/r)
    return {'sa':(None,None),'sb':(None,None),'ang':np.array([0,np.pi/2,m,-m]),'psi':np.array([np.cos(th),0,0,np.sin(th)])}
def zvec(st):
    Ax,By=ops(st); return zvec_ops(Ax,By,st['psi'])
# ---------- rigorous evaluation ----------
def _fr(x): return Fraction(float(x))
def zvec_arb(st):
    """Exact-ball z of the strategy: angles and psi are taken as the dyadic rationals stored (psi renormalised in arb)."""
    ang=[arb(Fraction(x).numerator)/Fraction(x).denominator for x in st['ang_q']]
    psi=[arb(Fraction(x).numerator)/Fraction(x).denominator for x in st['psi_q']]
    nrm=sum((p*p for p in psi),arb(0)).sqrt(); psi=[p/nrm for p in psi]
    def ob(spec,a):
        if spec is not None: return arb_mat([[spec,0],[0,spec]])
        return arb_mat([[a.cos(),a.sin()],[a.sin(),-a.cos()]])
    I=arb_mat([[1,0],[0,1]])
    def kron(M,N):
        return arb_mat([[M[i//2,k//2]*N[i%2,k%2] for k in range(4)] for i in range(4)])
    Ax=[kron(ob(st['sa'][x],ang[x]),I) for x in (0,1)]; By=[kron(I,ob(st['sb'][y],ang[2+y])) for y in (0,1)]
    v=arb_mat([[p] for p in psi]); vt=v.transpose()
    e=lambda M:(vt*M*v)[0,0]
    return [e(Ax[0]),e(Ax[1]),e(By[0]),e(By[1])]+[e(Ax[x]*By[y]) for x in (0,1) for y in (0,1)]
def freeze(st):
    """store angles and psi as exact dyadic rationals (strings) -> reproducible certificate data."""
    return {'sa':list(st['sa']),'sb':list(st['sb']),'ang_q':[str(_fr(a)) for a in st['ang']],'psi_q':[str(_fr(p)) for p in st['psi']]}

def phiplus_strategy(angles):
    """|Phi+> with x-z-plane observables at the given angles (xi0,xi1,zeta0,zeta1); E_xy = cos(xi_x - zeta_y)."""
    return {'sa':(None,None),'sb':(None,None),'ang':np.array(angles,float),'psi':np.array([1,0,0,1])/np.sqrt(2)}
