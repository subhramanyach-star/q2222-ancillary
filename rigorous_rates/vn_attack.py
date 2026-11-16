"""Explicit attacks giving rigorous UPPER bounds on the von Neumann entropy H(A_x|X=x,E).
Attack = finite mixture (label known to Eve) of two-qubit strategies with a MIXED state rho_AB purified by Eve, and
of deterministic points (H=0).  For one component with Alice measuring the rank-one projectors (1 + a A_x)/2:
   H(A_x|E) = h(p_+) + sum_a p_a S(rho_B^a) - S(rho_AB),   rho_B^a = tr_A[((1+aA_x)/2 (x) 1) rho]/p_a
(the conditional state of BE is pure).  Mixture: H(A|E K) = sum_k pi_k H_k, and H(A|E) <= H(A|E K)... is NOT the
direction we need; instead Eve HOLDS K, so the attack's entropy is H(A|E K) = sum_k pi_k H_k exactly.
Exact data: rho = O diag(lambda) O^T with O = (1-K)(1+K)^{-1} (Cayley, rational skew K -> rational orthogonal O),
rational lambda; observables cos a Z + sin a X with (cos,sin) = ((1-t^2)/(1+t^2), 2t/(1+t^2)), t rational.
So rho, observables and all correlators are EXACT rationals; only the entropies need ball arithmetic."""
import numpy as np, itertools
from fractions import Fraction as Fr
from scipy.optimize import minimize, linprog
from flint import arb
LOG2=arb(2).log()
I2=np.eye(2); Zm=np.diag([1.,-1]); Xm=np.array([[0.,1],[1,0]])
def cayley(k):
    K=np.zeros((4,4)); iu=np.triu_indices(4,1); K[iu]=k; K=K-K.T
    return (np.eye(4)-K)@np.linalg.inv(np.eye(4)+K)
def unpack(x):
    ang=x[:4]; O=cayley(x[4:10]); l=np.exp(x[10:14]-x[10:14].max()); l=l/l.sum()
    return ang,O,l
def obs(a): return np.cos(a)*Zm+np.sin(a)*Xm
def S(w): w=w[w>1e-300]; return float(-(w*np.log2(w)).sum())
def h2(p): return S(np.array([p,1-p]))
def evalf(x,gx):
    ang,O,l=unpack(x); rho=(O*l)@O.T
    A=[np.kron(obs(ang[i]),I2) for i in (0,1)]; B=[np.kron(I2,obs(ang[2+j])) for j in (0,1)]
    tr=lambda M: float(np.trace(rho@M))
    z=np.array([tr(A[0]),tr(A[1]),tr(B[0]),tr(B[1])]+[tr(A[i]@B[j]) for i in (0,1) for j in (0,1)])
    H=-S(l)
    Ax=obs(ang[gx]); 
    for a in (1,-1):
        P=np.kron((I2+a*Ax)/2,I2); M=P@rho@P; pa=np.trace(M)
        if pa<1e-15: continue
        rb=np.einsum('ijik->jk',M.reshape(2,2,2,2))/pa
        H+=pa*S(np.linalg.eigvalsh(rb)) 
    pp=float(np.trace(rho@np.kron((I2+Ax)/2,I2)))
    H+=h2(min(max(pp,0),1))
    return H,z
DETZ=[np.array([a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1],float) for a0,a1,b0,b1 in itertools.product([1,-1],repeat=4)]
def noisy_target_x(th,v):
    """parameters of the target strategy with white noise: rho = v|psi><psi| + (1-v)1/4 (Eve purifies)."""
    s=np.sin(2*th); r=np.sqrt(1+s*s); m=np.arccos(1/r)
    psi=np.array([np.cos(th),0,0,np.sin(th)])
    # orthonormal basis with psi first; O = cayley(k) cannot hit every O, so fit k by least squares
    Q,_=np.linalg.qr(np.column_stack([psi,np.eye(4)[:,1:3],np.array([np.sin(th),0,0,-np.cos(th)])]))
    if Q[:,0]@psi<0: Q[:,0]*=-1
    def res(k): return (cayley(k)[:,0]-psi)
    from scipy.optimize import least_squares
    k=least_squares(lambda k: res(k),np.zeros(6)).x
    lam=np.maximum(np.array([v+(1-v)/4]+[(1-v)/4]*3),1e-12)
    return np.concatenate([[0,np.pi/2,m,-m],k,np.log(lam)])
def colgen_vn(gx,rows_ge,seeds_x,iters=40,starts=10,rng=None,tol=1e-7,rhs_eq=None):
    rng=rng or np.random.default_rng(1)
    cols=[{'kind':'det','z':z,'H':0.0} for z in DETZ]
    for x in seeds_x:
        Hh,z=evalf(x,gx); cols.append({'kind':'q','x':x,'z':z,'H':Hh})
    hist=[]
    for it in range(iters):
        Z=np.array([c['z'] for c in cols]).T; Hs=np.array([c['H'] for c in cols])
        Aub=np.array([-(np.asarray(r)@Z) for r,w in rows_ge]) if rows_ge else None; bub=[-w for r,w in rows_ge] if rows_ge else None
        Aeq=np.ones((1,len(cols))); beq=[1.0]
        if rhs_eq is not None: Aeq=np.vstack([Aeq,Z]); beq=beq+list(rhs_eq)
        res=linprog(Hs,A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=(0,None),method='highs',
                    options={'primal_feasibility_tolerance':1e-10,'dual_feasibility_tolerance':1e-10})
        if res.status!=0: return None
        m=res.ineqlin.marginals if rows_ge else []; nu=res.eqlin.marginals[0]
        lin=-sum((mi*np.asarray(r) for mi,(r,w) in zip(m,rows_ge)),np.zeros(8))
        if rhs_eq is not None: lin=lin+res.eqlin.marginals[1:]
        # pricing: minimise H - lin.z - nu  (reduced cost)
        f=lambda x: (lambda Hz: Hz[0]-lin@Hz[1])(evalf(x,gx))
        best=[]
        x0s=[c['x'] for c in cols if c['kind']=='q']
        x0s=[x0s[i] for i in np.argsort([c['H']-lin@c['z'] for c in cols if c['kind']=='q'])[:3]]
        x0s+=[np.concatenate([rng.uniform(-np.pi,np.pi,4),rng.normal(0,1,6),rng.normal(0,2,4)]) for _ in range(starts)]
        for x0 in x0s:
            r=minimize(f,x0,method='L-BFGS-B',options={'maxiter':400})
            best.append((r.fun,r.x))
        best.sort(key=lambda t:t[0]); rc=best[0][0]-nu
        hist.append((res.fun,rc))
        if rc>-tol: break
        for val,x in best[:3]:
            if val-nu<0:
                Hh,z=evalf(x,gx); cols.append({'kind':'q','x':x,'z':z,'H':Hh})
    pi=res.x; keep=[k for k in range(len(pi)) if pi[k]>1e-13]
    return {'cols':[cols[k] for k in keep],'pi':pi[keep],'value':res.fun,'hist':hist}
# ---------------- exact / ball verification ----------------
def q(x): return Fr(float(x))
def _mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def _tr(A): return sum(A[i][i] for i in range(len(A)))
def _kron(A,B): return [[A[i//2][k//2]*B[i%2][k%2] for k in range(4)] for i in range(4)]
def freeze(c):
    """rational data of a component: tan-half-angles, skew K, lambdas."""
    if c['kind']=='det': return {'kind':'det','z':[int(round(t)) for t in c['z']]}
    x=c['x']; ang=np.mod(x[:4]+np.pi,2*np.pi)-np.pi
    t=[str(q(np.tan(a/2))) for a in ang]
    ang_,O,l=unpack(x)
    lam=[Fr(float(max(v,1e-14))).limit_denominator(10**15) for v in l]; s=sum(lam); lam=[v/s for v in lam]
    return {'kind':'q','t':t,'K':[str(q(v)) for v in x[4:10]],'lam':[str(v) for v in lam]}
from vn_exact import exact_component
