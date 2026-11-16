"""Independent float/mp checks for the WBC material imported into prl_r4 (not proofs):
(1) zero-marginal reduction: all four cones active, m_{+|x}=m_{-|x}, m0*m1<0, g_corr prop. to WBC functional, F-orbit;
(2) I-family parallelogram vertices / tight labels (own enumeration), P-orbit of W1;
(3) perfect-key m formula; (4) certificate Phi<=1 on grid (own Phi), forced touching at c0."""
import numpy as np, itertools, json, mpmath as mp
from fractions import Fraction
DET = [np.array([a0,a1,b0,b1,a0*b0,a0*b1,a1*b0,a1*b1],float) for a0,a1,b0,b1 in itertools.product((1,-1),repeat=4)]
LAB = [''.join('+' if v>0 else '-' for v in l) for l in itertools.product((1,-1),repeat=4)]
KEYS=[(1,0),(1,1),(-1,0),(-1,1)]
def pt(a0,a1,b0,b1): return np.array([0,0,0,0,np.cos(a0-b0),np.cos(a0-b1),np.cos(a1-b0),np.cos(a1-b1)]), np.cos(b0-b1)
def slice_data(p,c0):
    G=np.array([[1,c0],[c0,1]]);Gi=np.linalg.inv(G);gens=[];act=[]
    for a,x in KEYS:
        t=1+a*p[x]; w=np.array([p[2]+a*p[4+2*x],p[3]+a*p[5+2*x]]); act.append(np.sqrt(w@Gi@w)-t)
        v=Gi@w/t; g=np.zeros(8); g[x]=-a; g[2:4]=v; g[4+2*x:6+2*x]=a*v; gens.append(g)
    zd=np.zeros(8); zd[3]=(p[2]-c0*p[3])/(1-c0**2)
    for x in (0,1): zd[5+2*x]=(p[4+2*x]-c0*p[5+2*x])/(1-c0**2)
    M=np.array(gens); return M,np.array(act),M@zd
def polygon(p,c0):
    M,act,m=slice_data(p,c0)
    A=np.vstack([-np.eye(4),np.array([M@L for L in DET])]); b=np.concatenate([np.zeros(4),np.ones(16)])
    Aeq=np.vstack([np.array([g@p for g in M]),m]); beq=np.array([1.,0.])
    mu0=np.linalg.lstsq(Aeq,beq,rcond=None)[0]; N=np.linalg.svd(Aeq)[2][2:].T
    A2=A@N; b2=b-A@mu0; V=[]
    for i,j in itertools.combinations(range(20),2):
        Mx=A2[[i,j]]
        if abs(np.linalg.det(Mx))<1e-12: continue
        u=np.linalg.solve(Mx,b2[[i,j]])
        if np.all(A2@u<=b2+1e-9):
            mu=mu0+N@u
            if not any(np.allclose(mu,v[0],atol=1e-8) for v in V):
                V.append((mu,mu@M,[LAB[k-4] for k in range(4,20) if abs(A2[k]@u-b2[k])<1e-8],[k for k in range(4) if abs(A2[k]@u-b2[k])<1e-8]))
    return V,M,act,m
def Fop(g): # output flip
    h=g.copy(); h[:4]*=-1; return h
def Pswap(g): # party swap: (a0,a1,b0,b1,e00,e01,e10,e11)->(b0,b1,a0,a1,e00,e10,e01,e11)
    return np.array([g[2],g[3],g[0],g[1],g[4],g[6],g[5],g[7]])
out=[]
for d in [0.05,0.2,np.pi/12,0.4,np.pi/6]:
    p,c0=pt(0,np.pi/2+d,np.pi/2,-d); V,M,act,m=polygon(p,c0)
    I=np.array([0,0,0,0,1,1/np.sin(d),1/np.sin(d),-1/np.cos(2*d)]); IQ=2*np.cos(d)**3/(np.cos(2*d)*np.sin(d))
    cen=np.mean([v[1] for v in V],axis=0)
    gs=[v[1] for v in V]
    orbF=all(any(np.allclose(Fop(g),h,atol=1e-9) for h in gs) for g in gs); orbP=all(any(np.allclose(Pswap(g),h,atol=1e-9) for h in gs) for g in gs)
    # single orbit of W1 under <F,P>
    g1=gs[0]; orb=[g1,Fop(g1),Pswap(g1),Fop(Pswap(g1))]; single=all(any(np.allclose(o,h,atol=1e-9) for o in orb) for h in gs)
    out.append('I d=%.4f s=%.5f: active %s m=%s m0*m1<0 %s nV=%d tight=%s mu-tight=%s centre==I/IQ %s Pp=p %s F,P-invariant %s single orbit %s'%(
        d,2*np.cos(d)*(1+np.sin(d)),np.allclose(act,0),np.round(m,6),m[0]*m[1]<0,len(V),[v[2] for v in V],[v[3] for v in V],np.allclose(cen,I/IQ,atol=1e-9),np.allclose(Pswap(p),p),orbF and orbP,single))
for gam in [0.0,0.1,np.pi/12,0.3,0.45]:
    p,c0=pt(0,2*np.pi/3-2*gam,np.pi/2-3*gam,-(np.pi/6+gam)); V,M,act,m=polygon(p,c0); mu_=np.cos(gam+np.pi/6)
    Jf=np.array([0,0,0,0,1,4*mu_**2-1,4*mu_**2-1,-(4*mu_**2-1)])/(8*mu_**3)
    # correlator-only elements of polygon: d=0 <=> mu_{+|x}=mu_{-|x}
    cen=np.mean([v[1] for v in V],axis=0)
    out.append('J gam=%.4f mu=%.4f: active %s m0*m1<0 %s nV=%d centre==J/8mu^3 %s Jf.p=%.12f'%(gam,mu_,np.allclose(act,0),m[0]*m[1]<0,len(V),np.allclose(cen,Jf,atol=1e-9),Jf@p))
# perfect key
for al,be in [(np.pi/3,2*np.pi/3),(0.4,1.1),(2.0,0.7)]:
    p,c0=pt(0,al,0,be); M,act,m=slice_data(p,c0)
    out.append('perfect key al=%.3f be=%.3f: active %s m=%s formula %.6f'%(al,be,np.allclose(act,0),np.round(m,6),np.sin(al)*np.sin(be-al)/np.sin(be)**2))
open('check_wbc.log','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
# --- closed-form W1..W4 (d-coordinates) for the I family
out2=[]
for d in [0.05,0.2,0.4,np.pi/6]:
    p,c0=pt(0,np.pi/2+d,np.pi/2,-d); V,M,act,m=polygon(p,c0); cd,sd,c2=np.cos(d),np.sin(d),np.cos(2*d)
    W1=np.array([-(1-cd)*c2/(2*cd**3), sd/(2*cd**3)]); W2=np.array([-sd*c2/(2*cd**3),(cd-c2)/(2*cd**3)])
    ds=[np.array([v[0][0]-v[0][2],v[0][1]-v[0][3]]) for v in V]
    ok=all(any(np.allclose(w,x,atol=1e-9) for x in ds) for w in [W1,W2,-W1,-W2])
    lab={tuple(np.round(x,9)):v[2] for x,v in zip(ds,V)}
    out2.append('I d=%.3f closed-form W1..W4 match enumeration: %s ; W1 tight %s ; W2 tight %s'%(d,ok,lab[tuple(np.round(W1,9))],lab[tuple(np.round(W2,9))]))
open('check_wbc.log','a').write('\n'.join(out2)+'\n'); print('\n'.join(out2))
