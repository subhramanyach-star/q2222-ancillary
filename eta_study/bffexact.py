"""BFF certification with an EXACT rational PSD check (replaces the floating-point eigenvalue margin of
../rates_v11/bffcert.py).  Per Gauss-Radau node: the solver's dual S is PSD-repaired in the scaled frame, a margin
eps*1 is added, rounded to dyadic rationals St_q, and St_q >= 0 is verified by exact LDL^T (python-flint fmpq).
Then S_q = D^-1 St_q D^-1 with D = diag(beta^{#Z}) exact (beta dyadic), the residual r = c - <S_q,F_k> - (Az^T mu)_k is
exact, and   inf_i >= r_0 - sum_{k>0}|r_k| beta^{#Z_k} - mu.a0 + mu.z   for every quantum strategy (<S_q,Gamma> >= 0).
Only the Gauss-Radau nodes/weights (BFF Lemma 2.3) and the final summation are in double precision (error < 1e-12)."""
import sys; sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rates_v11')); sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rigorous_rates'))
import numpy as np
from fractions import Fraction
from flint import fmpq
from bffcert import BFFcert, Rate, nZ
from exact_psd import is_psd
class BFFexact(BFFcert):
    def certify(self, t, S, mode, g=None, omega=None, z=None):
        n=self.n; al=Fraction(self.beta(t)); Dq=[al**nZ(w) for w in self.S]; Df=np.array([float(d) for d in Dq])
        S=(S+S.T)/2; St=Df[:,None]*S*Df[None,:]
        lam,V=np.linalg.eigh(St); neg=lam.min(); St=(V*np.maximum(lam,0))@V.T; St=(St+St.T)/2
        sc=np.abs(St).max(); den=2**52
        for eps in (1e-13,1e-12,1e-11,1e-10,1e-9):
            E=St+eps*sc*np.eye(n)
            Q=[[Fraction(int(round(E[i,j]/sc*den)),den)*Fraction(sc) for j in range(n)] for i in range(n)]
            for i in range(n):
                for j in range(i+1,n): Q[j][i]=Q[i][j]
            ok,_=is_psd([[fmpq(x.numerator,x.denominator) for x in row] for row in Q])
            if ok: break
        if not ok: raise ValueError('no exact PSD certificate')
        Sq=[[Q[i][j]/(Dq[i]*Dq[j]) for j in range(n)] for i in range(n)]
        Sfl=np.array([[float(x) for x in row] for row in Sq])
        c=self.cvec(t); MS=self.M.T@Sfl.reshape(-1); rho=c-MS
        mask=np.ones(self.N,bool); mask[self.k0]=False
        if mode=='score':
            gh=np.asarray(g)@self.Az; lam_=max(0.0,float(np.dot(gh[mask],rho[mask])/np.dot(gh[mask],gh[mask]))); mu=lam_*np.asarray(g,float)
        else:
            mu=np.linalg.lstsq(self.Az[:,mask].T,rho[mask],rcond=None)[0]
        muf=[Fraction(x) for x in mu]; r0=None; res_pen=Fraction(0)
        for k in range(self.N):
            tot=Fraction(c[k])-sum((Sq[i][j] for (i,j) in self.pos[k]),Fraction(0))-sum((muf[q]*Fraction(self.Az[q,k]) for q in range(8) if self.Az[q,k]!=0),Fraction(0))
            if k==self.k0: r0=tot
            else: res_pen+=abs(tot)*al**nZ(self.wlist[k])
        constF=r0-res_pen-sum(muf[q]*Fraction(self.a0[q]) for q in range(8))
        const=float(constF)-1e-15*max(1.0,abs(float(constF)))
        if mode=='score': LB=const+lam_*float(Fraction(omega))-1e-15
        else: LB=const+float(sum(m*Fraction(x) for m,x in zip(muf,z)))-1e-15
        return dict(LB=LB,k0=const,kappa=mu,neg_eig=float(neg),psd_pen=0.0,res_pen=float(res_pen),eps=eps)
class RateExact(Rate):
    def __init__(self,m=8,xstar=1,extra=False):
        super().__init__(m,xstar,extra); self.B=BFFexact(xstar,extra)
if __name__=='__main__':
    import json,time
    d=[json.loads(l) for l in open('out/bff_v1.0.jsonl') if '"tilted"' in l][0]
    g=np.array([d['alpha'],0,0,0,1,1,1,-1.]); z=np.array(d['z'])
    t=time.time(); b=RateExact(8).bound('score',g=g,omega=g@z-1e-10); print(b['H_cert'],b['H_num'],d['H_cert'],time.time()-t)
