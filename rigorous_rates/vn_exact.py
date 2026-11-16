"""Exact (rational) / ball-arithmetic evaluation of a mixed two-qubit attack component (no solver needed).
H(A_x|E) = h(p_+) + sum_a p_a S(rho_B^a) - S(rho_AB) for Eve holding a purification of rho_AB."""
from fractions import Fraction as Fr
from flint import arb, fmpq_mat, fmpq
LOG2=arb(2).log()
def _mm(A,B): return [[sum(A[i][k]*B[k][j] for k in range(len(B))) for j in range(len(B[0]))] for i in range(len(A))]
def _tr(A): return sum(A[i][i] for i in range(len(A)))
def _kron(A,B): return [[A[i//2][k//2]*B[i%2][k%2] for k in range(4)] for i in range(4)]
def exact_component(fz,gx):
    """returns (z exact Fractions, H as arb ball)."""
    if fz['kind']=='det': return [Fr(v) for v in fz['z']],arb(0)
    t=[Fr(v) for v in fz['t']]
    cs=[((1-u*u)/(1+u*u),2*u/(1+u*u)) for u in t]
    ob=lambda c,s: [[c,s],[s,-c]]
    I=[[Fr(1),Fr(0)],[Fr(0),Fr(1)]]
    kv=[Fr(v) for v in fz['K']]; K=[[Fr(0)]*4 for _ in range(4)]; idx=0
    for i in range(4):
        for j in range(i+1,4): K[i][j]=kv[idx]; K[j][i]=-kv[idx]; idx+=1
    Mi=fmpq_mat(4,4,[fmpq((1 if i==j else 0)+K[i][j].numerator*0,1)+fmpq(K[i][j].numerator,K[i][j].denominator) for i in range(4) for j in range(4)])
    Mm=fmpq_mat(4,4,[fmpq(1 if i==j else 0,1)-fmpq(K[i][j].numerator,K[i][j].denominator) for i in range(4) for j in range(4)])
    Of=Mm*Mi.inv()
    O=[[Fr(int(Of[i,j].p),int(Of[i,j].q)) for j in range(4)] for i in range(4)]
    lam=[Fr(v) for v in fz['lam']]; assert all(v>0 for v in lam) and sum(lam)==1
    # check orthogonality exactly
    OtO=_mm([list(r) for r in zip(*O)],O); assert all(OtO[i][j]==(1 if i==j else 0) for i in range(4) for j in range(4))
    rho=[[sum(O[i][k]*lam[k]*O[j][k] for k in range(4)) for j in range(4)] for i in range(4)]
    A=[_kron(ob(*cs[i]),I) for i in (0,1)]; B=[_kron(I,ob(*cs[2+j])) for j in (0,1)]
    tr=lambda M:_tr(_mm(rho,M))
    z=[tr(A[0]),tr(A[1]),tr(B[0]),tr(B[1])]+[tr(_mm(A[i],B[j])) for i in (0,1) for j in (0,1)]
    def ent(ws):
        out=arb(0)
        for w in ws:
            if w.mid()>0: out-=w*w.log()/LOG2
        return out
    H=-ent([arb(v.numerator)/v.denominator for v in lam])
    c,s=cs[gx]; Ax=ob(c,s)
    for a in (1,-1):
        Pa=_kron([[ (I[i][j]+a*Ax[i][j])/2 for j in range(2)] for i in range(2)],I)
        M=_mm(_mm(Pa,rho),Pa); pa=_tr(M)
        if pa==0: continue
        rb=[[sum(M[2*k+i][2*k+j] for k in range(2))/pa for j in range(2)] for i in range(2)]
        det=rb[0][0]*rb[1][1]-rb[0][1]*rb[1][0]
        d=arb(1-4*det.numerator/Fr(det.denominator)) if False else (arb((1-4*det).numerator)/(1-4*det).denominator)
        sq=d.sqrt() if d>0 else arb(0)
        e1=(1+sq)/2; e2=(1-sq)/2
        H+=arb(pa.numerator)/pa.denominator*ent([e1,e2])
    pp=_tr(_mm(rho,_kron([[ (I[i][j]+Ax[i][j])/2 for j in range(2)] for i in range(2)],I)))
    H+=ent([arb(pp.numerator)/pp.denominator,arb((1-pp).numerator)/(1-pp).denominator])
    return z,H
