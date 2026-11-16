"""Exact rational positive-semidefiniteness test (LDL^T with symmetric pivoting) using python-flint fmpq."""

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

from flint import fmpq
def is_psd(M):
    """M: list of lists of fmpq (symmetric). Returns (True, min_pivot) if M is PSD, (False, info) otherwise.
    Exact: Gaussian elimination with diagonal pivoting; PSD iff all pivots >=0 and zero pivots have zero rows."""
    n=len(M); A=[row[:] for row in M]; act=list(range(n)); minpiv=None
    while act:
        # pick largest diagonal pivot among active (exact compare)
        k=max(act,key=lambda i:A[i][i])
        d=A[k][k]
        if d<0: return False,('negative pivot',k)
        if d==0:
            for i in act:
                if A[k][i]!=0 or A[i][k]!=0: return False,('zero pivot nonzero row',k)
            act.remove(k); continue
        minpiv=d if minpiv is None or d<minpiv else minpiv
        act.remove(k)
        rk=[A[k][j] for j in range(n)]
        for i in act:
            if A[i][k]!=0:
                f=A[i][k]/d
                for j in act: A[i][j]-=f*rk[j]
    return True,minpiv
