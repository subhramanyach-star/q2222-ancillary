"""Rigorous (Arb) version of tetra_check.py: Gram determinant of (L_i - p_theta*) for V1's three tight points at theta*,
x* = sin^2(2 theta*) the real root of x^3-8x^2+24x-16 in [0.9, 0.95] (sign change certified, then bisection)."""
from flint import arb, fmpq, arb_mat
f=lambda x: x**3-8*x**2+24*x-16
a,b=fmpq(9,10),fmpq(95,100); assert f(arb(a))<0 and f(arb(b))>0
for _ in range(80):
    m=(a+b)/2
    if f(arb(m))<0: a=m
    else: b=m
x=arb((a+b)/2,(b-a)/2); C=(1-x).sqrt(); r=(1+x).sqrt()
p=[C,arb(0),C/r,C/r,1/r,1/r,x/r,-x/r]
def L(A0,A1,B0,B1): return [A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1]
for name,pts in [('V1',[(1,1,1,1),(1,-1,-1,1),(-1,1,-1,-1)]),('V1p',[(1,1,1,-1),(1,-1,1,1),(-1,-1,-1,-1)])]:
    D=[[arb(L(*q)[i])-p[i] for i in range(8)] for q in pts]
    G=arb_mat([[sum((D[i][k]*D[j][k] for k in range(8)),arb(0)) for j in range(3)] for i in range(3)])
    d=G.det(); print(name,'Gram det =',d,' >0:',d>0)
