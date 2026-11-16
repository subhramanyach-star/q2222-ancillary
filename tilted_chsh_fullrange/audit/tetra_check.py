"""Audit: at theta*, p_theta* and the three tight deterministic points of V1 (++++, +--+, -+--) are affinely independent
(so the face is a genuine tetrahedron), and likewise for V1' (+++-, +-++, ----)."""
import numpy as np, mpmath as mp
mp.mp.dps=40
x=mp.findroot(lambda x: x**3-8*x**2+24*x-16, 0.91); s=mp.sqrt(x); C=mp.sqrt(1-x); r=mp.sqrt(1+x)
p=[C,0,C/r,C/r,1/r,1/r,x/r,-x/r]
def L(A0,A1,B0,B1): return [A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1]
for name,pts in [('V1',[(1,1,1,1),(1,-1,-1,1),(-1,1,-1,-1)]),('V1p',[(1,1,1,-1),(1,-1,1,1),(-1,-1,-1,-1)])]:
    M=mp.matrix([[L(*q)[i]-p[i] for i in range(8)] for q in pts])
    G=M*M.T; print(name,'Gram det of (L_i - p):',mp.nstr(mp.det(G),10))
print('theta*/pi =',mp.nstr(mp.asin(s)/2/mp.pi,15))
