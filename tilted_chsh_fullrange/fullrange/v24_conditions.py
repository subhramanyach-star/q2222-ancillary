# Exact sign check (comb.sign_open: norm-root isolation + Arb) of the V2, V4 certificate conditions.
# on the open interval s in (0,1), plus Arb at s=1:  kappa < 1/2, c1 < c0, (beta + sgn(a0) kappa)^2 < nu^2, A2 > 0,
# with u1 = |a0|/kappa, c1 = u1^2/2 - 1, c0 = (1-s^2)/(1+s^2).
import sympy as sp
from comb import sign_open, at_point, S_, C_, R_
from kfield import from_sympy
from flint import arb
s,C,r=S_,C_,R_; D=r+1; c0=(1-s**2)/(1+s**2)
V2=dict(a0=(r+C)/(2*D), beta=(C-1)/(2*D), kap=(C+r+2*s**2)/(4*D), nu=(r+2-C)/(4*D), sg=1)
V4=dict(a0=(C-r)/(2*D), beta=(C+1)/(2*D), kap=(-C+r+2*s**2)/(4*D), nu=(r+2+C)/(4*D), sg=-1)
ok=True
for n,V in (('V2',V2),('V4',V4)):
    a0,b,k,nu,sg=V['a0'],V['beta'],V['kap'],V['nu'],V['sg']
    c1=(sg*a0/k)**2/2-1
    conds={'sgn(a0)': sg*a0, '1/2-kappa': sp.Rational(1,2)-k, '2kappa^2(c0+1)-a0^2 [<=> c1<c0]': 2*k**2*(c0+1)-a0**2, 'nu^2-(beta+sgn kappa)^2': nu**2-(b+sg*k)**2,
           'A2': k**2-b**2+nu**2}
    for cn,e in conds.items():
        x=from_sympy(sp.together(e))
        try: so=sign_open(x, sp.Integer(0), sp.Integer(1))
        except ValueError as er: so='FAIL '+str(er)
        s1=at_point(x, arb(1))[0]
        print(n, cn, 'sign on (0,1):', so, ' at s=1:', s1, flush=True); ok &= (so==1 and s1==1)
print('V2/V4 CONDITIONS OK on (0,1]:', ok)
