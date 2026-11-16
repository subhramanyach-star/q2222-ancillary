# Exact check: Wa_s = sigma(Wa), Wb_s = sigma(Wb) (and V1s=sigma V1, V3s=sigma V3 if present), sigma: A1->-A1, B0<->B1:
# (a0,a1,b0,b1,e00,e01,e10,e11) -> (a0,-a1,b1,b0,e01,e00,-e11,-e10). Simplified with r^2 = 1+s^2, C^2 = 1-s^2.
import json, sympy as sp
s,C,r=sp.symbols('s C r',positive=True); loc={'s':s,'C':C,'r':r}
def red(e):
    e=sp.together(sp.expand(e)); n,d=sp.fraction(e)
    rr=lambda z: sp.expand(sp.Poly(sp.expand(z),r,C).as_expr()).subs({r**2:1+s**2}) 
    for _ in range(6): n=sp.expand(n).subs({r**2:1+s**2, C**2:1-s**2}); n=sp.expand(n)
    return sp.simplify(n)==0
def sig(g): a0,a1,b0,b1,e00,e01,e10,e11=g; return [a0,-a1,b1,b0,e01,e00,-e11,-e10]
O=json.load(open('oct_data.json'))['vertices']; T=json.load(open('thm3_data.json'))['vertices']
A={**{k:v for k,v in T.items()},**O}
print('available:', list(A.keys()))
for a,b in [('Wa','Wa_s'),('Wb','Wb_s'),('Wa','Wb_s'),('Wb','Wa_s'),('V1','V1s'),('V3','V3s')]:
    if a in A and b in A:
        ga=[sp.sympify(x,locals=loc) for x in A[a]['g']]; gb=[sp.sympify(x,locals=loc) for x in A[b]['g']]
        print(a,b,'sigma image exact:', all(red(x-y) for x,y in zip(sig(ga),gb)))
