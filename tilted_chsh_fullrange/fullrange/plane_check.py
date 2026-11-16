# Exact: every vertex g (hexagon and octagon) satisfies g = sum mu_{a|x} gamma_{a|x}, mu.m = 0, g.p_theta = 1,
# with gamma_{a|x}: a_x = -a, b = v_{a|x}, e_x = a v_{a|x} (oct_data 'gens'); reduction mod r^2 = 1+s^2, C^2 = 1-s^2.
import json, sympy as sp
s,C,r=sp.symbols('s C r',positive=True); loc={'s':s,'C':C,'r':r}
J=json.load(open('oct_data.json')); T=json.load(open('thm3_data.json'))
A={**T['vertices'],**J['vertices']}
S=lambda x: sp.sympify(x,locals=loc)
def zero(e):
    n,d=sp.fraction(sp.together(e)); n=sp.expand(n)
    for _ in range(8): n=sp.expand(n.subs({r**2:1+s**2, C**2:1-s**2}))
    return n==0 or sp.simplify(n)==0
p=[S(x) for x in J['p_theta']]; m=[S(x) for x in J['m']]
gam=[]
for (a,x),v in J['gens']:
    v=[S(z) for z in v]; g=[0]*8; g[x]=-a; g[2]=v[0]; g[3]=v[1]; g[4+2*x]=a*v[0]; g[5+2*x]=a*v[1]; gam.append(g)
allok=True
for n,V in A.items():
    g=[S(z) for z in V['g']]; mu=[S(z) for z in V['mu']]
    ok1=all(zero(g[i]-sum(mu[k]*gam[k][i] for k in range(4))) for i in range(8))
    ok2=zero(sum(a*b for a,b in zip(mu,m))); ok3=zero(sum(a*b for a,b in zip(g,p))-1)
    allok&=ok1 and ok2 and ok3; print(n,'g=sum mu gamma:',ok1,' mu.m=0:',ok2,' g.p=1:',ok3)
print('ALL PLANE IDENTITIES EXACT:',allok)
