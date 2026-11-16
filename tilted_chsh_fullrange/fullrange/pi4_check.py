# theta = pi/4 (s=1): do the 8 vertices of C_theta coincide with the Tsirelson octagon O of prl_r4 Thm 3
# (g = g(Lambda=1/2, d): e = E/2, a_x = -d_x, b_y = sum_x d_x E_xy, d on the regular octagon of radius rho = 1 - 1/sqrt2)?
import json, sympy as sp
s,C,r=sp.symbols('s C r',positive=True); loc={'s':s,'C':C,'r':r}
O=json.load(open('oct_data.json'))['vertices']; T=json.load(open('thm3_data.json'))['vertices']
A={**T,**O}; sub={s:1,C:0,r:sp.sqrt(2)}
E=sp.Matrix([[1,1],[1,-1]])/sp.sqrt(2); rho=1-1/sp.sqrt(2)
ok=True
for n in ['Wa','Wb','V2','Wa_s','Wb_s','V3','V4','V3s']:
    g=[sp.nsimplify(sp.radsimp(sp.sympify(x,locals=loc).subs(sub)),[sp.sqrt(2)]) for x in A[n]['g']]
    g=[sp.radsimp(sp.simplify(x)) for x in g]
    d=[-g[0],-g[1]]
    e_ok=all(sp.simplify(g[4+2*x+y]-E[x,y]/2)==0 for x in range(2) for y in range(2))
    b_ok=all(sp.simplify(g[2+y]-(d[0]*E[0,y]+d[1]*E[1,y]))==0 for y in range(2))
    nd=sp.simplify(sp.sqrt(d[0]**2+d[1]**2)-rho)==0
    ang=sp.atan2(d[1],d[0]); angd=float(ang/sp.pi*180)
    ok&= e_ok and b_ok and nd
    print('%-5s e=E/2:%s b=dE:%s |d|=rho:%s angle(d)=%.4f deg  g=%s'%(n,e_ok,b_ok,nd,angd,[sp.nsimplify(x) for x in g]))
print('ALL MATCH O:',ok)
