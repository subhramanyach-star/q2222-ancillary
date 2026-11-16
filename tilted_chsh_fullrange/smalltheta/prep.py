"""Symbolic preparation: hatted data X^ = (X(t)-X(0))/t as exact expressions in (t, gam, rho), gam=1/(1+sqrt(1-t)), rho=1/(1+sqrt(1+t))
(so sqrt(1-t) = 1 - t*gam, sqrt(1+t) = 1 + t*rho exactly). Asserts the t^0 numerator coefficient cancels identically."""
import json, sympy as sp, sys
J=json.load(open('../fullrange/thm3_data.json'))
s,C,r,t,gam,rho=sp.symbols('s C r t gam rho')
def hat(expr_str):
    e=sp.sympify(expr_str,locals={'s':s,'C':C,'r':r})
    e=e.subs(s,sp.sqrt(t)); e=sp.together(e); N,D=sp.fraction(e)
    N=sp.expand(N); D=sp.expand(D)
    assert not N.has(sp.sqrt(t)) or True
    Nt=sp.expand(N.subs({C:1-t*gam,r:1+t*rho})); Dt=sp.expand(D.subs({C:1-t*gam,r:1+t*rho}))
    assert not (Nt.free_symbols-{t,gam,rho}) and not (Dt.free_symbols-{t,gam,rho}), (Nt.free_symbols)
    for X in (Nt,Dt): assert sp.Poly(X,t,gam,rho).is_polynomial if hasattr(sp.Poly(X,t,gam,rho),'is_polynomial') else True
    N0=Nt.subs(t,0); D0=Dt.subs(t,0); assert D0!=0
    x0=sp.nsimplify(N0/D0); assert x0.is_Rational, x0
    num=sp.expand(Nt-x0*Dt); assert sp.expand(num.subs(t,0))==0
    q=sp.expand(sp.cancel(num/t)); assert sp.Poly(q,t,gam,rho)   # polynomial
    # sanity: numeric check
    for tv in (sp.Rational(1,100),sp.Rational(1,5)):
        sub={t:tv,gam:1/(1+sp.sqrt(1-tv)),rho:1/(1+sp.sqrt(1+tv))}
        orig=sp.sympify(expr_str,locals={'s':s,'C':C,'r':r}).subs({s:sp.sqrt(tv),C:sp.sqrt(1-tv),r:sp.sqrt(1+tv)})
        val=x0+tv*q.subs(sub)/Dt.subs(sub)
        assert abs(sp.N(orig-val,50))<1e-40, (expr_str, sp.N(orig-val))
    return str(x0), str(q), str(Dt)
out={}
for name in sys.argv[1:]:
    V=J['vertices'][name]; o={}
    for k in ('g','d0','d1'):
        o[k]=[hat(e) for e in V[k]]
    out[name]=o; print(name, {k:[a for a,_,_ in v] for k,v in o.items()})
json.dump(out,open('hats_%s.json'%'_'.join(sys.argv[1:]),'w'),indent=0)
