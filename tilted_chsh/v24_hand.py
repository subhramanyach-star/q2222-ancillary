"""Goal 2, axis vertices V2, V4: closed-form data, the explicit piecewise-algebraic certificate, and exact symbolic
checks of the identities used in the hand proof (valid for every theta in (0, pi/4))."""
import pickle, sympy as sp, numpy as np, sys
from arbeval import S_ as s, C_ as C, R_ as r
sys.path.insert(0, '../wbc_I_family'); from core import Phi, primal_hK
D = __import__('thm3_data').load()  # exact data from thm3_data.json
Dn = r + 1
def chk(e):  # zero modulo C^2 = 1 - s^2, r^2 = 1 + s^2  (substitute s^2 = r^2 - 1, C^2 = 2 - r^2)
    e = sp.together(e); n, d = sp.fraction(e)
    n = sp.expand(n.subs(s, sp.sqrt(r**2 - 1)))
    P = sp.Poly(n, C); out = 0
    for (i,), cf in P.terms(): out += cf * (2 - r**2)**(i // 2) * C**(i % 2)
    return sp.simplify(sp.expand(out)) == 0
for name, sg in (('V2', 1), ('V4', -1)):
    g = D['verts'][name]['g']
    a0, beta, kap, nu = g[0], g[2], g[4], g[6]
    # claimed closed forms (Cs = sg*C)
    Cs = sg * C
    cl = dict(a0=sg * (r + Cs) / (2 * Dn), beta=(Cs - 1) / (2 * Dn) * sg * (-1 if sg < 0 else 1),
              kap=(Cs + r + 2 * s**2) / (4 * Dn), nu=(r + 2 - Cs) / (4 * Dn))
    # V4: a0 = (C - r)/(2D), beta = (C+1)/(2D); V2: a0 = (r+C)/(2D), beta = (C-1)/(2D)
    cl['a0'] = (r + C) / (2 * Dn) if sg > 0 else (C - r) / (2 * Dn)
    cl['beta'] = (C - 1) / (2 * Dn) if sg > 0 else (C + 1) / (2 * Dn)
    print("==", name)
    for k, e in zip(('a0', 'beta', 'kap', 'nu'), (a0, beta, kap, nu)):
        print("   %s = %s   [closed form verified: %s]" % (k, cl[k], chk(e - cl[k])))
    print("   other entries: a1 = %s, b1-b0 = %s, c01-c00 = %s, c11+c10 = %s" % (g[1], sp.simplify(g[3] - g[2]), sp.simplify(g[5] - g[4]), sp.simplify(g[7] + g[6])))
    a0, beta, kap, nu = cl['a0'], cl['beta'], cl['kap'], cl['nu']
    u = sp.symbols('u'); u0 = 2 / r
    Q = (1 - kap * u)**2 - (beta * u + a0)**2 - nu**2 * (4 - u**2)
    A2 = kap**2 - beta**2 + nu**2
    print("   |a0| + 2 nu = 1:", chk(sg * a0 + 2 * nu - 1))
    print("   Q(u0) = 0:", chk(Q.subs(u, u0)), "  Q'(u0) = 0:", chk(sp.diff(Q, u).subs(u, u0)),
          "  Q - A2 (u-u0)^2 = 0:", all(chk(cf) for cf in sp.Poly(sp.expand(Q - A2 * (u - u0)**2), u).all_coeffs()))
    print("   forced data match: tau(c0) = -beta/2 - a0/u0 equals delta0:", chk(-beta / 2 - a0 / u0 - D['verts'][name]['d0'][0]),
          "; t* = -beta/2 - sg*kap")
    # numeric illustration of the certificate at several theta
    f = sp.lambdify((s, C, r), [a0, beta, kap, nu])
    for th in [np.pi / 16, np.pi / 8, 5 * np.pi / 32, 3 * np.pi / 16, 0.2 * np.pi]:
        sv = np.sin(2 * th); Cv = np.cos(2 * th); rv = np.sqrt(1 + sv * sv)
        A0, B, K, N = f(sv, Cv, rv); u1 = abs(A0) / K
        gg = np.array(sp.lambdify((s, C, r), D['verts'][name]['g'])(sv, Cv, rv), float)
        def tau(c):
            uu = np.sqrt(2 + 2 * c)
            return -B / 2 - A0 / uu if uu >= u1 else -B / 2 - sg * K
        cs = np.linspace(-1, 1, 4001)
        vals = np.array([Phi(gg, c, [tau(c), tau(c)]) for c in cs])
        print("   theta=%.4f: u1=%.4f (c1=%.4f) < u0=%.4f; max_c Phi_c(delta(c)) - 1 = %.2e; A2=%.4f; nu-|beta+sg*kap|=%.4f" % (
            th, u1, u1**2 / 2 - 1, 2 / rv, vals.max() - 1, K**2 - B**2 + N**2, N - abs(B + sg * K)))
