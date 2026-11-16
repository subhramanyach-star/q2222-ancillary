"""Goal 1. (a) Rigorous interval proof that for every theta in Theta the 6 points V_i(theta) are the vertices of the
section of C_theta (all 4 mu > 0 and all 14 non-incident deterministic constraints strict, for each vertex).
(b) theta*: the hexagon -> octagon transition, as a root of an explicit polynomial in s = sin 2theta."""
import pickle, sys, time
import sympy as sp, flint
from flint import arb, fmpq
from arbeval import compile_expr, theta_vars, S_, C_, R_
flint.ctx.prec = 120
D = __import__('thm3_data').load()  # exact data from thm3_data.json
syms = [S_, C_, R_]

def interval_proof(lo, hi, nbox=400, names=('V1', 'V2', 'V3', 'V4', 'V1s', 'V3s')):
    """Subdivide [lo,hi] (rational) into nbox equal theta-boxes (refining a box adaptively if needed);
    every condition must have a strictly positive lower bound on every box."""
    conds = []
    for n in names:
        V = D['verts'][n]
        for i, mu in enumerate(V['mu']): conds.append(((n, 'mu', i), compile_expr(mu, syms)))
        for key, e in V['cons'].items(): conds.append(((n,) + key, compile_expr(e, syms)))
    minlb = {}
    stack = [(lo + (hi - lo) * fmpq(k, nbox), lo + (hi - lo) * fmpq(k + 1, nbox)) for k in range(nbox)]
    nb = 0
    while stack:
        a, b = stack.pop(); nb += 1
        th = arb((a + b) / 2, (b - a) / 2)    # ball covering [a,b]
        v = theta_vars(th); ok = True
        for key, f in conds:
            val = f(v)
            if not (val > 0):
                ok = False; break
            lb = float(val.lower()); minlb[key] = min(minlb.get(key, 9), lb)
        if not ok:
            if (b - a) < fmpq(1, 10**7): return False, key, (float(a), float(b))
            m = (a + b) / 2; stack += [(a, m), (m, b)]
    return True, nb, minlb

def theta_star():
    """g_{V1}.L = 1 for L = (-1,1,-1,-1): eliminate C, r (C^2 = 1-s^2, r^2 = 1+s^2) by taking the norm."""
    e = D['verts']['V1']['cons'][('det', -1, 1, -1, -1)]
    num, den = sp.fraction(sp.together(e))
    P = sp.Poly(sp.expand(num), C_, R_)
    # reduce powers
    red = 0
    for (i, j), cf in P.terms():
        red += cf * (1 - S_**2)**(i // 2) * C_**(i % 2) * (1 + S_**2)**(j // 2) * R_**(j % 2)
    red = sp.expand(red)
    norm = 1
    conj = [red.subs({C_: sc * C_, R_: sr * R_}, simultaneous=True) for sc in (1, -1) for sr in (1, -1)]
    nm = sp.expand(conj[0] * conj[1] * conj[2] * conj[3])
    Pn = sp.Poly(nm, C_, R_)
    red2 = 0
    for (i, j), cf in Pn.terms():
        assert i % 2 == 0 and j % 2 == 0
        red2 += cf * (1 - S_**2)**(i // 2) * (1 + S_**2)**(j // 2)
    poly = sp.factor(sp.expand(red2))
    return red, poly

if __name__ == "__main__":
    lo, hi = fmpq(1, 8), fmpq(3, 16)        # theta/pi
    import math
    pi = arb.pi()
    t0 = time.time()
    # theta = pi * x; boxes in x then multiply by pi rigorously
    def run(xlo, xhi, nbox):
        global theta_vars
        orig = theta_vars
        return interval_proof(xlo, xhi, nbox)
    # wrap theta_vars to take x (theta/pi)
    import arbeval
    _tv = arbeval.theta_vars
    def tv_x(x): return _tv(x * arb.pi())
    globals()['theta_vars'] = tv_x
    ok, nb, minlb = interval_proof(lo, hi, 200)
    print("(a) Theta = [pi/8, 3pi/16]: all %d sign conditions strictly positive on every box: %s (%s boxes, %.1fs)" % (
        len(minlb) if ok else -1, ok, nb, time.time() - t0))
    if ok:
        worst = sorted(minlb.items(), key=lambda kv: kv[1])[:8]
        print("    smallest certified lower bounds:", [(k, round(v, 5)) for k, v in worst])
    red, poly = theta_star()
    print("(b) numerator of 1 - g_V1.L, L=(-1,1,-1,-1), reduced:  a + b C + c r + d C r with\n   ", red)
    print("    norm (product over C->+-C, r->+-r), factored in s:\n   ", poly)
    # locate root
    for fac, mult in sp.factor_list(poly)[1]:
        fp = sp.Poly(fac, S_)
        rts = [rt for rt in fp.real_roots() if 0 < rt.evalf() < 1]
        for rt in rts:
            sv = float(rt.evalf(30)); import math
            th = 0.5 * math.asin(sv)
            print("    factor deg %d root s=%.12f -> theta=%.12f = %.8f pi" % (fp.degree(), sv, th, th / math.pi))
    pickle.dump(dict(red=red, poly=poly), open('theta_star.pkl', 'wb'))
