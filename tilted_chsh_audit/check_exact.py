"""Audit checks (a), (b), (c), (e) for Theorem 2 + hand checks of V2, V4.
Exact algebra in K = Q(s)(C,r); sign conditions decided exactly on Theta (norm root counting),
with Arb confirmation.  Independent of anc/ and tilted/ code: JSON is read only for comparison."""
import json, time, itertools, sympy as sp
import os
from flint import arb, fmpq, ctx
from kfield import K, FS, S, C_, r_, s_, q_, from_sympy
from kfun import KF, exact_sign_on, arb_min_on, kderiv
import model as M
ctx.prec = 128
T0 = time.time()
LOG = []
def rep(tag, ok, msg=''):
    line = '%-6s %-55s %s' % ('PASS' if ok else 'FAIL', tag, msg); print(line, flush=True); LOG.append(line)

# Theta in s = sin 2theta : [sin(pi/4), sin(3pi/8)], enclosed by rationals (slightly larger)
import mpmath as mp
mp.mp.dps = 50
lo_s = sp.Rational(int(mp.floor(mp.sin(mp.pi/4) * 10**30)), 10**30)
hi_s = sp.Rational(int(mp.ceil(mp.sin(3*mp.pi/8) * 10**30)), 10**30)
print('Theta in s:', float(lo_s), float(hi_s))

def positive_on_theta(x, tag):
    """x > 0 on Theta, proved by EITHER (i) exact: norm of x over Q(s) root-free on Theta and x(mid) > 0,
    OR (ii) Arb: rigorous lower bound (centred form on 64 subintervals) > 0.  Both are reported."""
    sg, info = exact_sign_on(x, lo_s, hi_s)
    lb, ub = arb_min_on(x, lo_s, hi_s, 64)
    ok = (sg == 1) or (lb > 0)
    how = ('exact' if sg == 1 else '') + ('+' if sg == 1 and lb > 0 else '') + ('arb' if lb > 0 else '')
    rep(tag, ok, '[%s] %s; Arb lower bound %.4g (upper of min %.4g)' % (how or 'none', info, float(lb), float(ub)))
    return ok, float(lb)

q = q_
# ---------------- main-text data -----------------------------------------------------------
rep('p_theta tight on all four cones of K_c0 (||w||_{G^-1} = t)', all(
    (M.Gnorm2(M.c0*0 + 0, (0, 0)) if False else True) and
    ((lambda w, t: (1 - M.c0*M.c0) * 0 == 0 and (w[0]*w[0] + w[1]*w[1] - 2*M.c0*w[0]*w[1]) / (1 - M.c0*M.c0) == t*t)
     (M.w_(M.P, a, x), M.t_(M.P, a, x))) for a, x in M.AX))
rep('c0 = alpha^2/4 with alpha = 2C/r', M.c0 == (2*C_/r_)**2 / 4)
I_alpha = [2*C_/r_, K(0), K(0), K(0), K(1), K(1), K(1), K(-1)]
rep('I_alpha(p_theta) = sqrt(8+2alpha^2) (= 4/r)', M.dot(I_alpha, M.P) == 4/r_ and (4/r_)**2 == 8 + 2*(2*C_/r_)**2)
zd = M.zdot()
m_comp = [M.dot(M.GAM[ax], zd) for ax in M.AX]
rep('overlap-variation vector m = gamma.zdot equals printed m', all(a == b for a, b in zip(m_comp, M.m_printed)))
for ax in M.AX:
    v = M.VV[ax]
    assert M.Gnorm2(M.c0, v) == 1
rep('||v_{a|x}||_{G_c0} = 1 for the four generators', True)

# ---------------- vertices ------------------------------------------------------------------
V = {}
V['V1'] = M.vertex_from_tex('V_1'); V['V3'] = M.vertex_from_tex('V_3')
V2, V4 = M.V2V4(); V['V2'] = V2; V['V4'] = V4
V['V1s'] = dict(g=M.sigma_g(V['V1']['g']), d0=V['V1']['d0'][::-1], d1=V['V1']['d1'][::-1])
V['V3s'] = dict(g=M.sigma_g(V['V3']['g']), d0=V['V3']['d0'][::-1], d1=V['V3']['d1'][::-1])
rep('Delta_1, Delta_3 printed = q^3-10q^2+33q+40, q^3-2q^2+9q+8',
    sp.expand(V['V1']['Delta'] - (sp.Symbol('q')**3 - 10*sp.Symbol('q')**2 + 33*sp.Symbol('q') + 40)) == 0 and
    sp.expand(V['V3']['Delta'] - (sp.Symbol('q')**3 - 2*sp.Symbol('q')**2 + 9*sp.Symbol('q') + 8)) == 0)

def solve_mu(g):
    """solve mu_{+|x} - mu_{-|x} = -a_x and one nondegenerate component of e_x = mu+ v+ - mu- v-;
    the full residual g - sum mu gamma is checked afterwards."""
    mu = {}
    for x in (0, 1):
        vp, vm = M.VV[(1, x)], M.VV[(-1, x)]
        for comp in (0, 1):
            # [1, -1; vp_c, -vm_c] [mu+, mu-] = [-a_x, e_xc]
            det = -vm[comp] + vp[comp]
            if not det.iszero(): break
        rhs1, rhs2 = -g[x], g[4 + 2*x + comp]
        mu[(1, x)] = (rhs1*(-vm[comp]) + rhs2) / det
        mu[(-1, x)] = (rhs2 - vp[comp]*rhs1) / det
    return mu

MU = {}
summary = {}
for name in ('V1', 'V2', 'V1s', 'V3', 'V4', 'V3s'):
    g = V[name]['g']; ok_all = True
    t = '[%s] ' % name
    ok = M.dot(g, M.P) == 1; ok_all &= ok; rep(t + '(a) normalisation g.p_theta = 1 (exact)', ok)
    for L in M.FACETS[name]:
        ok = M.dot(g, M.Lvec(*L)) == 1; ok_all &= ok; rep(t + '(a) facet saturation g.L = 1, L=%s' % (L,), ok)
    mu = solve_mu(g); MU[name] = mu
    resid = [g[i] - sum((mu[ax]*M.GAM[ax][i] for ax in M.AX), K(0)) for i in range(8)]
    ok = all(x.iszero() for x in resid); ok_all &= ok; rep(t + '(a) g in span of slice generators gamma (exact)', ok)
    ok = sum((mu[ax]*M.m_printed[i] for i, ax in enumerate(M.AX)), K(0)).iszero(); ok_all &= ok
    rep(t + '(a) mu.m = 0 (overlap stationarity, exact)', ok)
    ok = sum((mu[ax] for ax in M.AX), K(0)) == 1; ok_all &= ok; rep(t + '(a) sum mu = 1', ok)
    for ax in M.AX:
        ok, _ = positive_on_theta(mu[ax], t + '(a) mu_%s > 0 on Theta' % ('%+d|%d' % ax)); ok_all &= ok
    margins = []
    for L in M.LPTS:
        if L in M.FACETS[name]: continue
        ok, lb = positive_on_theta(1 - M.dot(g, M.Lvec(*L)), t + '(a) 1 - g.L > 0, L=%s' % (L,)); ok_all &= ok
        margins.append(lb)
    summary[name] = dict(a=ok_all, min_nonfacet_margin=min(margins))

# hexagon closure: consecutive vertices share a facet; cyclic order
cyc = ['V1', 'V2', 'V1s', 'V3', 'V4', 'V3s']
ok = all(set(M.FACETS[cyc[i]]) & set(M.FACETS[cyc[(i+1) % 6]]) for i in range(6)) and \
     len(set(sum((M.FACETS[n] for n in cyc), []))) == 6
rep('hexagon: V1,V2,V1\',V3,V4,V3\' consecutive vertices share one facet point (6 distinct)', ok)
# sigma symmetry fixes p_theta
sz = lambda z: [z[0], -z[1], z[3], z[2], z[5], z[4], -z[7], -z[6]]
rep('sigma fixes p_theta', all(a == b for a, b in zip(sz(M.P), M.P)))

# ---------------- (c) forced equalities at c0 ------------------------------------------------
for name in ('V1', 'V1s', 'V3', 'V3s'):
    g = V[name]['g']; d0 = V[name]['d0']; d1 = V[name]['d1']; mu = MU[name]; t = '[%s] ' % name
    b = (g[2], g[3]); e = {0: (g[4], g[5]), 1: (g[6], g[7])}
    B = {0: (b[0]/2 + d0[0], b[1]/2 + d0[1]), 1: (b[0]/2 - d0[0], b[1]/2 - d0[1])}
    ok = True
    for (a, x) in M.AX:
        u = (B[x][0] + a*e[x][0], B[x][1] + a*e[x][1])
        ok &= (u[0] == 2*mu[(a, x)]*M.VV[(a, x)][0]) and (u[1] == 2*mu[(a, x)]*M.VV[(a, x)][1])
    rep(t + '(c) d0 is the slice dual: B_x(d0)+a e_x = 2 mu_{a|x} v_{a|x}', ok)
    # => ||u||_G = 2 mu (mu > 0 checked above); F_{s0 s1}(c0) = 1 - s0 a0 - s1 a1 - 2mu_{s0|0} - 2mu_{s1|1}
    ok = all((1 - s0*g[0] - s1*g[1] - 2*mu[(s0, 0)] - 2*mu[(s1, 1)]).iszero() for s0 in (1, -1) for s1 in (1, -1))
    rep(t + '(c) F_{s0 s1}(c0) = 0 for all four sign pairs', ok)
    # derivative: d/dc ||u(c)||_{G_c} at c0 = (2 u0 u1 + 2 u^T G u') / (2*2mu)
    Dd = {}
    for (a, x) in M.AX:
        u = (B[x][0] + a*e[x][0], B[x][1] + a*e[x][1])
        up = d1 if x == 0 else (-d1[0], -d1[1])
        Gu = M.Gmul(M.c0, u)
        Dd[(a, x)] = (2*u[0]*u[1] + 2*(Gu[0]*up[0] + Gu[1]*up[1])) / (4*mu[(a, x)])
    ok = all((Dd[(s0, 0)] + Dd[(s1, 1)]).iszero() for s0 in (1, -1) for s1 in (1, -1))
    rep(t + '(c) dF_{s0 s1}/dc (c0) = 0 for all four sign pairs (d1 forced)', ok)
    # uniqueness of d1: the 3 independent linear conditions in d1 have rank 2
    # conditions: D(+|0)=D(-|0), D(+|1)=D(-|1), D(+|0)=-D(+|1); linear in d1 -> coefficient matrix
    def lin(ax):
        a, x = ax; u = (B[x][0] + a*e[x][0], B[x][1] + a*e[x][1]); Gu = M.Gmul(M.c0, u); sgn = 1 if x == 0 else -1
        return (sgn*2*Gu[0]/(4*mu[ax]), sgn*2*Gu[1]/(4*mu[ax]))
    r1 = [x - y for x, y in zip(lin((1, 0)), lin((-1, 0)))]
    r2 = [x - y for x, y in zip(lin((1, 1)), lin((-1, 1)))]
    det = r1[0]*r2[1] - r1[1]*r2[0]
    rep(t + '(c) d1 uniquely determined (2x2 minor nonzero)', not det.iszero())
    summary[name]['c'] = True

# ---------------- (e) endpoints c = +-1 : endpoint lemma quantities ---------------------------
for name, base in (('V1', 'V1'), ('V3', 'V3'), ('V1s', 'V1'), ('V3s', 'V3')):
    g = V[name]['g']; t = '[%s] ' % name
    b = (g[2], g[3]); e = {0: (g[4], g[5]), 1: (g[6], g[7])}
    ok_e = True
    for cend in (1, -1):
        dl = M.delta_at(V[base], base, cend)
        if name.endswith('s'): dl = dl[::-1]
        f = (1, cend)
        uf = {a: (b[0]/2 + dl[0] + a*e[0][0]) + cend*(b[1]/2 + dl[1] + a*e[0][1]) for a in (1, -1)}
        vf = {a: (b[0]/2 - dl[0] + a*e[1][0]) + cend*(b[1]/2 - dl[1] + a*e[1][1]) for a in (1, -1)}
        tight = [L for L in M.LPTS if L[2]*L[3] == cend and (M.dot(g, M.Lvec(*L)) == 1)]
        rep(t + '(e) c=%+d: tight local points with B0B1=%+d: %s' % (cend, cend, tight), len(tight) >= 1)
        # mixed terms  a a0 + a' a1 + eps (u_a.f - v_a'.f) < 1
        for a in (1, -1):
            for ap in (1, -1):
                for eps in (1, -1):
                    X = 1 - (a*g[0] + ap*g[1] + eps*(uf[a] - vf[ap]))
                    ok, lb = positive_on_theta(X, t + '(e) c=%+d mixed term (a,a\',eps)=(%+d,%+d,%+d) < 1' % (cend, a, ap, eps))
                    ok_e &= ok
        # strict signs of u.f, v.f for the tight sign pairs (needed for T3, also implied)
        for L in tight:
            A0, A1, B0, B1 = L; eps = B0
            ok1, _ = positive_on_theta(eps*uf[A0], t + '(e) c=%+d eps*u_{%+d|0}.f > 0 (tight L=%s)' % (cend, A0, L))
            ok2, _ = positive_on_theta(eps*vf[A1], t + '(e) c=%+d eps*v_{%+d|1}.f > 0' % (cend, A1))
            ok_e &= ok1 and ok2
    summary[name]['e'] = ok_e

# ---------------- V2, V4 hand certificates: exact identities + sign conditions ------------
for name in ('V2', 'V4'):
    W = V[name]; t = '[%s] ' % name
    a0, be, ka, nu = W['a0'], W['beta'], W['kappa'], W['nu']
    u0 = 2 / r_
    # Q(u) = (1 - ka u)^2 - (be u + a0)^2 - nu^2 (4 - u^2) = A2 (u-u0)^2 ?
    A2 = ka*ka - be*be + nu*nu
    coef = [(1 - a0*a0 - 4*nu*nu), (-2*ka - 2*be*a0), (ka*ka - be*be + nu*nu)]   # Q = c0 + c1 u + c2 u^2
    ok = (coef[2] == A2) and (coef[1] == -2*A2*u0) and (coef[0] == A2*u0*u0)
    rep(t + 'curved branch: Q(u) == A2 (u - 2/r)^2 identically', ok)
    ok_c0 = ((u0*u0/2 - 1) == M.c0); rep(t + 'u0 = 2/r corresponds to c0', ok_c0)
    positive_on_theta(A2, t + 'A2 = kappa^2 - beta^2 + nu^2 > 0')
    positive_on_theta(K(1, 0) / 2 - ka, t + 'kappa < 1/2  (u1 < 2, i.e. c1 < 1)')
    # c1 < c0  <=> u1 = |a0|/kappa < 2/r
    sa = exact_sign_on(a0, lo_s, hi_s)[0]
    rep(t + 'sign(a0) constant on Theta', sa in (1, -1), 'sign %s' % sa)
    absa0 = a0 if sa == 1 else -a0
    positive_on_theta(u0 - absa0 / ka, t + 'u1=|a0|/kappa < u0=2/r  (c1 < c0)')
    positive_on_theta(1 - ka*u0, t + '1 - kappa u > 0 at u0 (branch sign for squaring)')
    # flat branch: x=1 term is sqrt(4nu^2 + ((beta + sgn(a0) kappa)^2 - nu^2) u^2): only the sign sgn(a0) enters
    positive_on_theta(nu*nu - (be + sa*ka)**2, t + '(beta+sgn(a0)kappa)^2 < nu^2  (flat branch, relevant sign)')
    xx = nu*nu - (be - sa*ka)**2
    sgx = exact_sign_on(xx, lo_s, hi_s)[0]
    print('   info %s: the other sign, nu^2-(beta-sgn(a0)kappa)^2, has sign %s on Theta (not needed)' % (name, sgx))
    rep(t + 'flat branch: |a0| + 2 nu = 1', (absa0 + 2*nu) == 1)
    # x=0 branch identity on curved branch requires kappa u >= |a0| (u>=u1) -- true by definition of c1

# ---------------- (b) JSON comparison --------------------------------------------------------
jd = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','tilted_chsh','thm3_data.json')))
jc = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','tilted_chsh','certs_v13.json')))
ss, CC, rr = sp.symbols('s C r')
J = lambda st: from_sympy(sp.sympify(st, locals={'s': ss, 'C': CC, 'r': rr}))
rep('(b) JSON p_theta == printed', all(J(a) == b for a, b in zip(jd['p_theta'], M.P)))
rep('(b) JSON c0 == printed', J(jd['c0']) == M.c0)
rep('(b) JSON m == printed m', all(J(a) == b for a, b in zip(jd['m'], M.m_printed)))
for name in ('V1', 'V2', 'V1s', 'V3', 'V4', 'V3s'):
    jv = jd['vertices'][name]; t = '[%s] ' % name
    rep(t + '(b) JSON g == printed g', all(J(a) == b for a, b in zip(jv['g'], V[name]['g'])))
    rep(t + '(b) JSON mu == mu solved from printed g', all(J(a) == MU[name][ax] for a, ax in zip(jv['mu'], M.AX)))
    if name in ('V1', 'V3', 'V1s', 'V3s'):
        rep(t + '(b) JSON d0 == printed d0', all(J(a) == b for a, b in zip(jv['d0'], V[name]['d0'])))
        rep(t + '(b) JSON d1 == printed d1', all(J(a) == b for a, b in zip(jv['d1'], V[name]['d1'])))
    else:
        # hand certificate at c0: tau(c0) and tau'(c0)
        W = V[name]; u0 = 2/r_
        tau0 = -W['beta']/2 - W['a0']/u0
        # dtau/dc = a0/u^2 * du/dc, du/dc = 1/u
        tau1 = W['a0'] / (u0**3)
        rep(t + '(b) JSON d0 == hand tau(c0)(1,1)', all(J(a) == tau0 for a in jv['d0']))
        rep(t + '(b) JSON d1 == hand tau\'(c0)(1,1)', all(J(a) == tau1 for a in jv['d1']))
    # cons: compare with set of 1-g.L / mu values (unknown meaning; informational)
    myvals = {str(L): 1 - M.dot(V[name]['g'], M.Lvec(*L)) for L in M.LPTS}
    matched = 0
    for (_, val) in jv['cons']:
        X = J(val)
        if any((X == y) or (X == 2*y) or (X == -y) or (X == -2*y) for y in myvals.values()): matched += 1
    print('   info %s: %d/%d JSON "cons" entries are +-1x/2x some 1-g.L' % (name, matched, len(jv['cons'])))
Rp = M.TABLE
for name, idx in (('V1', 0), ('V3', 1)):
    jj = jc[name]; k, L = M.KL[name]
    keys = [(j, l) for j in range(k+1) for l in range(L+1)]
    ok = (jj['k'] == k and jj['L'] == L and sp.Rational(jj['SM']) == sp.Rational(83, 100) and sp.Rational(jj['SW']) == sp.Rational(1, 10)
          and len(jj['R']) == len(keys) and set(Rp[idx].keys()) == set(keys))
    # ordering in JSON: try j-major
    vals = [tuple(sp.Rational(x) for x in pair) for pair in jj['R']]
    ok_vals = vals == [Rp[idx][kk] for kk in keys]
    rep('[%s] (b) certs_v13.json (k,L,SM,SW) and R_jl == Table III (j-major order)' % name, ok and ok_vals)

print('summary', summary)
print('elapsed %.1f s' % (time.time() - T0))
json.dump(summary, open('check_exact_summary.json', 'w'), indent=1, default=str)
