"""Exact (rational) data at a rational-trigonometry qubit point.
Device: state cos th|00>+sin th|11> with (cos 2th, sin 2th) = (C,S) rational; A_x = cos a_x Z + sin a_x X, B_y likewise,
all cos/sin rational (points on the unit circle via t = tan(angle/2) rational).  Then the ideal behaviour q is rational,
Bob's overlap c0 = cos(b1-b0) and sqrt(1-c0^2) are rational, and so are the slice gradients gamma_{a|x}, the overlap
derivative zdot, the candidate polygon C(q) (slice cone + stationarity + 16 deterministic cuts), its vertices, and the
forced certificate data d0, d1 at c0.  Coordinates z=(al0,al1,be0,be1,E00,E01,E10,E11); g=(a0,a1,b0,b1,e00,e01,e10,e11)."""
from fractions import Fraction as Fr
import itertools, math, json
def circ(t, den):
    u = Fr(math.tan(t/2)).limit_denominator(den); return ((1-u*u)/(1+u*u), 2*u/(1+u*u), u)
DET = [(A0, A1, B0, B1) for A0 in (1, -1) for A1 in (1, -1) for B0 in (1, -1) for B1 in (1, -1)]
def detz(L): A0, A1, B0, B1 = L; return [A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1]
def dot(g, z): return sum((x*y for x, y in zip(g, z)), Fr(0))
def zq(C, S, ca, sa, cb, sb):
    al = [C*ca[x] for x in (0, 1)]; be = [C*cb[y] for y in (0, 1)]
    E = [ca[x]*cb[y] + S*sa[x]*sb[y] for x in (0, 1) for y in (0, 1)]
    return al + be + E
def zdot_b1(C, S, ca, sa, cb, sb):
    """d z / d b1 (rotating Bob's B1 direction at fixed state and other settings): the overlap variation."""
    d = [Fr(0)]*8; d[3] = -C*sb[1]
    for x in (0, 1): d[5+2*x] = -ca[x]*sb[1] + S*sa[x]*cb[1]
    return d
ORDER = [(1, 0), (1, 1), (-1, 0), (-1, 1)]   # (a, x)
def point(par, den=10**6):
    """par = (theta, a0, a1, b0, b1) floats -> exact rational device data."""
    th, a0, a1, b0, b1 = par
    C, S, u2 = circ(2*th, den); A = [circ(a0, den), circ(a1, den)]; B = [circ(b0, den), circ(b1, den)]
    ca = [A[0][0], A[1][0]]; sa = [A[0][1], A[1][1]]; cb = [B[0][0], B[1][0]]; sb = [B[0][1], B[1][1]]
    return dict(C=C, S=S, ca=ca, sa=sa, cb=cb, sb=sb, u=[u2, A[0][2], A[1][2], B[0][2], B[1][2]])
def gram(c): return [[Fr(1), c], [c, Fr(1)]]
def ginv_apply(c, w):
    d = 1 - c*c; return [(w[0] - c*w[1])/d, (w[1] - c*w[0])/d]
def gnorm2(u, c): return u[0]*u[0] + u[1]*u[1] + 2*c*u[0]*u[1]
def analyse(dev):
    C, S, ca, sa, cb, sb = (dev[k] for k in ('C', 'S', 'ca', 'sa', 'cb', 'sb'))
    assert C*C + S*S == 1 and all(ca[i]**2 + sa[i]**2 == 1 and cb[i]**2 + sb[i]**2 == 1 for i in (0, 1))
    z = zq(C, S, ca, sa, cb, sb); zd = zdot_b1(C, S, ca, sa, cb, sb)
    c0 = cb[0]*cb[1] + sb[0]*sb[1]; assert abs(c0) < 1
    gam = {}; t = {}; w = {}; v = {}
    for (a, x) in ORDER:
        t[a, x] = 1 + a*z[x]; w[a, x] = [z[2] + a*z[4+2*x], z[3] + a*z[5+2*x]]
        assert t[a, x] > 0
        # all four cones active (pure conditional states): w^T G^{-1} w = t^2 exactly
        gi = ginv_apply(c0, w[a, x]); assert gi[0]*w[a, x][0] + gi[1]*w[a, x][1] == t[a, x]**2, 'cone not active'
        v[a, x] = [gi[0]/t[a, x], gi[1]/t[a, x]]
        g = [Fr(0)]*8; g[x] = Fr(-a); g[2] = v[a, x][0]; g[3] = v[a, x][1]; g[4+2*x] = a*v[a, x][0]; g[5+2*x] = a*v[a, x][1]
        gam[a, x] = g
        assert dot(g, z) == 1 and gnorm2(v[a, x], c0) == 1
    m = [dot(gam[k], zd) for k in ORDER]
    G = [gam[k] for k in ORDER]
    # polygon section: mu >= 0, sum mu = 1, mu.m = 0, (sum mu_k G_k).L <= 1 for all 16 L.
    rows = []   # constraints r.mu <= rhs ; tight means r.mu == rhs
    for i in range(4): rows.append(([Fr(-1) if j == i else Fr(0) for j in range(4)], Fr(0), ('mu', i)))
    for L in DET:
        zl = detz(L); rows.append(([dot(G[k], zl) for k in range(4)], Fr(1), ('L', L)))
    verts = []
    for i, j in itertools.combinations(range(len(rows)), 2):
        M = [[Fr(1)]*4, list(m), rows[i][0], rows[j][0]]; rhs = [Fr(1), Fr(0), rows[i][1], rows[j][1]]
        mu = solve(M, rhs)
        if mu is None: continue
        if all(dot(r, mu) <= rr for r, rr, _ in rows):
            if not any(mu == vv['mu'] for vv in verts):
                g = [sum((mu[k]*G[k][jj] for k in range(4)), Fr(0)) for jj in range(8)]
                tight = [lab for r, rr, lab in rows if dot(r, mu) == rr]
                verts.append(dict(mu=mu, g=g, tight=tight))
    return dict(z=z, zd=zd, c0=c0, gam=G, m=m, v=[v[k] for k in ORDER], t=[t[k] for k in ORDER], verts=verts)
def solve(M, b):
    n = len(M); A = [list(M[i]) + [b[i]] for i in range(n)]
    for col in range(n):
        piv = next((r for r in range(col, n) if A[r][col] != 0), None)
        if piv is None: return None
        A[col], A[piv] = A[piv], A[col]
        for r in range(n):
            if r != col and A[r][col] != 0:
                f = A[r][col]/A[col][col]; A[r] = [x - f*y for x, y in zip(A[r], A[col])]
    return [A[i][n]/A[i][i] for i in range(n)]
def forced(an, V):
    """Exact certificate data at c0 for vertex V = sum mu_k gamma_k (all mu > 0):
    d0 with B_x(d0) + a e_x = 2 mu_{a|x} v_{a|x}; d1 solving the first-order equalities of all four sign pairs."""
    c0 = an['c0']; mu = dict(zip(ORDER, V['mu'])); v = dict(zip(ORDER, an['v'])); g = V['g']
    b = g[2:4]; e = [g[4:6], g[6:8]]
    d0 = [mu[1, 0]*v[1, 0][i] + mu[-1, 0]*v[-1, 0][i] - b[i]/2 for i in range(2)]
    # check: B_0 + a e_0 = 2 mu v, B_1 + a e_1 = 2 mu v
    for (a, x) in ORDER:
        B = [b[i]/2 + (d0[i] if x == 0 else -d0[i]) for i in range(2)]
        assert [B[i] + a*e[x][i] for i in range(2)] == [2*mu[a, x]*v[a, x][i] for i in range(2)]
    # d/dc ||u||_c along u(c) = u0 + s*d1*(c-c0) (s = +1 for x=0, -1 for x=1), at u0 = 2 mu v (norm 2mu):
    #   = [u0_0 u0_1 + s (G u0).d1] / ||u0|| = [ (2mu)^2 v0 v1 + s 2mu (G v).d1 ] / (2 mu) = 2 mu v0 v1 + s (G v).d1
    def lin(a, x):
        s = 1 if x == 0 else -1; vv = v[a, x]; Gv = [vv[0] + c0*vv[1], c0*vv[0] + vv[1]]
        return (2*mu[a, x]*vv[0]*vv[1], [s*Gv[0], s*Gv[1]])   # const + coef.d1
    # need D(+,0) = D(-,0) = -D(+,1) = -D(-,1)
    eqs = []
    for (k1, k2, sg) in (((1, 0), (-1, 0), 1), ((1, 1), (-1, 1), 1), ((1, 0), (1, 1), -1)):
        c1, l1 = lin(*k1); c2, l2 = lin(*k2)
        eqs.append(([l1[i] - sg*l2[i] for i in range(2)], sg*c2 - c1))
    d1 = solve([eqs[0][0], eqs[1][0]], [eqs[0][1], eqs[1][1]])
    assert d1 is not None
    r3 = dot(eqs[2][0], d1) - eqs[2][1]
    return d0, d1, r3
def lossmap(eta):
    """z' = M z + m0 (exact)."""
    M = [[Fr(0)]*8 for _ in range(8)]; m0 = [Fr(0)]*8
    for x in (0, 1): M[x][x] = eta; m0[x] = 1 - eta
    for y in (0, 1): M[2+y][2+y] = eta; m0[2+y] = 1 - eta
    for x in (0, 1):
        for y in (0, 1):
            k = 4 + 2*x + y; M[k][k] = eta*eta; M[k][x] += eta*(1 - eta); M[k][2+y] += eta*(1 - eta); m0[k] = (1 - eta)**2
    return M, m0
def push_functional(V, eta):
    """g with g.(M z + m0) = V.z + const for all z:  M^T g = V  ->  g = M^{-T} V (exact)."""
    M, _ = lossmap(eta); MT = [[M[j][i] for j in range(8)] for i in range(8)]
    return solve(MT, V)
