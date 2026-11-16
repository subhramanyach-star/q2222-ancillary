"""Common tools for the DI randomness/key normal-cone workstream (internal working name).
Coordinates z=(alpha0,alpha1,beta0,beta1,E00,E01,E10,E11) as in ../core.py and prl.tex.
Qubit strategies with the maximally entangled state and xz-plane measurements:
observable at angle phi is cos(phi) Z + sin(phi) X, and E_xy = cos(a_x - b_y), marginals 0."""
import sys, os, itertools
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # ancillary copy: core.py is in this directory
import numpy as np
import cvxpy as cp
from core import Phi, dual_hK, gram_norm

def det_points():
    pts = []
    for A0, A1, B0, B1 in itertools.product((1, -1), repeat=4):
        pts.append(np.array([A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1], float))
    return np.array(pts)
DET = det_points()
DETLAB = list(itertools.product((1, -1), repeat=4))

def point_maxent(a0, a1, b0, b1):
    E = [np.cos(a0-b0), np.cos(a0-b1), np.cos(a1-b0), np.cos(a1-b1)]
    return np.array([0, 0, 0, 0] + E, float), np.cos(b0-b1)

# ---- anchor families (WBC 2205.00124) ----
def I_delta_angles(d):   # A0=Z, B0=X, A1=-sin d Z+cos d X, B1=cos d Z - sin d X
    return 0.0, np.pi/2 + d, np.pi/2, -d
def J_gamma_angles(g):
    return 0.0, 2*np.pi/3 - 2*g, np.pi/2 - 3*g, -(np.pi/6 + g)
def I_delta_functional(d):
    return np.array([0, 0, 0, 0, 1, 1/np.sin(d), 1/np.sin(d), -1/np.cos(2*d)])
def J_gamma_functional(g):
    cg = 4*np.cos(g + np.pi/6)**2 - 1
    return np.array([0, 0, 0, 0, 1, cg, cg, -cg])
CHSH = np.array([0, 0, 0, 0, 1, 1, 1, -1.])

def tw(z, a, x):
    t = 1 + a*z[x]
    w = np.array([z[2] + a*z[4+2*x], z[3] + a*z[5+2*x]])
    return t, w

KEYS = [(1, 0), (1, 1), (-1, 0), (-1, 1)]   # order (+|0),(+|1),(-|0),(-|1)

def slice_data(p, c0):
    """generators gamma_{a|x} of N_{K_c0}(p), activity residuals, overlap-variation vector m."""
    G = np.array([[1, c0], [c0, 1]]); Gi = np.linalg.inv(G)
    gens = []; act = []
    for (a, x) in KEYS:
        t, w = tw(p, a, x)
        act.append(np.sqrt(w @ Gi @ w) - t)
        v = Gi @ w / t
        g = np.zeros(8); g[x] = -a; g[2:4] = v; g[4+2*x:6+2*x] = a*v
        gens.append(g)
    M = np.array(gens)
    zd = np.zeros(8)
    zd[3] = (p[2] - c0*p[3])/(1 - c0**2)
    for x in (0, 1):
        zd[5+2*x] = (p[4+2*x] - c0*p[5+2*x])/(1 - c0**2)
    m = M @ zd
    return M, np.array(act), m, zd

def section_polygon(p, c0, tol=1e-12):
    """C = {sum mu gamma : mu>=0, sum mu = 1, mu.m = 0, g.L <= 1 for all L} as polygon in mu-space.
    Returns vertices (mu, g, tight constraints)."""
    M, act, m, zd = slice_data(p, c0)
    # constraints in mu (4 vars): mu>=0 (4), L: (M L) . mu <= 1 (16); equalities sum mu=1, m.mu=0
    A_ineq = np.vstack([-np.eye(4), (M @ DET.T).T]); b_ineq = np.concatenate([np.zeros(4), np.ones(16)])
    Aeq = np.vstack([np.ones(4), m]); beq = np.array([1.0, 0.0])
    # parametrise affine plane
    mu0 = np.linalg.lstsq(Aeq, beq, rcond=None)[0]
    N = np.linalg.svd(Aeq)[2][2:].T   # 4x2 null space
    A2 = A_ineq @ N; b2 = b_ineq - A_ineq @ mu0
    verts = []
    for i, j in itertools.combinations(range(len(b2)), 2):
        Mx = A2[[i, j]]
        if abs(np.linalg.det(Mx)) < 1e-12: continue
        u = np.linalg.solve(Mx, b2[[i, j]])
        if np.all(A2 @ u <= b2 + 1e-9):
            mu = mu0 + N @ u
            if not any(np.allclose(mu, v[0], atol=1e-8) for v in verts):
                tight = [k for k in range(len(b2)) if abs(A2[k] @ u - b2[k]) < 1e-8]
                verts.append((mu, mu @ M, tight))
    # order by angle
    if verts:
        U = np.array([np.linalg.lstsq(N, v[0]-mu0, rcond=None)[0] for v in verts]); cen = U.mean(0)
        ang = np.arctan2(U[:, 1]-cen[1], U[:, 0]-cen[0]); order = np.argsort(ang)
        verts = [verts[k] for k in order]
    return dict(M=M, act=act, m=m, zd=zd, verts=verts, mu0=mu0, N=N)

def tight_label(k):
    if k < 4: return 'mu%s|%d=0' % ('+' if KEYS[k][0] > 0 else '-', KEYS[k][1])
    L = DETLAB[k-4]; return 'L(' + ''.join('+' if s > 0 else '-' for s in L) + ')'

def hK_profile(g, cs):
    return np.array([dual_hK(g, c) for c in cs])

def cgrid(c0, n=201):
    g = list(np.linspace(-1, 1, n)) + [c0]
    for k in range(14):
        e = 0.05*0.6**k; g += [c0-e, c0+e]
    for e in (1e-2, 1e-3, 1e-4):
        g += [1-e, -1+e]
    return np.array(sorted(set(np.round([x for x in g if -1 <= x <= 1], 14))))
