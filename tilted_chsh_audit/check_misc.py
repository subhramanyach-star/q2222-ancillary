"""Extra audit checks: (1) face/exposure ingredients of Theorem 2 (Lemma tools(d) determinacy,
affine independence of triangles), exact; (2) sigma-symmetry of Phi for the primed vertices
(numerical spot check); (3) independent PRIMAL check (SLSQP over K_c) that h_{K_c}(g) <= Phi_c(delta(c))
and that the dual conventions B_{0,1} = b/2 +- delta are the right ones."""
import numpy as np, sympy as sp, itertools
from scipy.optimize import minimize
from flint import arb, fmpq
from kfield import K, C_, r_, s_, q_
from kfun import KF
import model as M

print('--- (1) faces / exposure ingredients (exact in K) ---')
P = M.P
uhat = {(a, x): tuple(wi / M.t_(P, a, x) for wi in M.w_(P, a, x)) for a, x in M.AX}
# unknowns z = (al0, al1, b0, b1, E00, E01, E10, E11); equations beta_y + a E_xy - (1 + a al_x) uhat_y = 0
rows = []
for a, x in M.AX:
    for y in (0, 1):
        row = [K(0)]*8
        row[x] = -a*uhat[(a, x)][y]; row[2 + y] = K(1); row[4 + 2*x + y] = K(a)
        rows.append(row)
# exact determinant by fraction-free elimination in K
Mx = [r[:] for r in rows]; det = K(1); n = 8
for col in range(n):
    piv = next((i for i in range(col, n) if not Mx[i][col].iszero()), None)
    if piv is None: det = K(0); break
    if piv != col: Mx[col], Mx[piv] = Mx[piv], Mx[col]; det = -det
    det = det * Mx[col][col]
    inv = Mx[col][col].inv()
    for i in range(col + 1, n):
        f = Mx[i][col] * inv
        if not f.iszero():
            Mx[i] = [Mx[i][j] - f*Mx[col][j] for j in range(n)]
print('Lemma tools(d): active directions determine z (8x8 system nonsingular):', not det.iszero(), ' det =', det.to_sympy())
for name in ('V1', 'V2', 'V1s', 'V3', 'V4', 'V3s'):
    L1, L2 = M.FACETS[name]
    m = (K(L1[0]) - P[0]) * (K(L2[1]) - P[1]) - (K(L1[1]) - P[1]) * (K(L2[0]) - P[0])
    print('triangle conv{p,L,L\'} for %-3s: (alpha0,alpha1)-minor =' % name, sp.simplify(m.to_sympy()))

print('--- (2)/(3) numerical spot checks ---')
def num(k, s): return float(KF(k)(arb(s)).mid())
V1 = M.vertex_from_tex('V_1'); V3 = M.vertex_from_tex('V_3')
def delta_num(V, tag, c, s):
    C = np.sqrt(1 - s*s); c0 = (1 - s*s) / (1 + s*s); sig = (s - 0.83) / 0.1
    R = M.TABLE[0] if tag == 'V1' else M.TABLE[1]
    P2 = np.zeros(2)
    for (j, l), (R0, R1) in R.items(): P2 += np.array([float(R0), float(R1)]) * c**j * sig**l
    d0 = np.array([num(V['d0'][0], s), num(V['d0'][1], s)]); d1 = np.array([num(V['d1'][0], s), num(V['d1'][1], s)])
    return d0 + d1*(c - c0) + (c - c0)**2 * P2
def Phi(gv, delta, c):
    G = np.array([[1, c], [c, 1]]); b = gv[2:4]; e = [gv[4:6], gv[6:8]]
    Bx = [b/2 + delta, b/2 - delta]; tot = 0
    for x in (0, 1):
        tot += max(a*gv[x] + np.sqrt(max(0, (Bx[x] + a*e[x]) @ G @ (Bx[x] + a*e[x]))) for a in (1, -1))
    return tot
def primal(gv, c, starts=30, rng=np.random.default_rng(0)):
    s2 = 1 - c*c
    def cons(z):
        out = []
        for a, x in M.AX:
            t = 1 + a*z[x]; p = z[2] + a*z[4+2*x]; q = z[3] + a*z[5+2*x]
            out += [s2*t*t - (q - c*p)**2 - s2*p*p, t]
        return np.array(out)
    best = -1e9
    for _ in range(starts):
        z0 = rng.uniform(-0.3, 0.3, 8)
        res = minimize(lambda z: -gv @ z, z0, constraints=[{'type': 'ineq', 'fun': cons}], method='SLSQP', options={'maxiter': 500, 'ftol': 1e-12})
        if res.success and cons(res.x).min() > -1e-9: best = max(best, -res.fun)
    return best
sg = lambda g: np.array([g[0], -g[1], g[3], g[2], g[5], g[4], -g[7], -g[6]])
worst_gap = 0; worst_sym = 0; maxPhi = -1
for tag, V in (('V1', V1), ('V3', V3)):
    for th in (np.pi/8, 5*np.pi/32, 3*np.pi/16):
        s = np.sin(2*th)
        gv = np.array([num(x, s) for x in V['g']])
        c0 = (1 - s*s) / (1 + s*s)
        for c in sorted(set(np.round(np.linspace(-0.98, 0.98, 9), 3)) | {round(c0, 6)}):
            d = delta_num(V, tag, c, s)
            ph = Phi(gv, d, c); pr = primal(gv, c)
            worst_gap = max(worst_gap, pr - ph); maxPhi = max(maxPhi, ph)
            worst_sym = max(worst_sym, abs(Phi(sg(gv), d[::-1], c) - ph))
            if abs(c - c0) < 1e-5: print('%s theta=%.4f pi c0=%.4f: primal h_Kc0 = %.10f, Phi(delta(c0)) = %.10f' % (tag, th/np.pi, c0, pr, ph))
print('max over spot grid of primal - Phi (must be <= ~1e-8):', worst_gap)
print('max Phi on grid (must be <= 1):', maxPhi)
print('max |Phi(sigma g, swap delta) - Phi(g, delta)| :', worst_sym)
