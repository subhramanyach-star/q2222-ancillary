"""Part 2: transparent i.i.d. confidence-region illustration for H_min = -log2 P_g (no EAT, no smoothing).
N rounds, inputs stratified (N/4 per pair xy).  Per round with inputs xy and outputs ab the 'estimator vector' is
v_u = (2a at A_x, 2b at B_y, 4ab at E_xy) (rates_v11/eat.py TO), z_hat = (1/N) sum_i v_{u_i} is unbiased for z.
 1-D (single score g):  S_hat = g.z_hat.  Hoeffding (independent rounds, ranges r_xy/N):
      P[g.z <= S_hat - t1] <= eps,  t1 = sqrt( ln(1/eps) * sum_xy r_xy^2 / (8N) ),  r_xy = max_ab g.v - min_ab g.v.
      P_g <= min_mu beta(mu g) - mu (S_hat - t1).
 8-D (full statistics): box from per-coordinate two-sided Hoeffding + union bound over the 8 coordinates (eps/8 each):
      |z_j - zhat_j| <= t_j,  t_j = sqrt(2 ln(16/eps)/N_j),  N_j = N/2 (marginals), N/4 (correlators).
      P_g <= max { P_g(z) : z in Q_NPA3, z in box }  (SDP).
 Gaussian heuristic (labelled heuristic): 1-D half-width z_eps sqrt(g'Sg), z_eps = Phi^-1(1-eps);
      8-D ellipsoid (z-zhat)' S^-1 (z-zhat) <= chi2_8(1-eps);  S = (1/4N) sum_xy Cov_xy(v) = exact covariance of z_hat.
Statistics evaluated at the honest expectation z_hat = p (typical run)."""
import sys, json; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, cvxpy as cp, fs_core as F
from scipy.optimize import minimize_scalar
from scipy.stats import norm, chi2
from npa_float import NPA
from core import z_honest
B = F.Beta('3'); NPA3 = NPA('3')
NS = [1e8, 1e9, 1e10, 1e11, 1e12]; EPS = [1e-6, 1e-10]
def vecs():
    V = {}
    for x in (0, 1):
        for y in (0, 1):
            L = []
            for a in (1, -1):
                for b in (1, -1):
                    v = np.zeros(8); v[x] = 2*a; v[2 + y] = 2*b; v[4 + 2*x + y] = 4*a*b; L.append((a, b, v))
            V[(x, y)] = L
    return V
VV = vecs()
def pab(z, x, y, a, b): return (1 + a*z[x] + b*z[2 + y] + a*b*z[4 + 2*x + y])/4
def cov(z):
    S = np.zeros((8, 8))
    for (x, y), L in VV.items():
        m = sum(pab(z, x, y, a, b)*v for a, b, v in L)
        S += sum(pab(z, x, y, a, b)*np.outer(v - m, v - m) for a, b, v in L)
    return S/4          # times 1/N
def one_d(g, z, half):
    om = g @ z - half
    r = minimize_scalar(lambda m: B(m*g) - m*om, bounds=(0, 40), method='bounded', options=dict(xatol=1e-9))
    return float(-np.log2(min(r.fun, 1.0)))
class Region:
    def __init__(s):
        cons = []; obj = 0; zt = 0; nm = 0; N = NPA3
        for e in (1, -1):
            y, G = N.block(); cons.append(G >> 0); zt = zt + N.z(y); nm = nm + y[N.k1]; obj = obj + (y[N.k1] + e*y[N.kA[1]])/2
        cons.append(nm == 1); s.zt = zt; s.cons = cons; s.obj = obj
    def box(s, lo, hi):
        p = cp.Problem(cp.Maximize(s.obj), s.cons + [s.zt >= lo, s.zt <= hi]); p.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
        return float(-np.log2(min(p.value, 1.0)))
    def ell(s, zh, Sinv, r2):
        p = cp.Problem(cp.Maximize(s.obj), s.cons + [cp.quad_form(s.zt - zh, cp.psd_wrap(Sinv)) <= r2]); p.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10)
        return float(-np.log2(min(p.value, 1.0)))
REG = Region(); A = json.load(open('asym.json')); out = []
for es in F.ETAS:
    eta = float(es); d = F.lossy_device(es); p = d['p']; g = d['g']; zc = z_honest(A[es]['chsh_x'], eta)
    for dev, z, wit in (('opt', p, g), ('chsh', zc, F.CHSH)):
        S = cov(z); Sinv = np.linalg.inv(S)
        rxy = np.array([max(wit @ v for *_, v in L) - min(wit @ v for *_, v in L) for L in VV.values()])
        for eps in EPS:
            for n in NS:
                t1 = np.sqrt(np.log(1/eps)*np.sum(rxy**2)/(8*n))
                tj = np.sqrt(2*np.log(16/eps)/(n*np.array([.5, .5, .5, .5, .25, .25, .25, .25])))
                zg = norm.isf(eps); r2 = chi2.isf(eps, 8)
                rec = dict(eta=es, device=dev, n=n, eps=eps,
                           H1_hoef=one_d(wit, z, t1), H8_box=REG.box(z - tj, z + tj),
                           H1_gauss=one_d(wit, z, zg*np.sqrt(wit @ S @ wit/n)), H8_ell=REG.ell(z, Sinv*n, r2),
                           t1=float(t1), zg2=float(zg**2), chi2_8=float(r2))
                out.append(rec); print(json.dumps(rec), flush=True)
        json.dump(out, open('out/cr.json', 'w'), indent=0)
