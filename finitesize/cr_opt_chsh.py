"""Addendum to cr.py: 1-D Hoeffding/Gaussian with the OPTIMAL single witness for the CHSH device (f* from full-stats dual)."""
import sys, json; sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import numpy as np, fs_core as F
from scipy.optimize import minimize_scalar
from scipy.stats import norm
from npa_float import pguess_full_dual
from core import z_honest
B = F.Beta('3'); A = json.load(open('asym.json'))
def vecs():
    V = {}
    for x in (0, 1):
        for y in (0, 1):
            L = []
            for a in (1, -1):
                for b in (1, -1):
                    v = np.zeros(8); v[x] = 2*a; v[2 + y] = 2*b; v[4 + 2*x + y] = 4*a*b; L.append(v)
            V[(x, y)] = L
    return V
VV = vecs()
def pab(z, x, y, a, b): return (1 + a*z[x] + b*z[2 + y] + a*b*z[4 + 2*x + y])/4
def cov(z):
    S = np.zeros((8, 8))
    for (x, y), L in VV.items():
        ab = [(a, b) for a in (1, -1) for b in (1, -1)]
        m = sum(pab(z, x, y, a, b)*v for (a, b), v in zip(ab, L))
        S += sum(pab(z, x, y, a, b)*np.outer(v - m, v - m) for (a, b), v in zip(ab, L))
    return S/4
def one_d(g, z, half):
    om = g @ z - half
    r = minimize_scalar(lambda m: B(m*g) - m*om, bounds=(0, 40), method='bounded', options=dict(xatol=1e-9))
    return float(-np.log2(min(r.fun, 1.0)))
out = []
for es in F.ETAS:
    zc = z_honest(A[es]['chsh_x'], float(es)); _, _, f = pguess_full_dual('3', ('A', 1), zc); f = np.array(f)
    f = max((f, -f), key=lambda g: one_d(g, zc, 0.0)); f = f/np.linalg.norm(f)
    S = cov(zc); rxy = np.array([max(f @ v for v in L) - min(f @ v for v in L) for L in VV.values()])
    for eps in (1e-6, 1e-10):
        for n in (1e8, 1e9, 1e10, 1e11, 1e12):
            t1 = np.sqrt(np.log(1/eps)*np.sum(rxy**2)/(8*n))
            rec = dict(eta=es, n=n, eps=eps, H1_hoef_fopt=one_d(f, zc, t1), H1_gauss_fopt=one_d(f, zc, norm.isf(eps)*np.sqrt(f @ S @ f/n)), H_asym=one_d(f, zc, 0.0))
            out.append(rec); print(json.dumps(rec), flush=True)
json.dump(out, open('out/cr_opt_chsh.json', 'w'), indent=0)
