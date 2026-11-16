"""Common pieces for the finite-size comparison (finitesize/).
Coordinates z=(a0,a1,b0,b1,E00,E01,E10,E11). Guess of Alice's A_1 (x*=1), as in lossycone.
Affine upper bounds on Eve's guessing probability valid for ALL quantum z (rigorous_rates/rates_core.py, Lemma 2):
    P_g(z) <= beta(f) - f.z ,   beta(f) = max_e h_Q(Pi_e + f),  Pi_e = (1 + e <A_1>)/2,  any f in R^8.
(single score g.z >= omega is the special case f = mu g, mu >= 0: Lemma 1).  Then for any s > 0 (tangent of -log2):
    H(A|E) >= H_min(A|E) >= -log2(beta - f.z) >= K0 + kappa.z,  kappa = f/(s ln2),  K0 = -log2 s - (beta - s)/(s ln2).
"""
import sys, json, warnings; warnings.simplefilter('ignore')
import os as _os; _H = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
for p in [_os.path.join(_H, d) for d in ('eta_study', 'rigorous_rates', 'lossycone', 'rates_v11')]:
    if p not in sys.path: sys.path.insert(0, p)
import numpy as np, cvxpy as cp
from fractions import Fraction as Fr
from npa_float import NPA
LN2 = np.log(2.0)
CHSH = np.array([0, 0, 0, 0, 1, 1, 1, -1.0])
ETAS = ['0.800', '0.850', '0.900']
# ---------------------------------------------------------------- devices
def lossy_device(eta_s):
    """honest lossy behaviour p (float) and certified optimal witness g (float, max|g|=1) of lossycone T4."""
    import exact as EX
    st = json.load(open(_os.path.join(_H, 'lossycone', 'data', 't4', '%s.json' % eta_s)))
    e = Fr(eta_s); q = [Fr(x) for x in st['point']['z']]; M, m0 = EX.lossmap(e)
    p = [sum((M[i][j]*q[j] for j in range(8)), Fr(0)) + m0[i] for i in range(8)]
    g = np.array([float(Fr(x)) for x in st['g']]); p = np.array([float(x) for x in p])
    return dict(p=p, g=g, omega=float(g @ p), H_lo=st['H_lo'], H_up=st['H_up'], H_cap=st['H_cap'], S_up=float(Fr(st['S_up'])))
def chsh_device(eta):
    """qubit device maximising the lossy CHSH value at efficiency eta (no-click -> +1)."""
    from core import z_honest
    from scipy.optimize import minimize
    best = None
    rng = np.random.default_rng(1)
    for k in range(40):
        x0 = rng.uniform(-np.pi, np.pi, 5); x0[0] = rng.uniform(0.05, np.pi/4)
        r = minimize(lambda x: -CHSH @ z_honest(x, eta), x0, method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-12, fatol=1e-14))
        if best is None or r.fun < best.fun: best = r
    p = z_honest(best.x, eta)
    return dict(p=p, g=CHSH.copy(), omega=float(CHSH @ p), x=list(best.x))
# ---------------------------------------------------------------- beta(f) via NPA (float; valid relaxation)
class Beta:
    """beta(f) = max_e max_{z in Q_NPA(level)} (1/2)(1 + e z_A1) + f.z  (upper bound on max_e h_Q(Pi_e+f))."""
    def __init__(s, level='2'):
        N = NPA(level); y, G = N.block(); s.c = cp.Parameter(8); z = N.z(y)
        s.prob = cp.Problem(cp.Maximize(s.c @ z), [G >> 0, y[N.k1] == 1]); s.cache = {}
    def hq(s, c):
        s.c.value = np.asarray(c, float)
        try: s.prob.solve(solver='CLARABEL', tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=300)
        except Exception: s.prob.solve(solver='SCS', eps=1e-9)
        return float(s.prob.value)
    def __call__(s, f):
        f = np.asarray(f, float); key = tuple(np.round(f, 13))
        if key in s.cache: return s.cache[key]
        v = max(0.5 + s.hq(f + e*0.5*np.eye(8)[1]) for e in (1, -1)); s.cache[key] = v; return v
def tangent(beta, f, s):
    f = np.asarray(f, float); return -np.log2(s) - (beta - s)/(s*LN2), f/(s*LN2)
# ---------------------------------------------------------------- EAT (vectorised copy of rates_v11/eat.py)
import eat as EAT0
TOE = np.array([e for (*_, e) in EAT0.TO])                # 16 x 8 estimator vectors e_u, s(u) = kappa.e_u
QNS = np.array([EAT0.q_of_z(v) for v in EAT0.NSV])        # 24 x 16
GAMMAS = np.logspace(-5, -0.3, 48); ALPHAS = 1 + np.logspace(-10, -0.5, 96); CFRAC = np.linspace(-2.0, 1.5, 36)
def eat_rate(K0, kappa, z_hon, n, eps_h, eps_eat, eps_c=1e-6, gammas=GAMMAS, alphas=ALPHAS, cfrac=CFRAC, full=False):
    """Identical formula to rates_v11/eat.py:eat_rate (Liu et al. Nat Phys 17, 448 (2021), App. A Thm 3 + spot checking),
    vectorised over (gamma, c_perp, alpha).  Returns best rate (bits per round, H_min^eps / n)."""
    kappa = np.asarray(kappa, float); su = TOE @ kappa; q_hon = EAT0.q_of_z(z_hon); h_hon = K0 + kappa @ z_hon
    if h_hon <= 0: return (0.0, {}) if full else 0.0
    log_err = np.log2(1/eps_eat) + np.log2(1/(eps_h**2/(1 + np.sqrt(1 - eps_h**2))))
    G = gammas[:, None, None]; C = cfrac[None, :, None]
    gu = (1 - G)*(K0 + su[None, None, :])                 # 48x1x16
    Maxg = gu.max(-1, keepdims=True); Ming = (QNS @ gu[:, 0, :].T).min(0)[:, None, None]   # 48x1x1
    c = Ming + C*(Maxg - Ming)                            # 48x36x1
    dev2 = (c - gu)**2                                    # 48x36x16
    Var = np.einsum('ku,gcu->gck', QNS, dev2).max(-1)[..., None]/G   # 48x36x1
    V = LN2/2*(EAT0.LOG2_9 + np.sqrt(2 + Var))**2
    Maxf = np.maximum(c, c + (Maxg - c)/G); D = 2 + Maxf - Ming
    fu = gu/G + (1 - 1/G)*c                               # 48x36x16
    pr_t = G*q_hon[None, None, :]; pr_g = (1 - G)          # weights (sum to 1)
    EF = pr_g*c + (pr_t*fu).sum(-1, keepdims=True)
    sig2 = pr_g*(c - EF)**2 + (pr_t*(fu - EF)**2).sum(-1, keepdims=True)
    bmax = np.maximum(np.abs(c - EF), np.abs(fu - EF).max(-1, keepdims=True))
    L = np.log(1/eps_c); bq = -2*bmax*L/3; cq = -2*sig2*L
    d = (-bq + np.sqrt(bq*bq - 4*n*cq))/(2*n); r = EF - d
    A = alphas[None, None, :]; am1 = A - 1
    with np.errstate(over='ignore', invalid='ignore'):
        Ka = np.where(D < 1000, 2**(am1*np.minimum(D, 1000))*np.log(2**np.minimum(D, 1000) + np.e**2)**3/(6*(2 - A)**3*LN2), np.inf)
    tot = r - (A/am1)*log_err/n - am1*V - am1**2*Ka
    tot = np.where(np.isfinite(tot), tot, -np.inf); j = np.unravel_index(np.argmax(tot), tot.shape); R = float(tot[j])
    R = max(R, 0.0)
    if not full: return R
    return R, dict(gamma=float(gammas[j[0]]), alpha=float(alphas[j[2]]), c_perp=float(c[j[0], j[1], 0]), delta=float(d[j[0], j[1], 0]),
                   V=float(V[j[0], j[1], 0]), Var=float(Var[j[0], j[1], 0]), D=float(D[j[0], j[1], 0]), h=float(h_hon))
