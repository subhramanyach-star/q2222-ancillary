"""Finite-size rates via the entropy accumulation theorem in the form of Liu et al., Nat. Phys. 17, 448 (2021)
[arXiv:1912.11159], App. A, Theorem 3 (EAT with improved second-order handling, built on Dupuis-Fawzi-Renner CMP 2020 and
Dupuis-Fawzi IEEE TIT 2019) combined with their spot-checking Lemma 3, generalised from the CHSH game to an arbitrary
affine min-tradeoff function of the correlators.  We use the conservative, statistics-agnostic simplifications:
Delta(f,p) >= 0 is dropped and V(f,p) is replaced by its supremum over all no-signalling test distributions.

Per-round entropy (generation rounds, input x*=1 for Alice): H(A|X=x*,E) >= K0 + kappa.z for ALL quantum z (BFF dual
certificate, or tangent of the analytic CHSH curve).  Rounds are tests with prob. gamma (uniform x,y; score register
U = (x,y,a,b)), so per-round g(q) = (1-gamma)(K0 + kappa.z(q)), with test-round value g(delta_u) = (1-gamma)(K0 + s(u)),
s(u) = 4 kappa_xy a b + 2 kappa_Ax a + 2 kappa_By b  (unbiased estimator of kappa.z).
Theorem 3:  H_min^{eps_h}(AB|XYTE)_{rho|Omega} > n r - a/(a-1) log2(1/(p_Omega (1 - sqrt(1-eps_h^2))))
                                             - n (a-1) V - n (a-1)^2 K_a,
with V = ln2/2 (log2(1+2 d_A) + sqrt(2 + Var))^2, d_A = 4 (A_i = (a_i,b_i)), Var = sup_q (1/gamma) sum_u q(u)(c - g(delta_u))^2,
K_a = 2^{(a-1)(log2 d_A + Max f - Min f)} ln^3(2^{log2 d_A + Max f - Min f} + e^2) / (6 (2-a)^3 ln 2),
Max f = max(c, c + (Max g - c)/gamma), Min f >= min_NS g;  p_Omega -> eps_EAT (Liu et al. A4).
Omega = {f(freq) >= r}, r = g(q_honest) - delta, delta from Bernstein so that the honest abort probability <= eps_c."""
import numpy as np, itertools
LOG2_9 = np.log2(9.0)
def ns_vertices():
    """Correlator vectors z of the 16 local deterministic and 8 PR boxes (vertices of the 2222 NS polytope)."""
    V = []
    for A0, A1, B0, B1 in itertools.product([1, -1], repeat=4):
        V.append(np.array([A0, A1, B0, B1, A0*B0, A0*B1, A1*B0, A1*B1], float))
    for s in itertools.product([1, -1], repeat=3):
        # PR boxes: zero marginals, E_xy = s_xy with s00 s01 s10 s11 = -1
        e = list(s) + [-s[0]*s[1]*s[2]]
        V.append(np.array([0, 0, 0, 0] + e, float))
    return V
NSV = ns_vertices()
def test_outcomes():
    """List of (x,y,a,b, prob weight factor 1/4) and the per-outcome estimator vector e_u with s(u) = kappa.e_u."""
    out = []
    for x in (0, 1):
        for y in (0, 1):
            for a in (1, -1):
                for b in (1, -1):
                    e = np.zeros(8); e[x] = 2*a; e[2 + y] = 2*b; e[4 + 2*x + y] = 4*a*b
                    out.append((x, y, a, b, e))
    return out
TO = test_outcomes()
def q_of_z(z):
    """Test-round distribution q(u) = p(ab|xy)/4 for correlators z."""
    q = []
    for (x, y, a, b, e) in TO:
        p = (1 + a*z[x] + b*z[2 + y] + a*b*z[4 + 2*x + y]) / 4
        q.append(p / 4)
    return np.array(q)
def hbin(p):
    p = np.clip(p, 1e-300, 1 - 1e-16); return -p*np.log2(p) - (1 - p)*np.log2(1 - p)
def eat_rate(K0, kappa, z_hon, n, eps_h=1.545e-12, eps_eat=3.09e-12, eps_c=1e-6,
             gammas=None, alphas=None, cfrac=None):
    """Best (over gamma, alpha, c_perp) EAT lower bound on H_min^{eps_h}/n. Returns dict."""
    kappa = np.asarray(kappa, float); z_hon = np.asarray(z_hon, float)
    if gammas is None: gammas = np.logspace(-5, -0.3, 48)
    if alphas is None: alphas = 1 + np.logspace(-10, -0.5, 96)
    if cfrac is None: cfrac = np.linspace(-2.0, 1.5, 36)
    su = np.array([kappa @ e for (*_, e) in TO])                # s(u)
    q_hon = q_of_z(z_hon)
    h_hon = K0 + kappa @ z_hon                                   # asymptotic per-gen-round rate
    best = dict(rate=-np.inf)
    if h_hon <= 0: return dict(rate=0.0, h=h_hon)
    qNS = [q_of_z(v) for v in NSV]
    log_err = np.log2(1/eps_eat) + np.log2(1/(eps_h**2 / (1 + np.sqrt(1 - eps_h**2))))   # log2 1/(p_Omega (1 - sqrt(1-eps^2)))
    for gam in gammas:
        gu = (1 - gam)*(K0 + su)                                  # g(delta_u)
        g_hon = (1 - gam)*h_hon
        Maxg = gu.max(); Ming = min(float(q @ gu) for q in qNS)
        for cf in cfrac:
            c = Ming + cf*(Maxg - Ming)                           # c_perp
            Var = max(float(q @ (c - gu)**2) for q in qNS)/gam
            V = np.log(2)/2*(LOG2_9 + np.sqrt(2 + Var))**2
            Maxf = max(c, c + (Maxg - c)/gam); D = 2 + Maxf - Ming
            # completeness: F = f(delta_U), honest i.i.d.
            fu = gu/gam + (1 - 1/gam)*c
            vals = np.concatenate([[c], fu]); probs = np.concatenate([[1 - gam], gam*q_hon*4/4])
            probs = probs/probs.sum(); EF = probs @ vals
            sig2 = probs @ (vals - EF)**2; bmax = np.abs(vals - EF).max()
            L = np.log(1/eps_c)
            # solve n d^2 = L (2 sig2 + 2 b d /3)
            aq = n; bq = -2*bmax*L/3; cq = -2*sig2*L
            d = (-bq + np.sqrt(bq*bq - 4*aq*cq))/(2*aq)
            r = EF - d
            # optimise alpha
            am1 = alphas - 1
            with np.errstate(over='ignore', invalid='ignore'):
                Ka = 2**(am1*D)*np.log(2**D + np.e**2)**3/(6*(2 - alphas)**3*np.log(2)) if D < 1000 else np.full_like(am1, np.inf)
            tot = r - (alphas/am1)*log_err/n - am1*V - am1**2*Ka
            j = int(np.nanargmax(tot)); R = tot[j]
            if R > best['rate']:
                best = dict(rate=float(R), gamma=float(gam), alpha=float(alphas[j]), c_perp=float(c), delta=float(d),
                            r=float(r), V=float(V), Var=float(Var), D=float(D), h=float(h_hon),
                            net=float(R - hbin(gam) - 2*gam))
    if best['rate'] < 0: best['rate'] = 0.0
    return best
def chsh_analytic(S):
    S = np.minimum(S, 2*np.sqrt(2))
    x = np.sqrt(np.maximum(S*S/4 - 1, 0)); return 1 - hbin(0.5 + 0.5*x)
def chsh_tangent(S0):
    """Affine lower bound H >= K0 + k*S tangent to the convex analytic CHSH curve at S0 (valid for all quantum S)."""
    h = 1e-6; F = chsh_analytic(S0); dF = (chsh_analytic(S0 + h) - chsh_analytic(S0 - h))/(2*h)
    g = np.array([0, 0, 0, 0, 1, 1, 1, -1.0]); return F - dF*S0, dF*g
if __name__ == '__main__':
    s2 = np.sqrt(2); zT = np.array([0, 0, 0, 0, 1, 1, 1, -1])/s2
    for v in (0.99, 0.98):
        z = v*zT; S = v*2*s2; K0, kap = chsh_tangent(S)
        for n in (1e8, 1e10, 1e12):
            b = eat_rate(K0, kap, z, n); print('CHSH Tsirelson v=%.2f n=%.0e asym %.4f rate %.4f gamma %.1e alpha-1 %.1e' % (v, n, b['h'], b['rate'], b.get('gamma', 0), b.get('alpha', 1) - 1))
