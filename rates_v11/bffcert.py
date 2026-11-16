"""Dual-certified Brown-Fawzi-Fawzi lower bounds on H(A_{x*}|X=x*,E) in the 2222 scenario.

BFF Lemma 2.3 (arXiv:2106.13692): H(A|X=x*,E) >= c_m + sum_{i<m} w_i/(t_i ln2) * inf_i, where
   inf_i = inf sum_a <M_a (Z_a + Z_a^* + (1-t_i) Z_a^* Z_a) + t_i Z_a Z_a^*>,  ||Z_a|| <= alpha_i = 1.5 max(1/t_i, 1/(1-t_i)),
and the infimum is over all quantum strategies (projective A_x, B_y, Z commuting with A, B) compatible with the constraint.
We relax each inf_i (speed-up (3) of BFF: separate infimum per node) to a real NPA moment matrix (monomials: NPA level 2
plus {A,B,AB} x {Z_a, Z_a^*}), solve with Clarabel, and then turn the solver's dual PSD matrix S into a RIGOROUS lower
bound valid for every strategy with ||Z_a|| <= alpha_i (see certify()):
   c.u = <S, Gamma(u)> + mu.(A u) + r.u          (exact identity defining the residual r)
       >= n_rows * min(0, lambda_min(D S D)) + mu.(z - a0) + r_0 - sum_{k>0} |r_k| B_k
with D = diag(row-norm bounds), B_k = bound on |<word_k>| (alpha^{#Z letters}), r computed exactly in rationals.
This yields an AFFINE function of the correlator vector z valid for ALL quantum strategies (an EAT min-tradeoff function).
z = (<A0>,<A1>,<B0>,<B1>,<A0B0>,<A0B1>,<A1B0>,<A1B1>) in +-1 observables; projector P_x = (1+A_x)/2 (outcome +1 <-> a=0)."""
import numpy as np, itertools, scipy.sparse as sps, cvxpy as cp
from fractions import Fraction
from numpy.polynomial import legendre as Lg

def gauss_radau(m):
    """m-point Gauss-Radau on [0,1] with fixed node t_m = 1; returns nodes ascending, weights."""
    c = np.zeros(m + 1); c[m - 1] = 1; c[m] = 1
    r = np.sort(np.real(Lg.legroots(c))); r = r[np.abs(r + 1) > 1e-10]
    Pm1 = Lg.legval(r, np.eye(m)[m - 1])
    w_int = (1 - r) / (m**2 * Pm1**2); w_end = 2 / m**2
    x = np.concatenate([[-1.0], r]); w = np.concatenate([[w_end], w_int])
    t = (1 - x) / 2; w = w / 2; o = np.argsort(t); return t[o], w[o]

def canon(word):
    A = [l for l in word if l[0] == 'A']; B = [l for l in word if l[0] == 'B']; Z = [l for l in word if l[0] in ('Z', 'Zs')]
    def red(seq):
        out = []
        for s in seq:
            if out and out[-1] == s: continue
            out.append(s)
        return out
    return tuple(red(A) + red(B) + Z)
def adj(word):
    A = [l for l in word if l[0] == 'A'][::-1]; B = [l for l in word if l[0] == 'B'][::-1]
    Z = [(('Zs' if l[0] == 'Z' else 'Z'), l[1]) for l in word if l[0] in ('Z', 'Zs')][::-1]
    return canon(tuple(A + B + Z))
def nZ(word): return sum(1 for l in word if l[0] in ('Z', 'Zs'))

class BFFcert:
    def __init__(self, xstar=1, extra=False):
        self.xs = xstar
        A = [('A', 0), ('A', 1)]; B = [('B', 0), ('B', 1)]
        S = [()] + [(a,) for a in A] + [(b,) for b in B] + [(a, b) for a in A for b in B] + [(A[0], A[1]), (A[1], A[0]), (B[0], B[1]), (B[1], B[0])]
        for k in range(2):
            for z in (('Z', k), ('Zs', k)):
                S += [(z,)] + [(a, z) for a in A] + [(b, z) for b in B] + [(a, b, z) for a in A for b in B]
        if extra:   # BFF-style 'AZZ' extras: A_{x*} Z^*Z, A_{x*} Z Z^*, Z^*Z, Z Z^*
            for k in range(2):
                for zz in ((('Zs', k), ('Z', k)), (('Z', k), ('Zs', k))):
                    S += [zz, (('A', xstar),) + zz]
        S = [canon(w) for w in S]; self.S = S; n = len(S); self.n = n
        words = {}; rows = []; cols = []
        for i in range(n):
            for j in range(i, n):
                wd = canon(adj(S[i]) + S[j]); rep = min(wd, adj(wd))
                if rep not in words: words[rep] = len(words)
                rows.append((i, j)); cols.append(words[rep])
        self.words = words; N = len(words); self.N = N
        I = [i * n + j for (i, j) in rows] + [j * n + i for (i, j) in rows if i != j]
        J = cols + [c for (i, j), c in zip(rows, cols) if i != j]
        self.M = sps.csr_matrix((np.ones(len(I)), (I, J)), shape=(n * n, N))
        self.pos = {}   # word -> list of (i,j) positions (full matrix)
        for (i, j), c in zip(rows, cols):
            self.pos.setdefault(c, []).append((i, j))
            if i != j: self.pos[c].append((j, i))
        self.wlist = [None] * N
        for w, k in words.items(): self.wlist[k] = w
        self.k0 = words[()]
        # z(u) = Az u + a0
        Az = np.zeros((8, N)); a0 = np.zeros(8)
        def e(w): v = np.zeros(N); wd = canon(w); v[words[min(wd, adj(wd))]] = 1; return v
        PA = [e((('A', x),)) for x in (0, 1)]; PB = [e((('B', y),)) for y in (0, 1)]
        for x in (0, 1): Az[x] = 2 * PA[x]; a0[x] = -1
        for y in (0, 1): Az[2 + y] = 2 * PB[y]; a0[2 + y] = -1
        for x in (0, 1):
            for y in (0, 1):
                Az[4 + 2 * x + y] = 4 * e((('A', x), ('B', y))) - 2 * PA[x] - 2 * PB[y]; a0[4 + 2 * x + y] = 1
        self.Az = Az; self.a0 = a0
        self.c0, self.c1 = self._obj(0.0), self._obj(1.0)   # unscaled objective c(t) = c0 + t (c1 - c0)
        self.nz = np.array([nZ(w) for w in self.wlist])
        self.diagw = [self.words[min(canon(adj(w) + w), adj(canon(adj(w) + w)))] for w in self.S]
        self._built = {}
    def _obj(self, t):
        N = self.N; c = np.zeros(N); xs = self.xs
        def add(w, val):
            wd = canon(w); c[self.words[min(wd, adj(wd))]] += val
        for a in (0, 1):
            terms = [(1, (('A', xs),))] if a == 0 else [(1, ()), (-1, (('A', xs),))]
            for cf, w in terms:
                add(w + (('Z', a),), 2 * cf)                    # <M Z> + <M Z^*> = 2 Re <M Z>
                add(w + (('Zs', a), ('Z', a)), cf * (1 - t))     # <M Z^* Z>
            add((('Z', a), ('Zs', a)), t)
        return c
    def alpha(self, t): return 1.5 * max(1 / t, 1 / (1 - t))
    zscale_pow = 0.5
    def zs(self, t): return self.alpha(t) ** self.zscale_pow          # Z = zs * Zhat
    def beta(self, t): return self.alpha(t) / self.zs(t)              # ||Zhat|| <= beta
    def cvec(self, t):
        """Objective in the scaled variables Zhat = Z/zs: coefficient of word k times zs^{#Z}."""
        return (self.c0 + t * (self.c1 - self.c0)) * self.zs(t) ** self.nz
    def build(self, mode):
        """mode 'score': g.z >= omega ; mode 'full': z == zfix."""
        if mode in self._built: return self._built[mode]
        u = cp.Variable(self.N); G = cp.reshape(self.M @ u, (self.n, self.n), order='C'); G = (G + G.T) / 2
        Wp = cp.Parameter((self.n, self.n), nonneg=True); Wp.value = np.ones((self.n, self.n))
        tt = cp.Parameter(self.N); psd = cp.multiply(Wp, G) >> 0     # row-scaled PSD constraint (same feasible set)
        dg = sorted(set(k for k, w in zip(self.diagw, self.S) if nZ(w) > 0))
        self.dg = dg; db = cp.Parameter(len(dg), nonneg=True)
        cons = [psd, u[self.k0] == 1, u[dg] <= db]        # valid: <w^* w> <= ||w||^2 <= beta^{2 #Z}
        if mode == 'score':
            gp = cp.Parameter(8); om = cp.Parameter(); lin = gp @ (self.Az @ u + self.a0) >= om; cons.append(lin); par = (gp, om)
        else:
            zf = cp.Parameter(8); lin = self.Az @ u + self.a0 == zf; cons.append(lin); par = (zf,)
        prob = cp.Problem(cp.Minimize(tt @ u), cons)
        self._built[mode] = (prob, u, tt, psd, lin, par, db, Wp); return self._built[mode]
    def solve_node(self, t, mode, g=None, omega=None, z=None, solver='SDPA', opts=None, D=None):
        prob, u, tt, psd, lin, par, db, Wp = self.build(mode); tt.value = self.cvec(t)
        Wp.value = np.ones((self.n, self.n)) if D is None else np.outer(1 / D, 1 / D)
        db.value = np.array([self.beta(t) ** (2 * 1.0) if self.nz[k] == 2 else self.beta(t) ** self.nz[k] for k in self.dg])
        if mode == 'score': par[0].value = np.asarray(g, float); par[1].value = float(omega)
        else: par[0].value = np.asarray(z, float)
        try:
            prob.solve(solver=solver, **(opts or {}))
        except Exception as ex:
            return dict(status='error:' + str(ex)[:80])
        if psd.dual_value is None: return dict(status=str(prob.status))
        S = np.array(psd.dual_value) * Wp.value                     # <St, W.G> = <St.W, G>
        diag = None
        if u.value is not None:
            Gv = (self.M @ u.value).reshape(self.n, self.n); diag = np.sqrt(np.maximum(np.diag(Gv), 1e-8))
        return dict(status=prob.status, val=prob.value, S=S, diag=diag)
    def certify(self, t, S, mode, g=None, omega=None, z=None):
        """Rigorous lower bound on inf_i (valid for all strategies with ||Z||<=alpha(t)) from a candidate dual S.
        Returns (LB at the given constraint, affine coefficients (k0, kappa) with inf_i >= k0 + kappa.z for all quantum z
        [score mode: kappa = lam*g, valid for all z with g.z >= omega ... we return the affine form in z anyway, which is
        valid for all quantum z because lam >= 0 only enters through lam*(g.z - omega) >= 0 is NOT used: see below])."""
        n = self.n; al = self.beta(t)                            # scaled variables: ||Zhat|| <= beta
        Drow = np.array([al ** nZ(w) for w in self.S])
        S = (S + S.T) / 2
        # PSD repair in the scaled frame: St = D S D ; clip negative eigenvalues
        St = Drow[:, None] * S * Drow[None, :]
        lam, V = np.linalg.eigh(St); neg = lam.min()
        St = (V * np.maximum(lam, 0)) @ V.T; St = (St + St.T) / 2
        lmin = np.linalg.eigvalsh(St).min()
        margin = 1e-14 * n * np.abs(St).max()                  # floating-point eigenvalue error allowance
        psd_pen = n * min(0.0, lmin - margin)                  # <St, D^-1 G D^-1> >= lmin * tr(...) , tr <= n
        S = St / Drow[:, None] / Drow[None, :]
        c = self.cvec(t)
        MS = self.M.T @ S.reshape(-1)                           # <S, F_k>
        rho = c - MS
        # least squares for multipliers mu (8-vector) on words k != 0 : rho ~ Az^T mu
        mask = np.ones(self.N, bool); mask[self.k0] = False
        if mode == 'score':
            gh = np.asarray(g) @ self.Az; lam_ = max(0.0, float(np.dot(gh[mask], rho[mask]) / np.dot(gh[mask], gh[mask])))
            mu = lam_ * np.asarray(g, float)
        else:
            mu = np.linalg.lstsq(self.Az[:, mask].T, rho[mask], rcond=None)[0]
        # exact residual in rationals: r = c - <S,F_k> - (Az^T mu)_k
        Sf = [[Fraction(x) for x in row] for row in S]
        muf = [Fraction(x) for x in mu]
        r0 = None; res_pen = Fraction(0)
        for k in range(self.N):
            tot = Fraction(c[k]) - sum((Sf[i][j] for (i, j) in self.pos[k]), Fraction(0)) - sum((muf[q] * Fraction(self.Az[q, k]) for q in range(8) if self.Az[q, k] != 0), Fraction(0))
            if k == self.k0: r0 = tot
            else: res_pen += abs(tot) * Fraction(al ** nZ(self.wlist[k]))
        # c.u = <S,G> + mu.(Az u) + r.u  ;  Az u = z - a0
        const = float(r0 - res_pen - sum(muf[q] * Fraction(self.a0[q]) for q in range(8))) + psd_pen
        # inf_i >= const + mu.z  for every quantum strategy (both modes; in score mode mu = lam g, lam>=0, so on g.z>=omega
        # it is >= const + lam*omega)
        if mode == 'score': LB = const + lam_ * float(omega)
        else: LB = const + float(np.dot(mu, z))
        return dict(LB=LB, k0=const, kappa=mu, neg_eig=float(neg), psd_pen=psd_pen, res_pen=float(res_pen))

class Rate:
    def __init__(self, m=8, xstar=1, extra=False):
        self.m = m; self.t, self.w = gauss_radau(m); self.B = BFFcert(xstar, extra)
        self.coef = self.w[:-1] / (self.t[:-1] * np.log(2)); self.cm = self.coef.sum()
    HARD_NODES = 2; ROWSCALE_ITERS = 2
    TRIALS = [(1.0, 'SDPA', None), (0.5, 'SDPA', None), (0.75, 'SDPA', None),
              (1.0, 'SDPA', dict(gammaStar=0.7, betaBar=0.3, epsilonStar=1e-8, epsilonDash=1e-8)),
              (0.5, 'SDPA', dict(gammaStar=0.7, betaBar=0.3, epsilonStar=1e-8, epsilonDash=1e-8))]
    def bound(self, mode, g=None, omega=None, z=None, accept_gap=1e-7):
        """Returns dict: H_num (solver primal value, NOT rigorous), H_cert (rigorous lower bound), and the affine
        min-tradeoff form (K0, kappa) with H(A|X=x*,E) >= K0 + kappa.z for ALL quantum strategies (score mode: kappa || g)."""
        Hn = self.cm; Hc = self.cm; K0 = self.cm; kap = np.zeros(8); info = []
        for i in range(self.m - 1):
            best = None
            trials = [(pw, so, op, None) for (pw, so, op) in self.TRIALS]
            nrs = self.ROWSCALE_ITERS if i < self.HARD_NODES else 0
            k = 0
            while k < len(trials):
                pw, solver, opts, D = trials[k]; k += 1
                self.B.zscale_pow = pw
                s = self.B.solve_node(self.t[i], mode, g, omega, z, solver=solver, opts=opts, D=None if D is None else D[0])
                if 'S' not in s: continue
                if nrs and s.get('diag') is not None and (D is None or D[1] < nrs):
                    trials.append((pw, solver, opts, (np.maximum(s['diag'], 1e-3), 1 if D is None else D[1] + 1)))
                try: cert = self.B.certify(self.t[i], s['S'], mode, g, omega, z)
                except Exception: continue
                if not np.isfinite(cert['LB']): continue
                val = s['val'] if (s['val'] is not None and np.isfinite(s['val'])) else np.nan
                if best is None or cert['LB'] > best[1]['LB']: best = (s, cert, val, '%s/%.2f/rs%d' % (solver, pw, 0 if D is None else D[1]))
            if best is None: return dict(H_num=np.nan, H_cert=np.nan, status='fail node %d' % i)
            s, cert, val, solver = best
            Hn += self.coef[i] * val; Hc += self.coef[i] * cert['LB']
            K0 += self.coef[i] * cert['k0']; kap += self.coef[i] * cert['kappa']
            info.append(dict(t=self.t[i], status=s['status'], solver=solver, val=val, LB=cert['LB'], gap=val - cert['LB'], neg=cert['neg_eig'], pp=cert['psd_pen'], rp=cert['res_pen']))
        return dict(H_num=Hn, H_cert=Hc, K0=K0, kappa=kap, nodes=info, status='ok')
