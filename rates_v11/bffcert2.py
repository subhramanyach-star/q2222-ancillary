"""Dual-certified BFF lower bounds on the GLOBAL (two-party) entropy H(AB|X=x*,Y=y*,E) in a 2 x nB scenario
(Alice 2 binary settings, Bob nB binary settings; nB=3 adds Bob's key setting B2 as in Acin-Massar-Pironio 2012 /
Brown-Fawzi-Fawzi 2021).  Same method as bffcert.py (BFF Lemma 2.3 Gauss-Radau, separate infimum per node, real NPA
relaxation with {1,A,B,AB} x {Z_k, Z_k^*} monomials, SDPA, exact-rational residual certification), generalised to the
four outcomes ab of (A_{x*},B_{y*}).
Coordinates z = (<A_x> (2), <B_y> (nB), <A_x B_y> (2 nB, x-major)).
Constraints: equalities z_k = zeq_k for k in Eq, and optionally one inequality g.z >= omega."""
import numpy as np, scipy.sparse as sps, cvxpy as cp
from fractions import Fraction
from bffcert import gauss_radau, canon, adj, nZ

class BFF2:
    def __init__(self, nB=3, xs=1, ys=2):
        self.nA, self.nB, self.xs, self.ys = 2, nB, xs, ys
        A = [('A', x) for x in range(2)]; B = [('B', y) for y in range(nB)]
        S = [()] + [(a,) for a in A] + [(b,) for b in B] + [(a, b) for a in A for b in B]
        S += [(A[0], A[1]), (A[1], A[0])] + [(b1, b2) for b1 in B for b2 in B if b1 != b2]
        self.K = 4
        for k in range(self.K):
            for z in (('Z', k), ('Zs', k)):
                S += [(z,)] + [(a, z) for a in A] + [(b, z) for b in B] + [(a, b, z) for a in A for b in B]
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
        self.pos = {}
        for (i, j), c in zip(rows, cols):
            self.pos.setdefault(c, []).append((i, j))
            if i != j: self.pos[c].append((j, i))
        self.wlist = [None] * N
        for w, k in words.items(): self.wlist[k] = w
        self.k0 = words[()]
        d = 2 + nB + 2 * nB; self.d = d
        Az = np.zeros((d, N)); a0 = np.zeros(d)      # letters are projectors P = (1+A)/2 (outcome +1 <-> a=0)
        def e(w): v = np.zeros(N); wd = canon(w); v[words[min(wd, adj(wd))]] = 1; return v
        PA = [e((('A', x),)) for x in range(2)]; PB = [e((('B', y),)) for y in range(nB)]
        for x in range(2): Az[x] = 2 * PA[x]; a0[x] = -1
        for y in range(nB): Az[2 + y] = 2 * PB[y]; a0[2 + y] = -1
        for x in range(2):
            for y in range(nB):
                k = 2 + nB + nB * x + y; Az[k] = 4 * e((('A', x), ('B', y))) - 2 * PA[x] - 2 * PB[y]; a0[k] = 1
        self.Az = Az; self.a0 = a0                   # z = Az u + a0
        self.c0, self.c1 = self._obj(0.0), self._obj(1.0)
        self.nz = np.array([nZ(w) for w in self.wlist])
        self.diagw = [self.words[min(canon(adj(w) + w), adj(canon(adj(w) + w)))] for w in self.S]
        self._built = {}
    def _obj(self, t):
        c = np.zeros(self.N)
        def add(w, val):
            wd = canon(w); c[self.words[min(wd, adj(wd))]] += val
        Ax = ('A', self.xs); By = ('B', self.ys)
        for k, (sa, sb) in enumerate([(1, 1), (1, -1), (-1, 1), (-1, -1)]):
            # M_ab = Pa Qb with P0 = P, P1 = 1 - P (projector letters)
            fa = [(1, (Ax,))] if sa == 1 else [(1, ()), (-1, (Ax,))]
            fb = [(1, (By,))] if sb == 1 else [(1, ()), (-1, (By,))]
            terms = [(ca * cb, wa + wb) for ca, wa in fa for cb, wb in fb]
            for cf, w in terms:
                add(w + (('Z', k),), 2 * cf)
                add(w + (('Zs', k), ('Z', k)), cf * (1 - t))
            add((('Z', k), ('Zs', k)), t)
        return c
    def alpha(self, t): return 1.5 * max(1 / t, 1 / (1 - t))
    zscale_pow = 0.5
    def zs(self, t): return self.alpha(t) ** self.zscale_pow
    def beta(self, t): return self.alpha(t) / self.zs(t)
    def cvec(self, t): return (self.c0 + t * (self.c1 - self.c0)) * self.zs(t) ** self.nz
    def build(self, Eq, has_g):
        key = (tuple(Eq), has_g)
        if key in self._built: return self._built[key]
        u = cp.Variable(self.N); G = cp.reshape(self.M @ u, (self.n, self.n), order='C'); G = (G + G.T) / 2
        Wp = cp.Parameter((self.n, self.n), nonneg=True); Wp.value = np.ones((self.n, self.n))
        tt = cp.Parameter(self.N); psd = cp.multiply(Wp, G) >> 0
        dg = sorted(set(k for k, w in zip(self.diagw, self.S) if nZ(w) > 0)); self.dg = dg
        db = cp.Parameter(len(dg), nonneg=True)
        cons = [psd, u[self.k0] == 1, u[dg] <= db]
        zeq = cp.Parameter(len(Eq)) if Eq else None
        if Eq: cons.append(self.Az[list(Eq)] @ u + self.a0[list(Eq)] == zeq)
        gp = om = None
        if has_g:
            gp = cp.Parameter(self.d); om = cp.Parameter(); cons.append(gp @ (self.Az @ u + self.a0) >= om)
        prob = cp.Problem(cp.Minimize(tt @ u), cons)
        self._built[key] = (prob, u, tt, psd, zeq, gp, om, db, Wp); return self._built[key]
    def solve_node(self, t, Eq, zeq, g=None, omega=None, solver='SDPA', opts=None, D=None):
        prob, u, tt, psd, zq, gp, om, db, Wp = self.build(Eq, g is not None)
        tt.value = self.cvec(t)
        Wp.value = np.ones((self.n, self.n)) if D is None else np.outer(1 / D, 1 / D)
        db.value = np.array([self.beta(t) ** self.nz[k] for k in self.dg])
        if Eq: zq.value = np.asarray(zeq, float)
        if g is not None: gp.value = np.asarray(g, float); om.value = float(omega)
        try: prob.solve(solver=solver, **(opts or {}))
        except Exception as ex: return dict(status='error:' + str(ex)[:80])
        if psd.dual_value is None: return dict(status=str(prob.status))
        S = np.array(psd.dual_value) * Wp.value
        diag = None
        if u.value is not None:
            Gv = (self.M @ u.value).reshape(self.n, self.n); diag = np.sqrt(np.maximum(np.diag(Gv), 1e-8))
        return dict(status=prob.status, val=prob.value, S=S, diag=diag)
    def certify(self, t, S, Eq, zeq, g=None, omega=None):
        n = self.n; al = self.beta(t); Eq = list(Eq)
        Drow = np.array([al ** nZ(w) for w in self.S]); S = (S + S.T) / 2
        St = Drow[:, None] * S * Drow[None, :]
        lam, V = np.linalg.eigh(St); neg = lam.min()
        St = (V * np.maximum(lam, 0)) @ V.T; St = (St + St.T) / 2
        lmin = np.linalg.eigvalsh(St).min(); margin = 1e-14 * n * np.abs(St).max()
        psd_pen = n * min(0.0, lmin - margin)
        S = St / Drow[:, None] / Drow[None, :]
        c = self.cvec(t); rho = c - self.M.T @ S.reshape(-1)
        mask = np.ones(self.N, bool); mask[self.k0] = False
        cols = [self.Az[k] for k in Eq] + ([np.asarray(g) @ self.Az] if g is not None else [])
        Cm = np.array(cols).T[mask]
        coef = np.linalg.lstsq(Cm, rho[mask], rcond=None)[0]
        if g is not None and coef[-1] < 0:
            coef = np.zeros(len(cols))
            if Eq: coef[:-1] = np.linalg.lstsq(Cm[:, :-1], rho[mask], rcond=None)[0]
        mu = np.zeros(self.d)                       # total multiplier vector on z
        for i, k in enumerate(Eq): mu[k] += coef[i]
        lam_g = coef[-1] if g is not None else 0.0
        if g is not None: mu = mu + lam_g * np.asarray(g, float)
        Sf = [[Fraction(x) for x in row] for row in S]; muf = [Fraction(x) for x in mu]
        r0 = None; res_pen = Fraction(0)
        nzAz = [[(q, Fraction(self.Az[q, k])) for q in range(self.d) if self.Az[q, k] != 0] for k in range(self.N)]
        for k in range(self.N):
            tot = Fraction(c[k]) - sum((Sf[i][j] for (i, j) in self.pos[k]), Fraction(0)) - sum((muf[q] * a for q, a in nzAz[k]), Fraction(0))
            if k == self.k0: r0 = tot
            else: res_pen += abs(tot) * Fraction(al ** nZ(self.wlist[k]))
        const = float(r0 - res_pen - sum(muf[q] * Fraction(self.a0[q]) for q in range(self.d))) + psd_pen
        # inf_i >= const + mu.z for every quantum strategy (mu = sum_eq + lam_g g, lam_g >= 0)
        mu_eq = mu - lam_g * (np.asarray(g, float) if g is not None else 0)
        LB = const + float(sum(mu_eq[k] * zeq[i] for i, k in enumerate(Eq))) + (lam_g * float(omega) if g is not None else 0.0)
        return dict(LB=LB, k0=const, kappa=mu, lam_g=lam_g, neg_eig=float(neg), psd_pen=psd_pen, res_pen=float(res_pen))

class Rate2:
    TRIALS = [(1.0, 'SDPA', None), (0.5, 'SDPA', None), (0.75, 'SDPA', None),
              (1.0, 'SDPA', dict(gammaStar=0.7, betaBar=0.3, epsilonStar=1e-8, epsilonDash=1e-8))]
    def __init__(self, m=8, nB=3, xs=1, ys=2):
        self.m = m; self.t, self.w = gauss_radau(m); self.B = BFF2(nB, xs, ys)
        self.coef = self.w[:-1] / (self.t[:-1] * np.log(2)); self.cm = self.coef.sum()
    def bound(self, Eq, zeq, g=None, omega=None, good_gap=2e-4, rowscale_nodes=2):
        Hn = Hc = K0 = self.cm; kap = np.zeros(self.B.d); info = []
        for i in range(self.m - 1):
            best = None; trials = [(pw, so, op, None) for (pw, so, op) in self.TRIALS]; k = 0
            while k < len(trials):
                pw, so, op, D = trials[k]; k += 1
                self.B.zscale_pow = pw
                s = self.B.solve_node(self.t[i], Eq, zeq, g, omega, solver=so, opts=op, D=D)
                if 'S' not in s: continue
                try: cert = self.B.certify(self.t[i], s['S'], Eq, zeq, g, omega)
                except Exception: continue
                if not np.isfinite(cert['LB']): continue
                val = s['val'] if (s['val'] is not None and np.isfinite(s['val'])) else np.nan
                if best is None or cert['LB'] > best[1]['LB']: best = (s, cert, val, '%s/%.2f/%s' % (so, pw, 'rs' if D is not None else ''))
                gap = val - best[1]['LB']
                if D is None and i < rowscale_nodes and s.get('diag') is not None and k <= len(self.TRIALS) and gap > good_gap:
                    trials.append((pw, so, op, np.maximum(s['diag'], 1e-3)))
                if np.isfinite(gap) and gap < good_gap: break     # good enough: stop trying
            if best is None: return dict(H_num=np.nan, H_cert=np.nan, status='fail node %d' % i)
            s, cert, val, so = best
            Hn += self.coef[i] * val; Hc += self.coef[i] * cert['LB']; K0 += self.coef[i] * cert['k0']; kap += self.coef[i] * cert['kappa']
            info.append(dict(t=float(self.t[i]), solver=so, val=float(val), LB=float(cert['LB']), gap=float(val - cert['LB'])))
        return dict(H_num=Hn, H_cert=Hc, K0=K0, kappa=kap, nodes=info, status='ok')
