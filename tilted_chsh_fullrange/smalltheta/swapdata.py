"""Exact data for the PARTY-SWAPPED problem (slices by Alice's overlap), in K = Q(s)(C,r) (fullrange/kfield.py).
swap: z=(a0,a1,b0,b1,E00,E01,E10,E11) -> (b0,b1,a0,a1,E00,E10,E01,E11); it maps Q onto Q (relabel parties), so
h_Q(g) = h_Q(swap g) and swap(N_Q(p)) = N_Q(swap p).  For V in {V1, V3}: g' = swap(g_V), target p' = swap(p_theta) whose
Bob' (= Alice) overlap is c0' = 0.  We recompute, with the vertex_data.py recipe (c0 = 0, Gram = identity), the slice-cone
generators at p', the multipliers mu' (g' = sum mu' gamma'), the forced dual delta'(0) = d0' and delta''(0) = d1', and
verify exactly: g' = swap(g_V); and the c0-identities of fullrange/exact_c0_all.py at c0' = 0."""
import sys, json, sympy as sp
sys.path.insert(0, '../fullrange')
from kfield import K, from_sympy
S_, C_, R_ = sp.symbols('s C r', positive=True); loc = {'s': S_, 'C': C_, 'r': R_}
J = json.load(open('../fullrange/thm3_data.json'))
E = lambda x: from_sympy(sp.sympify(x, locals=loc))
def swap(z): return [z[2], z[3], z[0], z[1], z[4], z[6], z[5], z[7]]
p = swap([E(x) for x in J['p_theta']])
KEYS = [(1, 0), (1, 1), (-1, 0), (-1, 1)]
FAC = {'V1': ((1, 1, 1, 1), (1, -1, -1, 1)), 'V3': ((1, -1, 1, 1), (-1, 1, 1, -1))}
def Lvec(A0, A1, B0, B1): return [K(x) for x in (A0, A1, B0, B1, A0 * B0, A0 * B1, A1 * B0, A1 * B1)]
def dot(u, v):
    s = K(0)
    for a, b in zip(u, v): s = s + a * b
    return s
def solve(A, b):
    n = len(A); M = [row[:] + [bb] for row, bb in zip(A, b)]
    for i in range(n):
        piv = next(r for r in range(i, n) if not M[r][i].iszero()); M[i], M[piv] = M[piv], M[i]
        inv = K(1) / M[i][i]; M[i] = [x * inv for x in M[i]]
        for r in range(n):
            if r != i and not M[r][i].iszero():
                f = M[r][i]; M[r] = [x - f * y for x, y in zip(M[r], M[i])]
    return [M[i][n] for i in range(n)]
gens = {}; vs = {}
for (a, x) in KEYS:
    t = K(1) + K(a) * p[x]; w = [p[2] + K(a) * p[4 + 2 * x], p[3] + K(a) * p[5 + 2 * x]]
    assert (w[0] * w[0] + w[1] * w[1] - t * t).iszero()
    v = [w[0] / t, w[1] / t]
    g = [K(0)] * 8; g[x] = K(-a); g[2], g[3] = v; g[4 + 2 * x] = K(a) * v[0]; g[5 + 2 * x] = K(a) * v[1]
    gens[(a, x)] = g; vs[(a, x)] = v
zd = [K(0)] * 8; zd[3] = p[2]; zd[5] = p[4]; zd[7] = p[6]       # c0 = 0
m = [dot(gens[k], zd) for k in KEYS]
out = {'c0': '0', 'p': [str(x.to_sympy(S_, C_, R_)) for x in p], 'vertices': {}}
for name in sys.argv[1:]:
    l1, l2 = FAC[name]; sw = lambda L: (L[2], L[3], L[0], L[1])
    L1, L2 = Lvec(*sw(l1)), Lvec(*sw(l2))
    A = [[K(1)] * 4, m, [dot(gens[k], L1) for k in KEYS], [dot(gens[k], L2) for k in KEYS]]
    mu = solve(A, [K(1), K(0), K(1), K(1)])
    g = [sum((mu[i] * gens[KEYS[i]][j] for i in range(4)), K(0)) for j in range(8)]
    gV = swap([E(x) for x in J['vertices'][name]['g']])
    assert all((a - b).iszero() for a, b in zip(g, gV)), 'swapped vertex mismatch'
    b = [g[2] / 2, g[3] / 2]
    eta = {k: [mu[i] * vs[k][0], mu[i] * vs[k][1]] for i, k in enumerate(KEYS)}
    d0 = [eta[(1, 0)][j] + eta[(-1, 0)][j] - b[j] for j in range(2)]
    rows = []; rhs = []
    for x in (0, 1):
        sx = K(1) if x == 0 else K(-1)
        u = {s: [2 * eta[(s, x)][0], 2 * eta[(s, x)][1]] for s in (1, -1)}
        nr = {s: 2 * mu[KEYS.index((s, x))] for s in (1, -1)}
        coef = [K(0), K(0)]; const = K(0)
        for s in (1, -1):
            G0 = [u[s][0], u[s][1]]   # c0 = 0
            coef = [coef[j] + K(s) / nr[s] * sx * G0[j] for j in range(2)]
            const = const + K(s) / nr[s] * u[s][0] * u[s][1]
        rows.append(coef); rhs.append(-const)
    d1 = solve(rows, rhs)
    # exact c0-identities (as fullrange/exact_c0_all.py), c0 = 0
    M = {k: mu[i] for i, k in enumerate(KEYS)}; ok = True
    for s0 in (1, -1):
        for s1 in (1, -1):
            u = [g[2 + i] / 2 + d0[i] + g[4 + i] * K(s0) for i in range(2)]
            w = [g[2 + i] / 2 - d0[i] + g[6 + i] * K(s1) for i in range(2)]
            Aq = u[0] * u[0] + u[1] * u[1]; Bq = w[0] * w[0] + w[1] * w[1]
            dA = K(2) * (u[0] * d1[0] + u[1] * d1[1]) + K(2) * u[0] * u[1]
            dB = K(-2) * (w[0] * d1[0] + w[1] * d1[1]) + K(2) * w[0] * w[1]
            T = K(1) - K(s0) * g[0] - K(s1) * g[1]; m0, m1 = M[(s0, 0)], M[(s1, 1)]
            ids = [Aq - K(4) * m0 * m0, Bq - K(4) * m1 * m1, T - K(2) * m0 - K(2) * m1, dA * m1 + dB * m0]
            ok &= all(z.iszero() for z in ids)
    print(name, 'c0-identities exact:', ok)
    f = lambda z: str(sp.simplify(z.to_sympy(S_, C_, R_)))
    out['vertices'][name] = {'g': [f(z) for z in g], 'mu': [f(z) for z in mu], 'd0': [f(z) for z in d0], 'd1': [f(z) for z in d1]}
    print(' mu =', out['vertices'][name]['mu']); print(' d0 =', out['vertices'][name]['d0']); print(' d1 =', out['vertices'][name]['d1'])
json.dump(out, open('swap_data.json', 'w'), indent=1)
