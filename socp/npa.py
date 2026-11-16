# copy of the authors' npa.py with tightened Clarabel tolerances (benchmark use)
"""NPA hierarchy (Navascues-Pironio-Acin 2008) for the 2222 scenario, real moment matrices.
Upper bounds on h_Q(g) = max over Q of g.z (coordinates as in core.py)."""
import itertools
import numpy as np
import cvxpy as cp

def _reduce(seq):
    out = []
    for s in seq:
        if out and out[-1] == s:
            out.pop()
        else:
            out.append(s)
    return tuple(out)

def canon(word):
    A = _reduce([w for w in word if w[0] == 'A'])
    B = _reduce([w for w in word if w[0] == 'B'])
    w1 = A + B
    w2 = tuple(reversed(A)) + tuple(reversed(B))   # adjoint (A,B commute)
    return min(w1, w2)

def operator_list(level):
    letters = [('A', 0), ('A', 1), ('B', 0), ('B', 1)]
    if level == '1+AB':
        S = [()] + [(l,) for l in letters] + [(('A', x), ('B', y)) for x in (0, 1) for y in (0, 1)]
        return S
    S = set([()])
    for L in range(1, level + 1):
        for w in itertools.product(letters, repeat=L):
            S.add(canon_nonsym(w))
    S = sorted(S, key=lambda w: (len(w), w))
    return S

def canon_nonsym(word):
    A = _reduce([w for w in word if w[0] == 'A'])
    B = _reduce([w for w in word if w[0] == 'B'])
    return A + B

def npa_value(g, level='1+AB', solver='CLARABEL'):
    g = np.asarray(g, float)
    S = operator_list(level)
    n = len(S)
    keys = {}
    for i in range(n):
        for j in range(i, n):
            w = canon(tuple(reversed(S[i])) + S[j])
            keys.setdefault(w, []).append((i, j))
    G = cp.Variable((n, n), symmetric=True)
    cons = [G >> 0]
    for w, lst in keys.items():
        i0, j0 = lst[0]
        for (i, j) in lst[1:]:
            cons.append(G[i, j] == G[i0, j0])
    cons.append(G[0, 0] == 1)
    idx = {S[k]: k for k in range(n)}
    def mom(word):
        w = canon(word)
        i, j = keys[w][0]
        return G[i, j]
    z = [mom((('A', 0),)), mom((('A', 1),)), mom((('B', 0),)), mom((('B', 1),)),
         mom((('A', 0), ('B', 0))), mom((('A', 0), ('B', 1))), mom((('A', 1), ('B', 0))), mom((('A', 1), ('B', 1)))]
    obj = sum(g[k] * z[k] for k in range(8))
    prob = cp.Problem(cp.Maximize(obj), cons)
    if solver == 'CLARABEL':
        prob.solve(solver=solver, tol_gap_abs=1e-10, tol_gap_rel=1e-10, tol_feas=1e-10, max_iter=1000)
    else:
        prob.solve(solver=solver)
    return prob.value
