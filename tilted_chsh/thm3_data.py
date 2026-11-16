"""Loader for thm3_data.json: returns the data as exact sympy expressions, in the same layout that the scripts use
(D['p'], D['c0'], D['m'], D['gens'], D['v'], D['verts'][name]['g'|'mu'|'d0'|'d1'|'cons'])."""
import json, sympy as sp
s, C, r = sp.symbols('s C r', positive=True)
_loc = {'s': s, 'C': C, 'r': r}
def _e(x): return sp.sympify(x, locals=_loc)
def load(path='thm3_data.json'):
    J = json.load(open(path))
    D = {'p': [_e(x) for x in J['p_theta']], 'c0': _e(J['c0']), 'm': [_e(x) for x in J['m']],
         'gens': [tuple(k) for k, _ in J['gens']], 'v': {tuple(k): [_e(x) for x in v] for k, v in J['gens']},
         'verts': {}}
    for name, V in J['vertices'].items():
        D['verts'][name] = {k: [_e(x) for x in V[k]] for k in ('g', 'mu', 'd0', 'd1')}
        D['verts'][name]['cons'] = {tuple(k): _e(v) for k, v in V['cons']}
    return D
