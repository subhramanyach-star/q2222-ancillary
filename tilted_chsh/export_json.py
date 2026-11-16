"""One-time export of the exact symbolic data (m2_symbolic.pkl, produced by symbolic_m2.py) to the documented,
pickle-free text format thm3_data.json.  All expressions are exact sympy strings in the symbols s, C, r with
rational coefficients and the relations C^2 = 1 - s^2, r^2 = 1 + s^2 (s = sin 2theta, C = cos 2theta, r = sqrt(1+s^2)).
Run:  python3 export_json.py   (requires m2_symbolic.pkl; the JSON is shipped, so this is only needed to regenerate)."""
import pickle, json, sympy as sp
D = pickle.load(open('m2_symbolic.pkl', 'rb'))
S = lambda e: sp.sstr(sp.sympify(e))
out = {
  "format": "thm3_data v1",
  "symbols": {"s": "sin(2*theta)", "C": "cos(2*theta)", "r": "sqrt(1+s**2)"},
  "relations": ["C**2 = 1 - s**2", "r**2 = 1 + s**2"],
  "coordinates": "z = (alpha0, alpha1, beta0, beta1, E00, E01, E10, E11); functional g = (a0, a1, b0, b1, e00, e01, e10, e11)",
  "ax_order": "(a|x) = (+|0), (+|1), (-|0), (-|1)",
  "p_theta": [S(x) for x in D['p']],
  "c0": S(D['c0']),
  "m": [S(x) for x in D['m']],
  "gens": [[list(k), [S(x) for x in D['v'][k]]] for k in D['gens']],
  "vertices": {}}
for name, V in D['verts'].items():
    out["vertices"][name] = {
        "g": [S(x) for x in V['g']], "mu": [S(x) for x in V['mu']],
        "d0": [S(x) for x in V['d0']], "d1": [S(x) for x in V['d1']],
        "cons": [[list(k), S(v)] for k, v in V['cons'].items()]}
json.dump(out, open('thm3_data.json', 'w'), indent=1)
print('wrote thm3_data.json')
