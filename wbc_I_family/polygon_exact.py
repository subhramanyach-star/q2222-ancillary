"""Exact polygon (outer cone C_s) for the I and J families: vertices as rational functions of the parameter,
exact sign analysis (real-root isolation) of all 108-type conditions: mu_{a|x} > 0 at vertices and
1 - g.L > 0 for non-tight local points, on a parameter interval. Also transitions."""
import sympy as sp, pickle, sys
from sym_family import build, vertex_data, local_lines, F
from exact_data import DET, lab
def sign_on(expr, x, lo, hi):
    """returns ('pos'|'neg'|'mixed'|'zero', roots in (lo,hi)) for rational function expr of x on open (lo, hi], exact."""
    num, den = sp.fraction(sp.factor(sp.cancel(expr)))
    if num == 0: return 'zero', []
    roots = []
    for P in (sp.Poly(num, x), sp.Poly(den, x)):
        for r in P.real_roots():
            if r > lo and r < hi: roots.append(r)
    mid = (lo + hi)/2 if not roots else None
    samples = [lo + (hi - lo)*sp.Rational(k, 7) for k in range(1, 8)]
    sg = set()
    pts = sorted(set(roots))
    cuts = [lo] + pts + [hi]
    for a, b in zip(cuts[:-1], cuts[1:]):
        q = (a + b)/2
        v = sp.nsimplify(expr.subs(x, q)); sg.add(sp.sign(sp.N(v, 50)))
    vhi = expr.subs(x, hi); 
    if not pts and sp.sign(sp.N(vhi, 50)) in sg | {0}: pass
    if pts: return ('touch' if len(sg) == 1 else 'mixed'), pts
    if sg == {1}: return 'pos', pts
    if sg == {-1}: return 'neg', pts
    return 'mixed', pts
def analyse(name, lo, hi, verts):
    B = build(name); x = B['x']; LL = local_lines(B); res = {}
    for vn, (l1, l2) in verts.items():
        V = vertex_data(B, l1, l2)
        conds = {}
        for k, mv in enumerate(V['mu']): conds['mu%d' % k] = sign_on(mv, x, lo, hi)
        for L in DET:
            if lab(L) in (l1, l2): continue
            conds['L' + lab(L)] = sign_on(F(1 - (V['g'].T*L)[0]), x, lo, hi)
        bad = {k: v for k, v in conds.items() if v[0] != 'pos'}
        print(name, vn, (l1, l2), 'consistency', V['consistency'], 'non-positive conditions:', bad, flush=True)
        res[vn] = dict(V=V, conds=conds)
    return B, res
if __name__ == '__main__':
    t, mu = sp.symbols('t mu', positive=True)
    VI = {'W1': ('-++-', '----'), 'W2': ('+--+', '----'), 'W3': ('++++', '+--+'), 'W4': ('++++', '-++-')}
    BI, RI = analyse('I', sp.Integer(0), 2 - sp.sqrt(3), VI)
    __import__('polygon_io').dump(dict(B=BI, R=RI), 'out/polygon_I.json')
