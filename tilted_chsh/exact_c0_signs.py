"""Interval side conditions for exact_c0.py on Theta = [pi/8, 3pi/16]: R > 0, A > 0, B > 0, R^2 - A - B > 0,
A'B' < 0 (or both enclose... must be strict).  Together with I1 = I2 = 0 (exact) this proves F(c0)=F_c(c0)=0."""
import pickle, flint
from flint import arb, fmpq
from arbeval import compile_expr, S_, C_, R_
flint.ctx.prec = 110
out = {}
import json, sympy as sp
for f in ('exact_c0_V1.json', 'exact_c0_V3.json'):
    for k, v in json.load(open(f)).items():
        out[tuple(k.split('|'))] = {kk: sp.sympify(vv, locals={'s': S_, 'C': C_, 'r': R_}) for kk, vv in v.items()}
PI = arb.pi()
res = {}
for key, v in out.items():
    fs = {k: compile_expr(v[k], [S_, C_, R_]) for k in ('R', 'A', 'B', 'dA', 'dB')}
    stack = [(fmpq(1, 8) + fmpq(j, 16 * 64), fmpq(1, 8) + fmpq(j + 1, 16 * 64)) for j in range(64)]
    nb = 0; mins = [1e9] * 5
    while stack:
        a, b = stack.pop(); x = arb((a + b) / 2, (b - a) / 2); th = x * PI
        sv = (2 * th).sin(); vv = [sv, (2 * th).cos(), (1 + sv * sv).sqrt()]
        E = {k: f(vv) for k, f in fs.items()}
        conds = [E['R'], E['A'], E['B'], E['R'] ** 2 - E['A'] - E['B'], -E['dA'] * E['dB']]
        if all(cnd > 0 for cnd in conds):
            nb += 1; mins = [min(m, float(cnd.lower())) for m, cnd in zip(mins, conds)]; continue
        if b - a < fmpq(1, 10**7): raise SystemExit('FAIL %s at %s' % (key, float(a)))
        m = (a + b) / 2; stack += [(a, m), (m, b)]
    res[key] = (nb, mins); print(key, 'OK, boxes', nb, 'min lower bounds [R, A, B, R^2-A-B, -A\'B\'] =', ['%.3g' % m for m in mins], flush=True)
