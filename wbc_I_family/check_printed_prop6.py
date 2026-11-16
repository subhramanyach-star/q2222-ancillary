"""Exact (sympy) comparison of the formulas printed in Proposition 6 / App. F with out/polygon_I.json, and the Gram
determinant used for exposure.  Checks: vertex d-coordinates W1, W2, -W1, -W2 and their tight deterministic points;
Lambda; g_corr = I_delta / I_delta^Q; c0 = -sin delta; Gram determinant of (L_i - p_delta), i = 1..4, equals
4096[2(1 - cos d) - sin^4 d] = 16384 sin^2(d/2)[1 - sin^2 d cos^2(d/2)].   Run: python3 check_printed_prop6.py"""
import sympy as sp, polygon_io
D = polygon_io.load('out/polygon_I.json'); B = D['B']; R = D['R']; t = B['x']
dl = sp.symbols('delta', positive=True); c, s, c2 = sp.cos(dl), sp.sin(dl), sp.cos(2*dl)
S = 2*t/(1 + t**2); C = (1 - t**2)/(1 + t**2)
def eq(a, b): return sp.simplify(a - b.subs({sp.sin(2*dl): 2*S*C, sp.cos(2*dl): C**2 - S**2}).subs({sp.sin(dl): S, sp.cos(dl): C})) == 0
W1 = (-(1 - c)*c2/(2*c**3), s/(2*c**3)); W2 = (-s*c2/(2*c**3), (c - c2)/(2*c**3))
want = {'W1': (W1, ('-++-', '----')), 'W2': (W2, ('+--+', '----')), 'W3': ((-W1[0], -W1[1]), ('++++', '+--+')), 'W4': ((-W2[0], -W2[1]), ('++++', '-++-'))}
ok = True
for k, (w, tight) in want.items():
    V = R[k]['V']; r = all(eq(V['d'][i], w[i]) for i in (0, 1)) and tuple(V['tight']) == tight; ok &= r
    print(k, '(printed as', {'W1': 'W1', 'W2': 'W2', 'W3': '-W1', 'W4': '-W2'}[k] + ')', 'd and tight points match:', r)
r = eq(B['Lam'][0], c2/(2*c**2)) and eq(B['Lam'][1], 1/(2*c**2)); ok &= r; print('Lambda matches:', r)
IQ = 2*c**3/(c2*s); I = [0, 0, 0, 0, 1, 1/s, 1/s, -1/c2]
r = all(eq(B['gcorr'][i], I[i]/IQ) for i in range(8)); ok &= r; print('g_corr = I/I^Q:', r)
r = eq(B['c0'], -s); ok &= r; print('c0 = -sin delta:', r)
p = sp.Matrix([0, 0, 0, 0, 0, c, c, -2*s*c])
def L(a0, a1, b0, b1): return sp.Matrix([a0, a1, b0, b1, a0*b0, a0*b1, a1*b0, a1*b1])
M = sp.Matrix.hstack(*[L(*q) - p for q in ((-1, 1, 1, -1), (-1, -1, -1, -1), (1, -1, -1, 1), (1, 1, 1, 1))])
G = sp.simplify((M.T*M).det())
u = sp.symbols('u', positive=True)  # u = delta/2
r1 = sp.simplify(G - 4096*(2*(1 - c) - s**4)) == 0
H = (G - 16384*sp.sin(dl/2)**2*(1 - s**2*sp.cos(dl/2)**2)).subs(dl, 2*u)
r2 = sp.simplify(sp.expand(sp.expand_trig(H)).subs(sp.cos(u)**2, 1 - sp.sin(u)**2)) == 0
ok &= r1 and r2; print('Gram determinant = 4096[2(1-cos d)-sin^4 d]:', r1, '; = 16384 sin^2(d/2)[1-sin^2 d cos^2(d/2)]:', r2)
print('ALL PRINTED PROP. 6 DATA MATCH:', ok)
