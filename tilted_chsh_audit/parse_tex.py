"""Parse the printed Appendix C data (V1, V3 formulas) and Table III directly from prl.tex."""
import re, sympy as sp
from sympy.parsing.sympy_parser import parse_expr, standard_transformations, implicit_multiplication_application
TR = standard_transformations + (implicit_multiplication_application,)
q, C, r, D1, D3 = sp.symbols('q C r D1 D3')
import os
_HERE = os.path.dirname(os.path.abspath(__file__))
# The paper source is expected two levels up (arXiv layout: <source>/anc/tilted_chsh_audit/); override with $PRL_TEX.
TEX = os.environ.get('PRL_TEX') or next((p for p in (os.path.join(_HERE, '..', '..', n) for n in ('prl_v14.tex', 'prl_v13.tex', 'prl_r4.tex', 'prl.tex')) if os.path.exists(p)), 'prl.tex')

def tex2sympy(body, Dsym):
    b = body
    b = b.replace(r'\bigl[', '(').replace(r'\bigr]', ')')
    b = re.sub(r'\\Delta_[13]', Dsym, b)
    b = re.sub(r'\\tfrac\{([^{}]*)\}\{([^{}]*)\}', r'((\1)/(\2))', b)
    b = re.sub(r'\^\{(\d+)\}', r'**\1', b)
    b = b.replace('Cr', 'C*r')
    return parse_expr(b, local_dict={'q': q, 'C': C, 'r': r, Dsym: sp.Symbol(Dsym)}, transformations=TR)

NAMES = {'a_0': 'a0', 'a_1': 'a1', 'b_0': 'b0', 'b_1': 'b1', 'e_{00}': 'e00', 'e_{01}': 'e01',
         'e_{10}': 'e10', 'e_{11}': 'e11', 'd_0^{(0)}': 'd0_0', 'd_0^{(1)}': 'd0_1',
         'd_1^{(0)}': 'd1_0', 'd_1^{(1)}': 'd1_1'}

def parse_vertex(tag):
    txt = open(TEX).read()
    start = txt.index(r'\paragraph{Data for $%s$.}' % tag)
    blk = txt[start: txt.index(r'\end{align*}', start)]
    Dtex = re.search(r'\\Delta_([13])=([^:$]*)\$?:', blk)
    Dsym = 'D' + Dtex.group(1)
    Dpoly = parse_expr(re.sub(r'(\d)q', r'\1*q', Dtex.group(2)).replace('^', '**'), local_dict={'q': q}, transformations=TR)
    out = {}
    for line in blk.splitlines():
        m = re.match(r'\s*([a-z]_[^&]*)&=(.*?)(,\\\\|\.)\s*$', line)
        if not m: continue
        name = NAMES[m.group(1).strip()]
        e = tex2sympy(m.group(2), Dsym)
        out[name] = sp.expand(e.subs(sp.Symbol(Dsym), Dpoly)) if False else e.subs(sp.Symbol(Dsym), Dpoly)
    return out, Dpoly

def parse_table():
    txt = open(TEX).read()
    start = txt.index(r'\label{tab:R}')
    blk = txt[txt.index(r'\hline', start): txt.index(r'\end{tabular}', start)]
    fr = lambda t: sp.Rational(*map(int, re.search(r'\\frac\{(\d+)\}\{(\d+)\}', t).groups())) * (-1 if t.strip().lstrip('$').startswith('-') else 1)
    V1, V3 = {}, {}
    for line in blk.splitlines()[1:]:
        if '&&' not in line: continue
        L, Rr = line.split('&&')
        for part, dst in ((L, V1), (Rr, V3)):
            cells = [c.strip() for c in part.split('&')]
            cells = [c for c in cells if c]
            if len(cells) < 3: continue
            jl = tuple(map(int, re.search(r'\((\d),(\d)\)', cells[0]).groups()))
            dst[jl] = (fr(cells[1]), fr(cells[2].replace('\\\\', '')))
    return V1, V3

if __name__ == '__main__':
    for t in ('V_1', 'V_3'):
        d, D = parse_vertex(t); print(t, D, len(d)); 
        for k, v in d.items(): print(' ', k, '=', v)
    V1, V3 = parse_table(); print(len(V1), len(V3)); print(V1); print(V3)
