"""Compile sympy expressions in (s, C, r) (and optionally other symbols) into rigorous arb evaluators.
Rationals are converted exactly (fmpq), sqrt -> arb.sqrt; no floating-point constants enter."""
import sympy as sp
import flint
from flint import arb, fmpq

def compile_expr(expr, syms):
    idx = {sy: i for i, sy in enumerate(syms)}
    def build(e):
        if e.is_Symbol:
            i = idx[e]; return lambda v: v[i]
        if e.is_Integer:
            k = int(e); return lambda v: arb(k)
        if e.is_Rational:
            q = fmpq(int(e.p), int(e.q)); return lambda v: arb(q)
        if e.is_Add:
            fs = [build(a) for a in e.args]
            def f(v):
                acc = fs[0](v)
                for g in fs[1:]: acc = acc + g(v)
                return acc
            return f
        if e.is_Mul:
            fs = [build(a) for a in e.args]
            def f(v):
                acc = fs[0](v)
                for g in fs[1:]: acc = acc * g(v)
                return acc
            return f
        if e.is_Pow:
            b = build(e.base); ex = e.exp
            if ex.is_Integer:
                n = int(ex)
                if n >= 0: return lambda v: b(v) ** n
                return lambda v: 1 / (b(v) ** (-n))
            if ex == sp.Rational(1, 2): return lambda v: b(v).sqrt()
            if ex == sp.Rational(-1, 2): return lambda v: 1 / b(v).sqrt()
            raise ValueError(ex)
        raise ValueError("unsupported: %s" % e)
    return build(sp.sympify(expr))

S_, C_, R_ = sp.symbols('s C r', positive=True)
def theta_vars(th):
    """th: arb (ball). Returns (s, C, r) enclosures, s = sin 2th, C = cos 2th, r = sqrt(1+s^2)."""
    s = (2 * th).sin(); C = (2 * th).cos()
    return [s, C, (1 + s * s).sqrt()]
