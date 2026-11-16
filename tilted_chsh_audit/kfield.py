"""Exact arithmetic in K = Q(s)(C, r), C^2 = 1 - s^2, r^2 = 1 + s^2.
Written for the independent audit (no code reused from anc/ or tilted/)."""
import sympy as sp
from sympy import QQ

s_sym = sp.Symbol('s')
FS = QQ.frac_field(s_sym)            # Q(s)
S = FS(s_sym)
ONE, ZERO = FS(1), FS(0)
C2 = 1 - S**2                         # C^2
R2 = 1 + S**2                         # r^2

class K:
    __slots__ = ('c',)                # coefficients of 1, C, r, Cr
    def __init__(self, a=0, b=0, c=0, d=0):
        self.c = tuple(x if not isinstance(x, (int, sp.Rational)) else FS(x) for x in (a, b, c, d))
        self.c = tuple(FS(x) if not hasattr(x, 'numer') else x for x in self.c)
    @staticmethod
    def lift(x):
        if isinstance(x, K): return x
        return K(FS(x) if not hasattr(x, 'numer') else x)
    def __add__(self, o):
        o = K.lift(o); return K(*[x + y for x, y in zip(self.c, o.c)])
    __radd__ = __add__
    def __neg__(self): return K(*[-x for x in self.c])
    def __sub__(self, o): return self + (-K.lift(o))
    def __rsub__(self, o): return K.lift(o) - self
    def __mul__(self, o):
        o = K.lift(o)
        a1, b1, c1, d1 = self.c; a2, b2, c2, d2 = o.c
        # basis 1, C, r, Cr ; C^2 = C2, r^2 = R2
        A = a1*a2 + b1*b2*C2 + c1*c2*R2 + d1*d2*C2*R2
        B = a1*b2 + b1*a2 + c1*d2*R2 + d1*c2*R2
        Cc = a1*c2 + c1*a2 + b1*d2*C2 + d1*b2*C2
        D = a1*d2 + d1*a2 + b1*c2 + c1*b2
        return K(A, B, Cc, D)
    __rmul__ = __mul__
    def conjC(self): a, b, c, d = self.c; return K(a, -b, c, -d)
    def conjr(self): a, b, c, d = self.c; return K(a, b, -c, -d)
    def inv(self):
        y = self * self.conjC()               # in Q(s)(r)
        assert y.c[1] == 0 and y.c[3] == 0
        z = y * y.conjr()                     # in Q(s)
        assert z.c[1] == 0 and z.c[2] == 0 and z.c[3] == 0
        if z.c[0] == 0: raise ZeroDivisionError('norm zero')
        return self.conjC() * y.conjr() * K(1 / z.c[0])
    def __truediv__(self, o): return self * K.lift(o).inv()
    def __rtruediv__(self, o): return K.lift(o) * self.inv()
    def __pow__(self, n):
        assert isinstance(n, int)
        if n < 0: return self.inv() ** (-n)
        out = K(1)
        for _ in range(n): out = out * self
        return out
    def iszero(self): return all(x == 0 for x in self.c)
    def __eq__(self, o): return (self - K.lift(o)).iszero()
    def norm(self):
        y = self * self.conjC(); z = y * y.conjr(); return z.c[0]
    def to_sympy(self, s=s_sym, C=sp.Symbol('C'), r=sp.Symbol('r')):
        out = 0
        for coef, mon in zip(self.c, (1, C, r, C*r)):
            out += FS.to_sympy(coef) * mon
        return out
    def __repr__(self): return 'K(%s)' % ', '.join(str(FS.to_sympy(x)) for x in self.c)

Ksym = {'C': K(0, 1), 'r': K(0, 0, 1), 's': K(S), 'q': K(S**2)}

def from_sympy(e):
    """Convert a sympy expression in s, q, C, r (rational ops only) to K."""
    if e.is_Symbol:
        return Ksym[e.name]
    if e.is_Rational:
        return K(FS(e))
    if e.is_Add:
        out = K(0)
        for t in e.args: out = out + from_sympy(t)
        return out
    if e.is_Mul:
        out = K(1)
        for t in e.args: out = out * from_sympy(t)
        return out
    if e.is_Pow:
        b, n = e.args
        assert n.is_Integer, e
        return from_sympy(b) ** int(n)
    raise ValueError('unsupported: %r' % e)

C_, r_, s_, q_ = Ksym['C'], Ksym['r'], Ksym['s'], Ksym['q']
