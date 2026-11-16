"""Theorem 2 objects rebuilt from the mathematics as printed in prl.tex (Sec. 'Tilted CHSH',
Appendix C, Table III).  Everything is an exact element of K = Q(s)(C,r), s = sin 2theta."""
import itertools, sympy as sp
from kfield import K, FS, S, C_, r_, s_, q_, from_sympy
from parse_tex import parse_vertex, parse_table

q = q_
# --- main text: p_theta, c0, m -------------------------------------------------------------
P = [C_, K(0), C_ / r_, C_ / r_, 1 / r_, 1 / r_, q / r_, -q / r_]
c0 = (1 - q) / (1 + q)
m_printed = [r_*r_/4 * x for x in (1 + C_, -q, 1 - C_, -q)]          # order (+|0),(+|1),(-|0),(-|1)
AX = [(+1, 0), (+1, 1), (-1, 0), (-1, 1)]

def alpha(z, x): return z[x]
def beta(z): return (z[2], z[3])
def Ex(z, x): return (z[4 + 2*x], z[5 + 2*x])
def t_(z, a, x): return 1 + a*z[x]
def w_(z, a, x):
    b0, b1 = beta(z); E0, E1 = Ex(z, x); return (b0 + a*E0, b1 + a*E1)

def Gmul(c, u): return (u[0] + c*u[1], c*u[0] + u[1])
def Gnorm2(c, u): return u[0]*u[0] + u[1]*u[1] + 2*c*u[0]*u[1]
def Ginv_mul(c, w):
    d = 1 - c*c; return ((w[0] - c*w[1]) / d, (-c*w[0] + w[1]) / d)

def dot(g, z): return sum((gi*zi for gi, zi in zip(g, z)), K(0))

# deterministic points L = (A0,A1,B0,B1)
def Lvec(A0, A1, B0, B1): return [K(A0), K(A1), K(B0), K(B1), K(A0*B0), K(A0*B1), K(A1*B0), K(A1*B1)]
LPTS = list(itertools.product((1, -1), repeat=4))

def gamma_and_v():
    v = {}; gam = {}
    for (a, x) in AX:
        vv = Ginv_mul(c0, w_(P, a, x)); vv = (vv[0] / t_(P, a, x), vv[1] / t_(P, a, x))
        v[(a, x)] = vv
        g = [K(0)]*8; g[x] = K(-a); g[2] = vv[0]; g[3] = vv[1]; g[4+2*x] = a*vv[0]; g[5+2*x] = a*vv[1]
        gam[(a, x)] = g
    return gam, v
GAM, VV = gamma_and_v()

def zdot():
    d = 1 - c0*c0; z = [K(0)]*8
    z[3] = (P[2] - c0*P[3]) / d
    z[5] = (P[4] - c0*P[5]) / d
    z[7] = (P[6] - c0*P[7]) / d
    return z

# --- vertices --------------------------------------------------------------------------------
def vertex_from_tex(tag):
    d, Dpoly = parse_vertex(tag)
    conv = {k: from_sympy(v) for k, v in d.items()}
    g = [conv[k] for k in ('a0', 'a1', 'b0', 'b1', 'e00', 'e01', 'e10', 'e11')]
    return dict(g=g, d0=(conv['d0_0'], conv['d0_1']), d1=(conv['d1_0'], conv['d1_1']), Delta=Dpoly)

def V2V4():
    D = r_ + 1
    a0 = (r_ + C_) / (2*D); be = (C_ - 1) / (2*D); ka = (C_ + r_ + 2*q) / (4*D); nu = (r_ + 2 - C_) / (4*D)
    V2 = [a0, K(0), be, be, ka, ka, nu, -nu]
    # V4: kappa,nu with C -> -C ; a0 = (C-r)/(2D), beta = (C+1)/(2D)
    a0_4 = (C_ - r_) / (2*D); be4 = (C_ + 1) / (2*D); ka4 = (-C_ + r_ + 2*q) / (4*D); nu4 = (r_ + 2 + C_) / (4*D)
    V4 = [a0_4, K(0), be4, be4, ka4, ka4, nu4, -nu4]
    return dict(g=V2, a0=a0, beta=be, kappa=ka, nu=nu), dict(g=V4, a0=a0_4, beta=be4, kappa=ka4, nu=nu4)

def sigma_g(g):
    a0, a1, b0, b1, e00, e01, e10, e11 = g
    return [a0, -a1, b1, b0, e01, e00, -e11, -e10]
def sigma_L(L):
    A0, A1, B0, B1 = L; return (A0, -A1, B1, B0)

FACETS = {'V1': [(1, 1, 1, 1), (1, -1, -1, 1)], 'V2': [(1, -1, -1, 1), (1, 1, 1, -1)],
          'V3': [(1, -1, 1, 1), (-1, 1, 1, -1)], 'V4': [(-1, 1, 1, -1), (-1, -1, -1, 1)]}
FACETS['V1s'] = [sigma_L(L) for L in FACETS['V1']]
FACETS['V3s'] = [sigma_L(L) for L in FACETS['V3']]

TABLE = parse_table()
KL = {'V1': (1, 2), 'V3': (2, 3)}
SIG = 10*s_ - K(sp.Rational(83, 10))          # sigma = (sin2theta - 83/100)/(1/10)

def Rpoly_at(tag, cval):
    """sum_{j,l} R_jl c^j sigma^l at exact c (K), returns 2-vector in K"""
    R = TABLE[0] if tag == 'V1' else TABLE[1]
    out = [K(0), K(0)]
    for (j, l), (R0, R1) in R.items():
        mon = (K.lift(cval) ** j) * (SIG ** l)
        out[0] = out[0] + K(R0) * mon; out[1] = out[1] + K(R1) * mon
    return out

def delta_at(V, tag, cval):
    """certificate (eq. V13cert) at exact c (an int / K)."""
    cv = K.lift(cval); h = cv - c0
    Rv = Rpoly_at(tag, cv)
    return tuple(V['d0'][i] + V['d1'][i]*h + h*h*Rv[i] for i in range(2))
