"""Exact rational loss-map identities for Bell functionals in the 2222 scenario.

Solver-free: this module uses only the standard library (fractions).  It is the
common source of loss_functional_F and G_CHSH_F for both certificate generation
and stand-alone verification, so that verification does not import SciPy, CVXPY
or search modules through cert_eta.py.
"""

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

from fractions import Fraction as F

G_CHSH_F = [F(x) for x in (0, 0, 0, 0, 1, 1, 1, -1)]


def loss_functional_F(g, eta):
    """Exact identity: g · z_eta = c0 + g_eta · z_ideal for the binned behaviour.

    eta and the entries of g may be int or Fraction; the returned c0 and the
    eight-vector g_eta are Fraction.
    """
    eta = F(eta)
    g = [F(x) for x in g]
    if len(g) != 8:
        raise ValueError('Bell functional must have eight coefficients')
    E = g[4:]
    c0 = (1 - eta) * (g[0] + g[1] + g[2] + g[3]) + (1 - eta)**2 * sum(E)
    ga = [eta * g[x] + eta * (1 - eta) * (E[2 * x] + E[2 * x + 1]) for x in (0, 1)]
    gb = [eta * g[2 + y] + eta * (1 - eta) * (E[y] + E[2 + y]) for y in (0, 1)]
    return c0, ga + gb + [eta**2 * e for e in E]


def class2_functional(eta):
    """Loss-transformed three-negative CHSH variant, signs (-,-,-,+), no-click as +1."""
    eta = F(eta)
    sv = [-1, -1, -1, 1]
    sa = [sv[0] + sv[1], sv[2] + sv[3]]
    sb = [sv[0] + sv[2], sv[1] + sv[3]]
    c0 = (1 - eta)**2 * sum(sv)
    g = [eta * (1 - eta) * sa[0], eta * (1 - eta) * sa[1],
         eta * (1 - eta) * sb[0], eta * (1 - eta) * sb[1]] + [eta**2 * x for x in sv]
    return c0, g
