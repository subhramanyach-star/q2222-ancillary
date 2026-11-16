"""Certified enclosures of the Gigena et al. canonical CHSH maximum under equal-efficiency loss.

For eta in (2/3, 1], Theorem smax of the manuscript states
    S_max(eta) = 2(1-eta)^2 + eta^2 c_Q(alpha, alpha),
    alpha = 2(1-eta)/eta,
where c_Q(alpha, alpha) is the largest real root of the degree-4 polynomial
of Ref. [Gigena2025, Thm. 1, Eq. (18)]:

    f(lambda) = lambda^4 + (4-alpha^2) lambda^3 + (11/4 alpha^4 - 12 alpha^2 - 4) lambda^2
              + (2 alpha^6 - alpha^4 - 20 alpha^2 - 32) lambda
              + 5 alpha^6 - 21 alpha^4 + 16 alpha^2 - 32.

Coefficients are exact rationals when eta is rational.  The largest real root is
isolated by a Sturm sequence and bisection in exact rationals, so the returned
interval [L, U] satisfies L < c_Q < U with no floating-point rounding.
"""
from fractions import Fraction as Fr
from flint import fmpq, fmpq_poly

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')


def _fr_to_fmpq(x):
    x = Fr(x)
    return fmpq(x.numerator, x.denominator)


def gigena_quartic(alpha):
    """Exact Eq. (18) as an fmpq_poly, low-degree coefficient first."""
    a = _fr_to_fmpq(alpha)
    a2 = a * a
    a4 = a2 * a2
    a6 = a4 * a2
    c0 = 5 * a6 - 21 * a4 + 16 * a2 - 32
    c1 = 2 * a6 - a4 - 20 * a2 - 32
    c2 = fmpq(11, 4) * a4 - 12 * a2 - 4
    c3 = 4 - a2
    c4 = fmpq(1)
    return fmpq_poly([c0, c1, c2, c3, c4])


def _sturm_sequence(p):
    seq = [p, p.derivative()]
    while not seq[-1].is_zero() and seq[-1].degree() >= 0:
        r = seq[-2] % seq[-1]
        if r.is_zero():
            break
        seq.append(-r)
        if len(seq) > p.degree() + 4:
            raise ValueError('Sturm sequence did not terminate')
    return seq


def _sign(q):
    if q == 0:
        return 0
    return 1 if q > 0 else -1


def _eval_fmpq(p, x):
    return p(_fr_to_fmpq(x))


def _sign_changes(seq, x, end='finite'):
    """Number of sign changes of the Sturm sequence at a finite point or at +/- infinity."""
    signs = []
    if end == 'plus_inf':
        for s in seq:
            if s.is_zero():
                continue
            lc = s.leading_coefficient()
            signs.append(_sign(lc))
    elif end == 'minus_inf':
        for s in seq:
            if s.is_zero():
                continue
            lc = s.leading_coefficient()
            deg = s.degree()
            signs.append(_sign(lc) * ((-1) ** deg))
    else:
        xf = _fr_to_fmpq(x)
        for s in seq:
            v = s(xf)
            sv = _sign(v)
            if sv != 0:
                signs.append(sv)
    n = 0
    for a, b in zip(signs, signs[1:]):
        if a != b:
            n += 1
    return n


def n_real_roots(p, a=None, b=None):
    """Number of distinct real roots in (a, b); omit an endpoint for an infinite ray."""
    seq = _sturm_sequence(p)
    if a is None:
        sa = _sign_changes(seq, None, 'minus_inf')
    else:
        if _eval_fmpq(p, a) == 0:
            raise ValueError('Sturm count requires a non-root left endpoint')
        sa = _sign_changes(seq, a, 'finite')
    if b is None:
        sb = _sign_changes(seq, None, 'plus_inf')
    else:
        if _eval_fmpq(p, b) == 0:
            raise ValueError('Sturm count requires a non-root right endpoint')
        sb = _sign_changes(seq, b, 'finite')
    return sa - sb


def isolate_largest_root(p, lo=Fr(0), hi=Fr(8), max_den_bits=80):
    """Return Fractions (L, U) with exactly one real root of p in (L, U), the largest one.

    Proves there is no real root in [U, +infty) and exactly one in (L, U).
    """
    lo, hi = Fr(lo), Fr(hi)
    n_hi = n_real_roots(p, hi, None)
    if n_hi != 0:
        # widen
        hi = Fr(hi * 2)
        if n_real_roots(p, hi, None) != 0:
            raise ValueError('largest root is not below the search ceiling')
    n_int = n_real_roots(p, lo, hi)
    if n_int < 1:
        raise ValueError('no real root in the search interval')
    # shrink the left end until a single root remains in (L, hi)
    L, U = lo, hi
    while n_real_roots(p, L, U) > 1:
        mid = (L + U) / 2
        if _eval_fmpq(p, mid) == 0:
            mid = (L + 2 * U) / 3
        if n_real_roots(p, mid, U) >= 1:
            L = mid
        else:
            U = mid
        if (U - L) < Fr(1, 10 ** 12) and n_real_roots(p, L, U) > 1:
            raise ValueError('failed to isolate the largest root')
    # now exactly one root in (L, U); bisect to the requested tightness
    # using Sturm counts at every step rather than a sign-change assumption
    target = Fr(1, 1 << max_den_bits)
    while (U - L) > target:
        mid = (L + U) / 2
        if _eval_fmpq(p, mid) == 0:
            # mid is the root, which we are not allowed to return as an endpoint
            L, U = (L + mid) / 2, (mid + U) / 2
            continue
        n_right = n_real_roots(p, mid, U)
        if n_right == 1:
            if n_real_roots(p, L, mid) != 0:
                raise ValueError('bisection lost the unique largest-root isolation')
            L = mid
        elif n_right == 0:
            if n_real_roots(p, L, mid) != 1:
                raise ValueError('bisection lost the unique root')
            U = mid
        else:
            raise ValueError('bisection lost the unique largest-root isolation')
    if n_real_roots(p, L, U) != 1 or n_real_roots(p, U, None) != 0:
        raise ValueError('bisection lost the unique largest-root isolation')
    return L, U


def cQ_bounds(alpha, max_den_bits=80):
    alpha = Fr(alpha)
    if not (0 <= alpha < 1):
        raise ValueError('alpha = 2(1-eta)/eta must lie in [0, 1)')
    p = gigena_quartic(alpha)
    # For alpha=0 the root is 2*sqrt(2) in (2, 3); for alpha -> 1- it stays below 4.
    return isolate_largest_root(p, lo=Fr(2), hi=Fr(5), max_den_bits=max_den_bits)


def S_max_bounds(eta, max_den_bits=80):
    """Exact-rational bounds L < S_max(eta) < U."""
    eta = Fr(eta)
    if not (Fr(2, 3) < eta <= 1):
        raise ValueError('eta must lie in (2/3, 1]')
    alpha = 2 * (1 - eta) / eta
    cL, cU = cQ_bounds(alpha, max_den_bits=max_den_bits)
    pre = 2 * (1 - eta) ** 2
    w = eta * eta
    return pre + w * cL, pre + w * cU
