"""Rigorous Arb enclosures for acceptance tests.

Keep python-flint arb balls through every comparison.  Convert to a printed
decimal only after the rounding direction has been certified by an Arb
comparison.  Do not use binary float mid/rad conversions as outward bounds.
"""
from fractions import Fraction as Fr
import math
from flint import arb, ctx

if not __debug__:
    raise RuntimeError('Verification must not run with Python optimization enabled')

if ctx.prec < 200:
    ctx.prec = 200

LOG2 = arb(2).log()


def aq(q):
    q = Fr(q)
    return arb(q.numerator) / q.denominator


def strictly_above(b):
    """An arb ball lying entirely strictly above every value in b."""
    u = b.upper()
    if u > b:
        return u
    # Inflate by a positive dyadic in case upper() is inconclusive.
    mid, rad, exp = map(int, b.mid_rad_10exp())
    bump = arb(1) / arb(10)**max(0, 8 - exp) if exp < 8 else arb(10)**(exp - 8)
    u = b + abs(b).upper() * arb(0) + bump
    u = (b + bump).upper()
    if not (u > b):
        raise ValueError('could not construct a strict upper enclosure')
    return u


def strictly_below(b):
    """An arb ball lying entirely strictly below every value in b."""
    lo = b.lower()
    if lo < b:
        return lo
    mid, rad, exp = map(int, b.mid_rad_10exp())
    bump = arb(1) / arb(10)**max(0, 8 - exp) if exp < 8 else arb(10)**(exp - 8)
    lo = (b - bump).lower()
    if not (lo < b):
        raise ValueError('could not construct a strict lower enclosure')
    return lo


def certified_max_upper(balls):
    """An arb ball strictly above every value in every input ball."""
    if not balls:
        raise ValueError('empty collection')
    u = strictly_above(balls[0])
    for b in balls[1:]:
        bu = strictly_above(b)
        if bu > u:
            u = bu
        elif u > bu:
            pass
        else:
            u = strictly_above(u.max(bu))
    return u


def certified_min_lower(balls):
    """An arb ball strictly below every value in every input ball."""
    if not balls:
        raise ValueError('empty collection')
    lo = strictly_below(balls[0])
    for b in balls[1:]:
        bl = strictly_below(b)
        if bl < lo:
            lo = bl
        elif lo < bl:
            pass
        else:
            lo = strictly_below(lo.min(bl))
    return lo


def definitely_gt(a, b):
    return bool(a > b)


def definitely_ge(a, b):
    return bool(a > b) or bool(a >= b)


def definitely_lt(a, b):
    return bool(a < b)


def hmin_lower_from_P_upper(P_up):
    """Certified lower-bound ball for -log2(P_guess), given P_up strictly above P_guess."""
    if not (P_up < 1):
        return arb(0)
    if not (P_up > 0):
        raise ValueError('P_up is not certified positive and less than one')
    return strictly_below(-(P_up.log() / LOG2))


def hmin_upper_from_P_lower(P_lo):
    """Certified upper-bound ball for -log2(P_guess), given P_lo strictly below P_guess."""
    if not (P_lo > 0):
        return None
    return strictly_above(-(P_lo.log() / LOG2))


def chsh_hmin_upper(S):
    """Certified upper bound on f(S) = -log2(1/2 + 1/2 sqrt(2-S^2/4)) for exact rational S.

    f is increasing in S on [2, 2*sqrt(2)].  An NPA upper bound on the recorded
    CHSH value may be passed in and still yields a valid upper bound on the
    CHSH min-entropy cap.
    """
    S = Fr(S)
    Sa = aq(S)
    inside = arb(2) - Sa * Sa / 4
    if inside < 0:
        raise ValueError('CHSH argument 2-S^2/4 is negative')
    # Enclosure of the exact guessing probability; take a strict lower bound
    # so that -log2 of it is a strict upper bound on f(S).
    P = arb(1) / 2 + inside.sqrt() / 2
    P_lo = strictly_below(P)
    if not (P_lo > 0):
        raise ValueError('CHSH guessing probability lower enclosure is not positive')
    return strictly_above(-(P_lo.log() / LOG2))


def decimal_floor(b, digits):
    """Largest integer k such that k / 10^digits is definitely < every value in b.

    Returns (k, scale) with scale = 10^digits.  Raises if a unique floor cannot
    be isolated.
    """
    if digits < 0:
        raise ValueError('digits must be nonnegative')
    scale = 10 ** digits
    s = strictly_below(b) * scale
    k = s.floor().unique_fmpz()
    if k is None:
        k = strictly_below(s).floor().unique_fmpz()
    if k is None:
        raise ValueError('could not isolate a certified decimal floor')
    k = int(k)
    val = arb(k) / scale
    while not (val < b):
        k -= 1
        val = arb(k) / scale
        if k < -10 * scale:
            raise ValueError('decimal floor search failed')
    # maximality: k+1 should not be definitely below b
    return k, scale


def decimal_ceil(b, digits):
    """Smallest integer k such that k / 10^digits is definitely > every value in b."""
    if digits < 0:
        raise ValueError('digits must be nonnegative')
    scale = 10 ** digits
    s = strictly_above(b) * scale
    k = s.ceil().unique_fmpz()
    if k is None:
        k = strictly_above(s).ceil().unique_fmpz()
    if k is None:
        raise ValueError('could not isolate a certified decimal ceiling')
    k = int(k)
    val = arb(k) / scale
    while not (val > b):
        k += 1
        val = arb(k) / scale
        if k > 10 * scale:
            raise ValueError('decimal ceiling search failed')
    return k, scale


def format_lower(b, digits):
    """Decimal string, rounded down, certified <= every value in b."""
    k, scale = decimal_floor(b, digits)
    return _fmt_fixed(k, digits)


def format_upper(b, digits):
    """Decimal string, rounded up, certified >= every value in b."""
    k, scale = decimal_ceil(b, digits)
    return _fmt_fixed(k, digits)


def _fmt_fixed(k, digits):
    sign = '-' if k < 0 else ''
    k = abs(k)
    scale = 10 ** digits
    ip, fp = divmod(k, scale)
    return sign + str(ip) + '.' + str(fp).zfill(digits)


def certifies_at_least(b, claimed):
    """True iff every value in b is definitely >= the exact rational claimed."""
    return definitely_ge(b, aq(claimed))


def certifies_greater(b, claimed):
    return definitely_gt(b, aq(claimed))


def certified_float_upper(b):
    """Python float guaranteed to lie strictly above every value in b.

    Zero balls may return 0.0.  After the decimal-to-binary conversion, the
    exact binary-rational value of the returned float is compared with the
    original ball, and the float is moved outward with math.nextafter until
    the comparison is certified.  Non-finite floats raise.
    """
    if b.is_zero():
        return 0.0
    x = float(format_upper(b, 17))
    while True:
        if not math.isfinite(x):
            raise ValueError('Non-finite float bound')
        if aq(Fr.from_float(x)) > b:
            return x
        x = math.nextafter(x, math.inf)


def certified_float_lower(b):
    """Python float guaranteed to lie strictly below every value in b.

    Zero balls may return 0.0.  The returned float itself is certified, not
    merely the decimal used to construct it.  Non-finite floats raise.
    """
    if b.is_zero():
        return 0.0
    x = float(format_lower(b, 17))
    while True:
        if not math.isfinite(x):
            raise ValueError('Non-finite float bound')
        if aq(Fr.from_float(x)) < b:
            return x
        x = math.nextafter(x, -math.inf)
