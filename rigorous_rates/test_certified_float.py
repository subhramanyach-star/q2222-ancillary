#!/usr/bin/env python3
"""Unit check that certified_float_{upper,lower} do not round the wrong way in binary.

Counterexamples: b = 1 + 2^{-60} and b = 1 - 2^{-60}.  A 17-place decimal of
the first converts to the float 1.0, which lies below the ball; the second
converts to 1.0, which lies above the ball.  The helpers must move outward
with math.nextafter until Fraction.from_float of the returned value is a
certified bound.
"""
import math
import sys
from fractions import Fraction as Fr
from flint import arb

from arb_bounds import aq, certified_float_upper, certified_float_lower, format_upper, format_lower


def _check_pair():
    b_up = arb(1) + arb(1) / arb(2)**60
    naive_up = float(format_upper(b_up, 17))
    if not (aq(Fr.from_float(naive_up)) < b_up):
        raise AssertionError('expected naive upper conversion of 1+2^{-60} to fall below the ball')
    x_up = certified_float_upper(b_up)
    if not math.isfinite(x_up):
        raise AssertionError('non-finite certified upper float')
    if not (aq(Fr.from_float(x_up)) > b_up):
        raise AssertionError('certified_float_upper(1+2^{-60}) is not strictly above the ball')

    b_lo = arb(1) - arb(1) / arb(2)**60
    naive_lo = float(format_lower(b_lo, 17))
    if not (aq(Fr.from_float(naive_lo)) > b_lo):
        raise AssertionError('expected naive lower conversion of 1-2^{-60} to lie above the ball')
    x_lo = certified_float_lower(b_lo)
    if not math.isfinite(x_lo):
        raise AssertionError('non-finite certified lower float')
    if not (aq(Fr.from_float(x_lo)) < b_lo):
        raise AssertionError('certified_float_lower(1-2^{-60}) is not strictly below the ball')

    z = arb(0)
    if certified_float_upper(z) != 0.0 or certified_float_lower(z) != 0.0:
        raise AssertionError('zero balls must return 0.0')


if __name__ == '__main__':
    _check_pair()
    print('certified_float unit check: ALL OK')
    sys.exit(0)
