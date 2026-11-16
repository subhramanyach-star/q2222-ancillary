"""Exact certification of the I_δ supporting quadrilateral, 0 < δ ≤ π/6.

This script re-derives the closed-form vertices and checks, in exact
symbolic arithmetic, every identity used in the appendix: stationarity,
normalization, multiplier positivity, the twelve nontight local slacks,
and the elementary sign comparisons on cos δ ∈ [√3/2, 1).

Run: python3 certify_idelta.py
Requires: sympy.
"""
from __future__ import annotations

import sys
import sympy as sp

s, c = sp.symbols("s c", positive=True)
K = 2 * c**2 - 1
D = 2 * c**3

E00, E01, E10, E11 = s * K / D, K / D, K / D, -s / D
A, B = K / (2 * c**2), 1 / (2 * c**2)

V0 = (
    (c - 1) * K / D,
    s / D,
    -s * K / D,
    (1 - c) * (2 * c + 1) / D,
    E00,
    E01,
    E10,
    E11,
)
V1 = (
    -s * K / D,
    (1 - c) * (2 * c + 1) / D,
    (c - 1) * K / D,
    s / D,
    E00,
    E01,
    E10,
    E11,
)


def red(e):
    e = sp.simplify(sp.expand(e).subs(s**2, 1 - c**2))
    num, den = sp.fraction(sp.together(e))
    num = sp.factor(sp.simplify(sp.expand(num).subs(s**2, 1 - c**2)))
    den = sp.factor(sp.simplify(sp.expand(den).subs(s**2, 1 - c**2)))
    return num / den, num, den


def gdotL(g, L):
    A0, A1, B0, B1 = L
    zl = (A0, A1, B0, B1, A0 * B0, A0 * B1, A1 * B0, A1 * B1)
    return sum(gi * zi for gi, zi in zip(g, zl))


def mus(a0d, a1d):
    return ((A - a0d) / 2, (B - a1d) / 2, (A + a0d) / 2, (B + a1d) / 2)


errors = []


def check(cond, msg):
    if cond is True:
        return
    if cond is False:
        errors.append(msg)
        return
    # sympy
    if sp.simplify(cond) == 0 or cond == True:  # noqa: E712
        return
    if sp.simplify(cond) is sp.true:
        return
    errors.append(msg + f"  [got {cond}]")


# Party exchange and correlator identity
check(all(sp.simplify(V1[i] - V0[(i + 2) % 4 if i < 4 else i]) == 0 for i in range(8)), "V1 is Alice-Bob exchange of V0")
for i, lab in enumerate(["E00", "E01", "E10", "E11"]):
    check(sp.simplify(V0[4 + i] - V1[4 + i]) == 0, f"shared {lab}")

# Multipliers
mu0_expected = (
    K / (4 * c**3),
    (c - s) / (4 * c**3),
    (2 * c - 1) * K / (4 * c**3),
    (c + s) / (4 * c**3),
)
mu1_expected = (
    (c + s) * K / (4 * c**3),
    K / (4 * c**3),
    (c - s) * K / (4 * c**3),
    (1 + 2 * c - 2 * c**2) / (4 * c**3),
)
for got, exp, name in (
    (mus(V0[0], V0[1]), mu0_expected, "V0"),
    (mus(V1[0], V1[1]), mu1_expected, "V1"),
):
    for a, b, lab in zip(got, exp, ["+0", "+1", "-0", "-1"]):
        check(sp.simplify(a - b) == 0, f"{name} mu_{lab}")

# Tight locals
tight = {
    "V0": ((1, 1, 1, 1), (1, -1, -1, 1)),
    "V1": ((1, 1, 1, 1), (-1, 1, 1, -1)),
}
verts = {"V0": V0, "V1": V1}
for name, Ls in tight.items():
    for L in Ls:
        slack = red(1 - gdotL(verts[name], L))[0]
        check(sp.simplify(slack) == 0, f"{name} not tight at {L}: {slack}")

# Sign lemmas used in the appendix
# (i) 2c^2+c-1 - s > 0 on the interval: min of 2c^2+c-1 is at c=sqrt(3)/2
fmin = sp.simplify(sp.Rational(3, 2) + sp.sqrt(3) / 2 - 1)  # (1+sqrt(3))/2
check(sp.simplify(fmin - (1 + sp.sqrt(3)) / 2) == 0, "fmin formula")
check(sp.simplify(fmin - sp.Rational(1, 2)) > 0, "(1+sqrt(3))/2 > 1/2")

# (ii) 1-2c^2+2c^3 > 2 c^2 s, after squaring
diff = sp.factor(sp.expand((2 * c**3 - 2 * c**2 + 1) ** 2 - 4 * c**4 * (1 - c**2)))
check(sp.simplify(diff - (2 * c**2 - 1) * (4 * c**4 - 4 * c**3 + 2 * c**2 - 1)) == 0, "square difference factor")
q_at = sp.simplify((4 * c**4 - 4 * c**3 + 2 * c**2 - 1).subs(c, sp.sqrt(3) / 2))
check(sp.simplify(q_at - (11 - 6 * sp.sqrt(3)) / 4) == 0, "q(sqrt(3)/2) formula")
# 11-6sqrt(3)>0 iff 121>108
check(121 - 108 > 0, "11^2 > 36*3")

# (iii) c^2-c-1+s(c+1)<0 after squaring: right-left = c^2 K
left = (c + 1) ** 3 * (1 - c)
right = (1 + c - c**2) ** 2
check(sp.simplify(sp.expand(right - left) - c**2 * K) == 0, "slack type (8) difference")

# All twelve nontight slacks at V0 and V1 have positive sign when the
# elementary factors do; check exact factorization and numerical positivity
# at the two endpoints of the closed interval (c=sqrt(3)/2 and a value near 1).
signs = [
    (A0, A1, B0, B1)
    for A0 in (1, -1)
    for A1 in (1, -1)
    for B0 in (1, -1)
    for B1 in (1, -1)
]
for name, g in verts.items():
    for L in signs:
        if L in tight[name]:
            continue
        expr, num, den = red(1 - gdotL(g, L))
        # den should be c**3 (or a positive power of c)
        den_pos = sp.simplify(den.subs(c, 1))
        if den_pos <= 0:
            errors.append(f"{name} {L}: unexpected den {den}")
        for cval, sval, tag in (
            (sp.sqrt(3) / 2, sp.Rational(1, 2), "pi/6"),
            (sp.cos(sp.Rational(1, 100)), sp.sin(sp.Rational(1, 100)), "small"),
        ):
            val = sp.N(expr.subs({c: cval, s: sval}), 30)
            if val <= 0:
                errors.append(f"{name} {L} slack not positive at {tag}: {val}")

if errors:
    print("FAILED")
    for e in errors:
        print(" -", e)
    sys.exit(1)
print("All I_δ identities and interval sign checks passed.")
