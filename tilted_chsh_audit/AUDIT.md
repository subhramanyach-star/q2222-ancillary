# Independent audit of Theorem 2 (tilted-CHSH normal cones), 21 Oct 2026 (HKT)

Written from the mathematics as printed in the paper (Theorem 1, Sec. Tilted CHSH, App. C, Table III); the source file is
`$PRL_TEX` if set, else `../../prl_v13.tex`, `../../prl_r4.tex` or `../../prl.tex`. No code from `../tilted_chsh/` or earlier scripts was imported or
copied. The JSON files in `../tilted_chsh/` are read only by `check_exact.py`, part (b), for comparison.

Parametrisation: s = sin 2theta in place of theta. On Theta this map is a bijection onto
[sin(pi/4), sin(3pi/8)]; C = sqrt(1-s^2), r = sqrt(1+s^2). Everything is exact in K = Q(s)(C, r) or uses
outward-rounded Arb (python-flint 0.9, 128-bit).

| file | content | runtime |
|---|---|---|
| `kfield.py` | exact arithmetic in Q(s)(C,r) | - |
| `parse_tex.py` | parses the V1/V3 data and Table III **directly from prl.tex** | - |
| `model.py` | p_theta, c0, m, generators, vertices, sigma, certificate | - |
| `kfun.py` | Arb evaluation, exact sign decision (norm root counting) | - |
| `check_exact.py` -> `check_exact.log` | checks (a), (b), (c), (e) and the V2/V4 hand proof: 303 PASS, 0 FAIL | 227 s |
| `bnb.py`, `run_full.sh` -> `bnb_V1.jsonl`, `bnb_V3.jsonl` (stdout duplicate `.log` files omitted) | checks (d), (f): 2-D branch-and-bound | V1 162 s, V3 205 s wall, ~65 MB |
| `summarize_bnb.py` | coverage/tiling summary of the B&B output | - |
| `check_misc.py` -> `check_misc.log` | exposure ingredients (exact); sigma symmetry; primal SLSQP spot check | 112 s |
| `check_hexrange2.py` -> `check_hexrange2.log` | 108 hexagon inequalities on Theta and [pi/16, 0.2022pi] | 12 s |
| `spot_v24.py` | float spot check of the V2/V4 hand certificates | 5 s |
| `endpoint_lemma.tex/.pdf` | precise endpoint lemma, hypotheses (H1), (H2), proof | - |

Commands: `python3 check_exact.py`, `./run_full.sh` (checkpointed: re-running skips finished strips),
`python3 summarize_bnb.py V1 V3`, `python3 check_misc.py`, `python3 check_hexrange2.py`.

B&B tests per (c-box I, s-strip S, sign pair):
- T0: enclosure of F > 0. Radicands are clamped at 0, so zero radicands are allowed.
- T1: mean value in c.
- T2: F_cc > 0 on hull(I, c0(S)), using the exact F = F_c = 0 at c0.
- T3: monotone in c toward c = +-1.
- T3b: second-order bound from c = +-1, using the exact F(+-1, s) >= 0 from the endpoint lemma.

T1, T2, T3 and T3b each require certified A, B > 0 on the box. Coverage: 4096 exact rational s-strips tile
[floor(sin(pi/4)), ceil(sin(3pi/8))] at 1e-30 resolution, which slightly over-covers Theta. Within each strip,
sub-strips and c-leaves tile the strip and [-1,1]; this is checked in exact rational arithmetic.
Quadratic forms are evaluated in the G_c-orthogonal basis, ||u||^2 = (1+c)/2 (u.f+)^2 + (1-c)/2 (u.f-)^2,
to avoid cancellation. This matters because the V3 radicand at c = -1 is only ~1.25e-4.
