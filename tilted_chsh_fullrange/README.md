# tilted_chsh_fullrange/: Theorem 2 on the whole tilted-CHSH family, theta in (0, pi/4]

Claim: N_Q(p_theta) = C_theta for every theta in (0, pi/4]. The cone is a hexagon on (0, theta*] and an octagon on
(theta*, pi/4]. Here theta* = 0.2022405998 pi, and s*^2 = sin^2(2 theta*) is the real root of x^3 - 8x^2 + 24x - 16.
At pi/4 the octagon equals the Tsirelson octagon of Theorem 3(a) (pi4_check.log).
Notation: theta = pi x, s = sin 2theta, C = cos 2theta, r = sqrt(1+s^2), t = s^2. **The scripts' t is the paper's q = sin^2 2theta** (prl_v13, App. D, "Small theta").
`fullrange/REPORT.md` gives the full proof structure. `audit/AUDIT.md` is the independent audit, with its verdict and caveats.

The scripts load data by relative path, so run each command from inside its own directory.
They need `thm3_data.json` (shipped in `fullrange/`; identical to `../tilted_chsh/thm3_data.json`, md5 179df94e...).
Software is as in the top-level README: sympy, python-flint, mpmath, numpy. Arb working precision is 120 bits by default (`--prec`), as stated in App. D of prl_v13.

## fullrange/ (theta-strip certificates and exact combinatorics)
Exact checks (each takes seconds to about 2 min):
* `python3 comb.py`: polygon combinatorics. Every mu > 0 and every non-facet margin is decided exactly over Q(s)(C,r) on (0,s*) and (s*,1). Output: `comb.log`.
* `python3 v24_conditions.py`: V2 and V4 side conditions, including A2 > 0.
* `python3 v24_datacheck.py` (added for v13): the V2, V4 data in `thm3_data.json` equal the closed forms printed in App. D, exactly in Q(s)(C,r). (The archived `v24_conditions.log` of the earlier bundle ended with two "data matches closed form" lines that the shipped `v24_conditions.py` does not print; they came from an unshipped step. `v24_conditions.log` was regenerated from the shipped script, and that check now lives in `v24_datacheck.py` / `v24_datacheck.log`.)
* `python3 sigma_check.py`, `pi4_check.py`, `plane_check.py`, `exact_c0_all.py`: these check the sigma images, the pi/4 = Tsirelson octagon, the plane/stationarity identities and the exact forced data at c0.

Arb branch and bound: `python3 vgen.py CERT.json [--strips j0 j1] [--jobs N]`.

| cert | vertex | x = theta/pi | strips | boxes | log |
|---|---|---|---|---|---|
| c_Wa.json | Wa | [632/3125, 1/4] | 24 | 71114 | v_Wa_0_19_final.log, v_Wa_20_23.log |
| c_Wb.json | Wb | [632/3125, 1/4] | 24 | 43190 | v_Wb.log |
| c_V1h.json | V1 | [3/16, 809/4000] | 8 | 25192 | v_V1h.log |
| c_V3h.json | V3 | [3/16, 1/4] | 24 | 81099 | v_V3h.log |
| c_V1l.json | V1 | [1/16, 1/8] | 8 | 764840 | v_V1l.log |
| c_V3l.json | V3 | [1/16, 1/8] (mid-pin at c=-1) | 16 | 433681 | v_V3l.log |
| c_V1ll.json, c_V3ll_par.json | V1, V3 | [0.03, 1/16] (redundant second cover) | 16 | | v_V1ll.log, v_V3ll.log |

[1/8, 3/16] is covered by `../tilted_chsh/` (certs_v13.json). Wa and Wb use an endpoint pin at c=+1.
Strips take from seconds to a few minutes each. For example, `python3 vgen.py c_Wa.json --strips 23 24` takes about 5 s.

Negative controls must report `ALL OK = False`:
* `c_Wa_neg*.json` (v_neg*.log);
* `c_Wa_negpin_[a-d].json`, the pin path;
* `c_V3l_negmid_[01].json`, the mid-pin path.

`contain_pin.py` and `test_pin.py` are unit tests of the pin enclosure. `neg_*_design.py` generated the negative controls.

## smalltheta/ (uniform in t = sin^2 2theta, t in [0, 3/20], overlapping [pi/16, ...])
* V1, Alice-overlap slicing (party swap; exact swapped data `swap_data.json` from `swapdata.py`):
  * `python3 sgen.py c_V1s.json` covers t in [0, 1/63] and takes about 10 s;
  * `python3 sgen.py c_V1s_ext.json` covers [1/63, 3/20];
  * logs: `s_V1s_final.log`, `s_V1s_ext_final.log`;
  * negative controls: `c_V1s_neg[a-e].json` (s_neg*.log).
* V3, Bob-overlap slicing with exact division by t and Taylor coefficients at t=0 (`hats_V3_V1.json` from `prep.py`):
  * `python3 tgen3.py c_V3t.json [--tstrips j0 j1]` covers t in [0, 1/63] (t3_V3t.log, t_V3t_*.log);
  * `python3 tgen3.py c_V3m.json` covers [1/63, 3/20] (t3_V3m.log);
  * negative controls: `c_V3t_neg{ep,E*,F*}.json` (t2_negep.log, t3_negF*.log).
* `test_sgen.py`, `test_tgen3.py` (with `g3.py`) and `negE*.py` are float cross-checks and control generators. They are not acceptance tests.

## audit/
* `corner_check.py`: P = M = 0 exactly at (c,t) = (-1,0), together with the corner Taylor coefficients.
* `swap_uniqueness.py`: exact uniqueness of the swapped slice data.
* `tetra_arb.py` and `tetra_check.py`: the tetrahedron at theta*.
* `neg_T4_design.py`: wrote `c_V3t_negT4_[01].json`, which must fail under `cd ../smalltheta; python3 tgen3.py ../audit/c_V3t_negT4_0.json`.
* `test_tgen3_ep.py`, plus the spot-rerun logs `a_*.log` and `n_*.log` written by `run_audit.sh`.

## Changes relative to the working copies (paths only)
Absolute paths in the research copies were replaced with relative ones. No logic was changed. As a result, the `.py`/`.sh`
entries in `fullrange/HASHES.md5` and `smalltheta/HASHES_smalltheta.md5` no longer match for these files:
* `comb.py`, `vgen.py`, `exact_c0_all.py`, `pi4_check.py`, `plane_check.py`, `sigma_check.py`, `g3.py`, `prep.py`, `swapdata.py`, `tgen3.py`;
* the audit scripts.

`PATCHES.md5` lists the old and new hashes. All JSON certificates and data are byte-identical to the originals (`JSON.md5`).
`tgen3.py` used to read the cert field `hats` as an absolute path. It now reads its basename from the working directory.
The design-only scripts `cert.py`, `poly.py` and `branches.py` depend on unshipped research modules and are omitted, as are pickles.

## Earlier bundle smoke test (21 Oct 2026, HKT, from a fresh unpack of anc_submission_v7.tar.gz; the v13 smoke test is in ../SMOKETEST.log)
Each command below matched its archived output:
* `pi4_check.py`: ALL MATCH O: True;
* `sigma_check.py`: all images exact;
* `vgen.py c_Wa.json --strips 23 24`: ALL OK, 1667 boxes;
* `vgen.py c_Wa_negpin_a.json`: ALL OK = False, as required;
* `sgen.py c_V1s.json`: ALL OK;
* `tgen3.py c_V3t.json --tstrips 0 1`: ALL OK, 33758 boxes;
* `tgen3.py c_V3t_negep.json` and `../audit/c_V3t_negT4_0.json`: both False, as required;
* `audit/corner_check.py`: P and M vanish at the corner;
* `audit/swap_uniqueness.py`: det nonzero.

## Reference logs added for v13 (22 Oct 2026, HKT)
`audit/corner_check.log`, `audit/swap_uniqueness.log`, `audit/tetra_arb.log` (Gram determinant 213.25944433770 +/- 4e-12 > 0 for V1 and V1' at theta*), `audit/tetra_check.log`.
