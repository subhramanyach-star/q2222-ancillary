# Referee checklist: computer-assisted steps of prl_v14

"Exact Optimization of Bell Witnesses for Device-Independent Randomness", R. Ramanathan and Y. Wu (HKU).
Numbering as in prl_v14. Run each command inside the directory named in its row. Pass = the printed marker shown.
Times: 4-core box. **[fast]** = under about 2 min.

## 0. Integrity
- [ ] `md5sum -c --quiet HASHES.md5` (top level) prints nothing. Certificates/data: every `*.json` and `thm2a_even_polys.txt` is listed.

## 1. Theorem 2 covers all of 0 < theta <= pi/4 (x = theta/pi; q = sin^2 2theta)

| range | vertex | where / command | pass | time |
|---|---|---|---|---|
| polygon combinatorics, all (0, pi/4] | all | `tilted_chsh_fullrange/fullrange`: `python3 comb.py` | `ALL OK (open intervals): True` | [fast] 10 s |
| q in [0, 3/20] (theta <= 0.0634 pi) | V1 (party swap, **Alice-overlap** slices) | `.../smalltheta`: `python3 sgen.py c_V1s.json` (q <= 1/63); `python3 sgen.py c_V1s_ext.json` (1/63..3/20) | `"ok": true` on every strip | [fast] 10 s; ext ~1 min |
| q in [0, 3/20] | V3 (Bob slices, exact /q) | `python3 tgen3.py c_V3t.json` (q <= 1/63; `--tstrips 0 1` = 50 s); `python3 tgen3.py c_V3m.json` (1/63..3/20) | `ALL OK = True` | 3 min; 10 min |
| x in [1/16, 1/8] | V1 / V3 | `.../fullrange`: `python3 vgen.py c_V1l.json` / `c_V3l.json` (mid-pin at c=-1), `--strips j0 j1` | `ALL OK = True` | 7 min / 4 min |
| x in [1/8, 3/16] = Theta | V1 / V3 | `tilted_chsh`: `python3 verify_thm3.py V1 --jobs 4` / `V3 --jobs 4` (`--strips j0 j1` = strips j0..j1-1 of 64) | `ALL OK = True` | 60 s / 10 min |
| x in [3/16, 809/4000] | V1 (pinned near theta*) | `.../fullrange`: `python3 vgen.py c_V1h.json` | `ALL OK = True` | 25 s |
| x in [3/16, 1/4] | V3 | `python3 vgen.py c_V3h.json` | `ALL OK = True` | 50 s |
| x in [632/3125, 1/4] | Wa, Wb (pin at c=+1) | `python3 vgen.py c_Wa.json`, `c_Wb.json` | `ALL OK = True` | 1.5 min / 40 s |
| V2, V4 | closed form, all theta | `.../fullrange`: `python3 v24_conditions.py`, `python3 v24_datacheck.py`; `tilted_chsh`: `python3 v24_hand.py` | all `True` | [fast] |

- [ ] Overlaps: q = 3/20 > sin^2(pi/8) (theta = 0.0634 pi > pi/16); 809/4000 = 0.20225 > theta*/pi = 0.2022406; 632/3125 = 0.20224 < theta*/pi. Each strip claims only the part inside its vertex's regime.
- [ ] Primed vertices V1', V3', Wa', Wb' follow by Sigma: `python3 sigma_check.py` -> `sigma image exact: True`. [fast]
- [ ] Redundant second cover of x in [0.03, 1/16]: `c_V1ll.json`, `c_V3ll_par.json` (optional).

## 2. Endpoints
- [ ] theta = pi/4: octagon = BB octagon vertex by vertex: `.../fullrange`: `python3 pi4_check.py` -> `ALL MATCH O: True`. [fast]
- [ ] theta = pi/4, analytic counterpart Theorem 3(a): `tsirelson_octagon`: `python3 exact_certificates.py`, `exact_certificates_odd.py`, `verify_thm2a_even.py` (32 polynomials), `step5_equality_sets.py` -> `SUMMARY ... True`. [fast]
- [ ] theta -> 0: degenerate corner (c,q) = (-1,0): `.../audit`: `python3 corner_check.py` (P = M = 0 exactly; Taylor coefficients). Uniform-in-q certificates reach q = 0 (rows above). [fast]
- [ ] c = +-1 endpoint lemma (Lemma 5) hypothesis (H2) on Theta: `tilted_chsh_audit`: `PRL_TEX=<path>/prl_v14.tex python3 check_exact.py` part (e) (all 2x8 mixed terms < 1). Off Theta, the strip verifiers check the endpoint terms themselves (or use the pin). [4 min]

## 3. Transition theta*
- [ ] sin^2 2theta* is the real root of x^3 - 8x^2 + 24x - 16, theta* = 0.2022406 pi: `tilted_chsh`: `python3 hexagon_m2.py`; `comb.py` (exactly the transition margins vanish at s*). [fast]
- [ ] Tetrahedral faces of V1, V1' at theta*: `.../audit`: `python3 tetra_arb.py` -> `Gram det = [213.25944433770 +/- 4e-12] >0: True` (paper: 213.259...). [fast]

## 4. Exact identities
- [ ] Forced data F = dF/dc = 0 at c0 for V1, V3, Wa, Wb, exact in Q(s)(C,r): `.../fullrange`: `python3 exact_c0_all.py` -> `ALL c0 IDENTITIES EXACT: True`; on Theta also `tilted_chsh`: `python3 exact_c0_signs.py` (interval side conditions, every line `OK`). [fast]
- [ ] Plane / stationarity identities (g = sum mu gamma, mu.m = 0, g.p = 1): `python3 plane_check.py` -> `ALL PLANE IDENTITIES EXACT: True`. [fast]
- [ ] V2/V4 closed-form certificates and side conditions (incl. A2 > 0): `.../fullrange`: `python3 v24_conditions.py` -> `V2/V4 CONDITIONS OK on (0,1]: True`; shipped data = printed closed forms: `python3 v24_datacheck.py` -> `... (exact): True`. [fast]
- [ ] Printed formulas = shipped data (Table III, V1/V3 formulas, normalisations, saturations): `tilted_chsh_audit`: `PRL_TEX=<path>/prl_v14.tex python3 check_exact.py` -> 303 PASS, 0 FAIL. [4 min]

## 5. Strictness, faces, exposure
- [ ] Strict inequalities: 108 (hexagon, 0 < s < s*) and 144 (octagon, s* < s <= 1), all mu_{a|x} > 0: `comb.py`. [fast]
- [ ] Equality sets in {-1, c0, 1}: built into every B&B acceptance test (F > 0 off these points); for P_T, `step5_equality_sets.py`.
- [ ] Party-swapped V1 face data: `.../audit`: `python3 swap_uniqueness.py` (det(d0,d1) = -(1+s^2)/2 != 0). [fast]
- [ ] Faces/exposure ingredients (direction determinacy, affine independence): `tilted_chsh_audit`: `PRL_TEX=<path>/prl_v14.tex python3 check_misc.py` -> `max Phi on grid (must be <= 1): 0.99999...`. [2 min]
- [ ] **Negative controls must FAIL** (`ALL OK = False` / `"ok": false`): `.../fullrange`: `vgen.py c_Wa_neg.json`, `c_Wa_negpin_[a-d].json` (pin path), `c_V3l_negmid_[01].json` (mid-pin path); `.../smalltheta`: `sgen.py c_V1s_neg[a-e].json`, `tgen3.py c_V3t_negep.json`, `c_V3t_negF*.json`, `../audit/c_V3t_negT4_[01].json`. [fast each]

## 6. Proposition 7 (I_phi family; App. H)
- [ ] `wbc_I_family`: `python3 rerun_I_proofs.py` -> `bnb_ok: True`, 84 741 and 128 791 boxes; `python3 polygon_exact.py`; `python3 check_polygon_endpoints.py`; `python3 check_printed_prop6.py` (= Prop. 7) -> `ALL PRINTED PROP. 6 DATA MATCH: True`. [fast]

## Not a proof
## 7. Results 2-3, Fig. 2, Table IV, Lemma 10 (App. I): white-noise min-entropies
- [ ] `rigorous_rates`: `python3 verify_all.py` -> 299 rows `H in [lo, hi]`, none with `FAIL` (299 cases, 315 exact certificates: 262 level 1+AB, 39 level 2, 14 level 3; largest enclosure 1.7e-6 bits); writes `out/verified.json`; `python3 make_table.py` rebuilds `TABLE.md`. [fast] 2 s
- [ ] Headline numbers (pi/8, v = 0.999): sup over hexagon 0.832443, tilted 0.776634, CHSH 0.226659, full [0.835087, 0.835088]: `TABLE.md`. Ratios 3.6727 / 2.0778 / 1.0627 (vs CHSH), 1.0719 / 1.1425 / 1.1147 (vs tilted), 0.99683 (vs full): `PROOFS.md` "Derived numbers". [fast]
- [ ] von Neumann attack upper bounds for every g in the hexagon: `python3 verify_vn.py` -> writes `out/vn_verified.json`. [fast] <1 s

## 8. Result 4, Fig. 3, Proposition 8, Table V (App. J): detection loss
- [ ] `eta_study`: `python3 verify_eta.py` -> one row per certificate file (56); 54 rows contain `PROVED`: tilted and cone at every v = 1 grid point 0.69-0.995 (31) and every v = 0.99 grid point 0.75-0.995 (17), cone in the six AB11 rows (two-party). Not proved: eta = 0.68 (v = 1); `eta0.755_v0.99.json` is a CHSH-only file. Writes `out/verified_eta.json`. [fast] 2 s
- [ ] Tampering: `python3 neg_tests.py` -> 7 of 8 tampered certificates rejected; the 8th (alpha changed) yields a valid, much weaker bound. [fast]
- [ ] vN numbers (Table 3 of `REPORT.md`): exact-PSD BFF records in `out/vn_v1.0.jsonl` (regeneration `vn_opt.py`, needs cvxpy; slow). Quadrature nodes are double precision (App. J).

## Not a proof
`rates_numerics/` (floating-point NPA scripts that first located the rate certificates); the numerical curves in Fig. 3(a).
