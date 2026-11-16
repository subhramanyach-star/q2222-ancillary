# Ancillary files for prl_v19 (R. Ramanathan and Y. Wu, HKU), 27 Oct 2026 (HKT)

**Start with `MANIFEST.md`**: it maps every claim of prl_v18 (Theorems S1-S6, Propositions S1-S2, Tables S1-S4, Results 1-5)
to its verifier, inputs, expected output and dependencies, and lists what is numerical only.
This bundle is the prl_v17 bundle with: (1) new `chsh_variants/` (Proposition S2(iii)--(iv): other CHSH sign variants and
recordings of missed detections are no better than the canonical CHSH, certified at all 112 grid efficiencies);
(2) source-noise class-2 certificates (`noise_class2_*`, `verify_noise_variants.py`); (3) `verify_extra/verify_lower.py`, a stand-alone verifier for the certified lower bounds and CHSH caps of
`wapcone/` (Table S2 comparison columns) and `noise/` (Table S4); (3) absolute paths in `noise/`, `wapcone/` and
`mkfig_v16.py` replaced by relative ones. Repository URL: still a placeholder (submission blocker).
The sections below document earlier versions with their own numbering.

---

# Ancillary files for "Exact Optimization of Bell Witnesses for Device-Independent Randomness" (manuscript prl_v16)

R. Ramanathan and Y. Wu (The University of Hong Kong). Bundle for manuscript version **prl_v16** (23 Oct 2026, HKT).
It contains the prl_v14 bundle unchanged (documented below, with prl_v14 numbering) plus four new folders:
`lossycone/`, `noise/`, `wapcone/`, `fritz/`, and `mkfig_v16.py` (draws Fig. 1 of prl_v16 from `lossycone/t4_summary.json`).

## New in prl_v16 (prl_v16 numbering)

| directory | supports | kind |
|---|---|---|
| `lossycone/` | **Results 1 and 3, Fig. 1, App. B**: Theorem 6 (exact normal cones at the lossy optima), Theorem 7 (certified rates), Theorem 8 (every eta in [0.75, 0.98]), Table I. Device data `data/points.json` (six main devices of the first pass) and `data/t4/<eta>.json` (the 112 certified efficiencies used in the paper: rational device, exact polygon, per-vertex certificates delta(c), witness weights, NPA certificates, CHSH-cap certificate, attack). Verifiers: `vpt_t4.py` (cone side at all 112 points; reference output `vpt_t4.log`), `verify_t4_rates.py` (rate side, CHSH cap and the interval guarantees; reference output `verify_t4_rates.log`), `vpt.py` (cone side for `data/points.json`), `tamper.py` (negative controls). Theorem 5 and Propositions 4-5 are pen-and-paper; `theorem.tex/pdf` is the working note. | exact + Arb |
| `wapcone/` | asymmetric-CHSH and C_theta' columns of Table I (certificates `certs/`, `cert.jsonl`; `REPORT.md`) | exact + Arb (certificates produced by `cert_wc.py`; no separate stand-alone verifier) |
| `noise/` | **Result 4, App. M, Table VII** (source noise, experiment-matched parameters): certificates `certs/`, tables `out/tables.md`, scripts `nz.py` (optimiser), `ncert.py` (certification), `REPORT.md`, `STATUS.md` (experiment parameters and sources) | exact + Arb (certified witness side); numerical best CHSH protocol under dephasing |
| `fritz/` | **App. N, Proposition 11 and the explicit example**: `cert_q14.json`, exact-rational verifier `verify_q14.py` (Python `fractions` only), builder `make_cert.py`, notes `LINK.md` | exact |

Quick start for the new parts (about 6 minutes, 1 core; needs python-flint and numpy):
```
cd lossycone && python3 vpt_t4.py > my_vpt_t4.log && tail -1 my_vpt_t4.log      # RESULT ALL OK = True (112 points, 656 vertices)
python3 verify_t4_rates.py | tail -8                                              # gains, interval guarantees, ALL OK = True
python3 tamper.py && cd ..                                                         # perturbed certificates rejected
cd fritz && python3 verify_q14.py && cd ..                                         # VERIFIED
```

### Numbering map prl_v14 -> prl_v16 (for the folders documented below)
Theorems 1-4 unchanged; Theorem 5 -> 9. Lemmas 1-10 -> 3-12 (so Lemma 5 endpoint -> 7, Lemma 6 scope -> 8, Lemma 10 two-sided bounds -> 12).
Propositions 4, 5, 6, 7, 8 -> 6, 7, 8, 9, 10 (Prop. 7 WBC I-family -> 9; Prop. 8 detection loss -> 10). Corollary 4 -> 5.
Tables I-V -> II-VI. Figures 2, 3, 4, 5 -> 4, 5, 2, 3. Appendices B-H -> C-I, I -> K, J -> L, K -> O.
Results 2-3 of prl_v14 (white noise) are now App. J; Result 4 of prl_v14 (earlier witnesses under loss) is now App. L.

---

# Ancillary files for "Exact Optimization of Bell Witnesses for Device-Independent Randomness"

R. Ramanathan and Y. Wu (The University of Hong Kong). Bundle for manuscript version prl_v14 (22 Oct 2026, HKT). It contains the prl_v13 bundle unchanged (apart from the
notes below) plus `rigorous_rates/`, `eta_study/` and `rates_v11/`.

These files support the computer-assisted steps of the paper and the numerical rate illustration. All theorem,
proposition, lemma, corollary, appendix, table and figure numbers below are those **of prl_v14** (taken from its
LaTeX labels). A one-page referee checklist is in `VERIFY.md`; md5 hashes of every shipped file except `HASHES.md5` and `SMOKETEST.log` are in `HASHES.md5`
(`md5sum -c HASHES.md5`); the smoke test of this bundle from a fresh unpack is in `SMOKETEST.log`.

## Quick start (about 5 minutes in total, 4-core machine)

```
pip install sympy python-flint mpmath numpy scipy cvxpy matplotlib
md5sum -c --quiet HASHES.md5                          # silent = all files intact
cd tsirelson_octagon && python3 exact_certificates.py && python3 exact_certificates_odd.py \
   && python3 verify_thm2a_even.py && python3 step5_equality_sets.py && cd ..      # Thm 3(a), ~4 s
cd tilted_chsh_fullrange/fullrange && python3 comb.py && python3 exact_c0_all.py && python3 plane_check.py \
   && python3 sigma_check.py && python3 pi4_check.py && python3 v24_conditions.py && python3 v24_datacheck.py   # Thm 2 exact parts, ~1 min
python3 vgen.py c_Wa.json --strips 23 24 && python3 vgen.py c_Wa_negpin_a.json && cd ../..   # Arb strip + negative control
cd tilted_chsh_fullrange/smalltheta && python3 sgen.py c_V1s.json && cd ../..                # small-theta V1 (Alice overlap)
cd tilted_chsh && python3 verify_thm3.py V1 --jobs 4 && python3 verify_thm3.py V3 --strips 63 64 && cd ..   # Theta proof
cd wbc_I_family && python3 rerun_I_proofs.py && python3 check_printed_prop6.py && cd ..   # Prop 7, ~40 s
```

Expected: every exact check prints `True` / `ALL OK = True` / `ALL MATCH O: True`; every negative control prints
`ALL OK = False`. Reference output for every script sits next to it (`*.log`, `*.jsonl`). Run each command from inside its
own directory (data are loaded by relative path). The full runs (V3 Theta proof about 10 min; full-range strips
from seconds to about 7 min each; small-theta V3 about 3 and 10 min) are listed in `VERIFY.md`.

## Map: directory -> statement of prl_v14

| directory | supports (prl_v14 numbering) | kind |
|---|---|---|
| `tsirelson_octagon/` | **Theorem 3(a)** (Tsirelson octagon, App. C): exact Q(sqrt2) root isolation for the 32 even-vertex conditions, odd-vertex certificates, equality sets E_k; Table I saturation pattern | exact |
| `tilted_chsh/` | **Theorem 2** on Theta = [pi/8, 3pi/16] (App. D, "Proof of validity"): certificates of V1, V3 (Table III = `certs_v13.json`), stand-alone Arb verifier, search B&B, **Lemma 5** (endpoint lemma) hypotheses (H2) via exact c0 data; closed-form V2/V4 certificates (`v24_hand.py`) | Arb + exact |
| `tilted_chsh_audit/` | independent re-implementation of the Theta proof of Theorem 2 **from the formulas printed in the paper source** (`prl_v13.tex`), incl. faces/exposure ingredients (App. D "Faces and exposure") and Lemma 5 proof (`endpoint_lemma.tex/pdf`) | exact + Arb |
| `tilted_chsh_fullrange/fullrange/` | **Theorem 2 on all of (0, pi/4]** (App. D): exact polygon combinatorics and transition theta* ("Exact combinatorics", 108 hexagon + 144 octagon strict inequalities; **Lemma 4** vertex exhaustion), exact forced data at c0, plane/stationarity identities, Sigma images, pi/4 endpoint = BB octagon ("Endpoints"), V2/V4 side conditions, theta-strip certificates for V1, V3, Wa, Wb ("Certificates on the rest of (0, pi/4]") with negative controls | exact + Arb |
| `tilted_chsh_fullrange/smalltheta/` | **Theorem 2, small theta** (App. D "Small theta", q = sin^2 2theta in [0, 3/20]): V1 by party swap / Alice-overlap slices (`sgen.py`), V3 by Bob-overlap slices with exact division by q (`tgen3.py`), negative controls | exact + Arb |
| `tilted_chsh_fullrange/audit/` | corner (c,q) = (-1,0) Taylor data, uniqueness of swapped slice data (det = -(1+s^2)/2, App. D "Vertices"), **tetrahedral faces at theta*** (Gram det 213.259..., App. D "Dimensions"), spot re-runs | exact + Arb |
| `wbc_I_family/` | **Proposition 7** (certified normal cone along I_phi, App. H; Ref. [ancWBC]): W1 certificate B&B on upsilon in [1/40, 3/20] and [3/20, 5359/20000], exact polygon data, printed-formula and Gram-determinant checks | Arb + exact |
| `wbc_I_family/` (cont.) | App. H closed forms for I_phi and J_gamma used with **Lemma 9** and **Corollary 1(b)** (Lambda, g_corr, active cones, m_{+|0}m_{+|1} < 0): derived symbolically in `sym_family.py`, cross-checked numerically in `check_wbc.py` ("re-checked numerically at several parameter values") | symbolic + float |
| `rigorous_rates/` | **Results 2-3**, **Fig. 2**, **Table IV**, **Lemma 10** and App. I (certification of every white-noise min-entropy number: exact-rational NPA certificates at levels 1+AB/2/3 + explicit qubit attacks, supremum over the hexagon, fixed-source comparison, P_T); von Neumann attack upper bounds for the hexagon (App. I "Von Neumann entropy"). Stand-alone verifiers `verify_all.py`, `verify_vn.py`; proofs and tables `PROOFS.md`, `TABLE.md`; `mkfig2_v14.py` draws Fig. 2 from `out/verified.json` | exact + Arb |
| `eta_study/` | **Result 4**, **Fig. 3**, **Proposition 8**, **Table V** and App. J (detection loss): 56 certificate files in `certs/`, stand-alone verifier `verify_eta.py`, tampering tests `neg_tests.py`, CHSH-variant check `variants_check.py`, exact-PSD BFF (`bffexact.py`, `vn_opt.py`) and EAT numbers (`out/vn_v1.0.jsonl`), report `REPORT.md`, tables `out/tables.md`; `mkfig_eta_v14.py` draws Fig. 3 | exact + Arb (H_min); exact-PSD with double-precision quadrature nodes (vN) |
| `rates_v11/` | BFF hierarchy with exact-rational dual certification and the EAT routine (Liu et al. 2021 form), imported by `eta_study/` and the source of the BFF lower bounds for tilted CHSH quoted in App. I | exact-PSD + float |
| `rates_numerics/` | floating-point NPA scripts with which the certificates of `rigorous_rates/` were first located (the prl_v13 Fig. 2 / Table IV data); no statement of prl_v14 depends on them | float |

Items of prl_v14 with no computer-assisted step (proved analytically in the text): Theorem 1, Theorem 4, Theorem 5,
Propositions 1-6, Corollaries 1-4, Lemmas 1-3 and 6-9 (apart from the App. H closed forms noted above). Corollary 3 (noise window) is analytic (App. K); its thresholds
are only drawn in Fig. 2.

## Numbering changes prl_v13 -> prl_v14
New: Result 4, Fig. 3 (detection loss), Lemma 10, Proposition 8, Table V, App. I (now "Certification of the randomness
numbers"), App. J (detection loss). Shifted: the appendix figures (old Fig. 3, 4 -> Fig. 4, 5) and the proof of
Corollary 3 (old App. J -> App. K). All other theorem, lemma, proposition, corollary and table numbers are unchanged, so
the notes and logs of the other folders, written against prl_v13, keep their meaning. App. D box count corrected in
prl_v14: V1 on the rest of (0, pi/4] accepts 790 032 boxes (`fullrange/v_V1h.log` 25 192 + `v_V1l.log` 764 840).

## Quick start for the rate certificates (under 10 s; needs python-flint and numpy only)
```
cd rigorous_rates && python3 verify_all.py && python3 verify_vn.py && cd ..   # 299 cases, ALL OK; vN attacks
cd eta_study && python3 verify_eta.py && python3 neg_tests.py && cd ..         # 56 files; PROVED rows; 7 of 8 tamperings rejected, 1 valid-but-weaker
```

## Legacy file names and notation

Files keep the names under which they were produced and audited (renaming would change hashes and logs):
* `thm2a` (`verify_thm2a_even.py`, `thm2a_even_polys.txt`) = **Theorem 3(a)**.
* `thm3` (`verify_thm3.py`, `thm3_data.json`, `thm3_data.py`) = **Theorem 2**; `certs_v13.json` / `bnb_v13.py`: "v13" is an
  internal certificate version, unrelated to manuscript version prl_v13.
* `check_printed_prop6.py` and "Proposition 6" in `wbc_I_family/` = **Proposition 7** (props were renumbered after an earlier
  Proposition 1 was deleted).
* Inside older notes (`tilted_chsh_audit/AUDIT.md`, `fullrange/REPORT.md`, `tilted_chsh_fullrange/audit/AUDIT.md`): "App. C" for the
  tilted-CHSH appendix = **App. D**; "Theorem 3" for the tilted family = **Theorem 2**; "prl_r4" = the draft they were written against.
  The mathematics is unchanged; only the numbering moved.
* `t` in `tilted_chsh_fullrange/smalltheta/` = the paper's **q = sin^2 2theta**; `t` in `wbc_I_family/` = the paper's **upsilon = tan(phi/2)**.
* Box counts: the paper's 508 798 (V1) and 4 853 223 (V3) on Theta are from the search B&B `bnb_v13.py` (`bnb_v13.log`); the
  stand-alone verifier `verify_thm3.py` uses its own covering and accepts 508 800 and 4 853 214 boxes (`verify_V*_full.log`).

## Software and precision
Python 3.13.5, sympy 1.14.0, python-flint 0.9.0 (Arb), mpmath 1.3.0, numpy 2.4.6; scipy 1.18.1 and cvxpy 1.9.3 are imported
only by design scripts, floating-point cross-checks and `rates_numerics/` (never in an acceptance test).
Arb working precision: 110 bits (`tilted_chsh/`, `wbc_I_family/`), 120 bits (`tilted_chsh_fullrange/`, as in App. D),
128 bits (`tilted_chsh_audit/`). No pickles are shipped; exact data are JSON (sympy expressions as strings; loaders
`thm3_data.py`, `wbc_I_family/polygon_io.py`). Some regeneration scripts write a local pickle cache
(`symbolic_m2.py` -> `m2_symbolic.pkl`, read by `export_json.py`; `hexagon_m2.py` -> `theta_star.pkl`).

## Commands per directory (reference times, 4-core box)
`tsirelson_octagon/` (each about 1 s): `exact_certificates.py` (even vertices), `exact_certificates_odd.py`,
`verify_thm2a_even.py` (stand-alone even-vertex check; writes `thm2a_even_polys.txt`, the 32 polynomials), `step5_equality_sets.py`.

`tilted_chsh/`: `verify_thm3.py V1 --jobs 4` (about 60 s), `verify_thm3.py V3 --jobs 4` (about 10 min); stand-alone, reads only
`thm3_data.json` and `certs_v13.json`; `--strips j0 j1` runs strips j0..j1-1 of 64. `bnb_v13.py V1|V3` (search B&B),
`hexagon_m2.py` (1 s; theta* polynomial), `exact_c0_signs.py` (1 s), `v24_hand.py` (2 s; imports `../wbc_I_family/core.py`),
`symbolic_m2.py` then `export_json.py` (regenerate `thm3_data.json`). `hexagon_extend.py` (2.5 min) is the earlier
interval-arithmetic hexagon check on [pi/16, 0.2022pi]; App. D states it is **superseded** by the exact `comb.py`; its log also
records where it fails near theta = 0, which is why the exact method replaced it.

`tilted_chsh_audit/` (needs the paper source: `$PRL_TEX`, default `../../prl_v14.tex`, i.e. the arXiv layout
`<source>/anc/tilted_chsh_audit/`): `check_exact.py` (about 4 min; 303 PASS, 0 FAIL against prl_v13.tex, `check_exact_prl_v13.log`), `./run_full.sh`
(audit B&B, about 7 min, 65 MB), `summarize_bnb.py V1 V3`, `check_hexrange2.py` (12 s), `check_misc.py` (2 min), `spot_v24.py`. Without the paper source these scripts stop with `FileNotFoundError: 'prl.tex'`;
set `PRL_TEX=/path/to/prl_v14.tex` (App. D of prl_v14 is that of prl_v13 except the corrected box count).

`tilted_chsh_fullrange/`: see its `README.md` (certificate table, strips, box counts, negative controls).

`wbc_I_family/`: `rerun_I_proofs.py` (about 40 s; 84 741 and 128 791 boxes, endpoint lemma at c = +-1), `polygon_exact.py`
(20 s; writes `out/polygon_I.json`), `check_polygon_endpoints.py`, `check_printed_prop6.py` (printed W1, W2, Lambda, g_corr, c0 and the
Gram determinant, exact). `check_wbc.py`, `check_I_touch.py`: floating-point / 60-digit cross-checks. Design only: `sym_family.py`,
`cert_design.py`, `prove_family_auto.py`, `common.py`, `core.py`, `uniform_cert.py`.

`rates_numerics/`: see its `README.md`.

## What was removed relative to the previous bundle (anc_submission_v7)
* `numerics/` (floating-point stability tables and a search history for statements that are not results of prl_v13).
* From `rates_numerics/`: data and scripts not used in prl_v13 (`rates_Id1`, `dualfun`, `mkfig3.py`); `mkfig2_v13.py` added.
* `__pycache__/`. No certificate, data file or acceptance script was changed. Code edits: the default paper path in
  `tilted_chsh_audit/parse_tex.py` (now tries `prl_v13.tex` first); new `tilted_chsh_fullrange/fullrange/v24_datacheck.py`
  (an exact check whose output lines were in the old `v24_conditions.log` without shipped code; that log was regenerated);
  new reference logs for `tilted_chsh_fullrange/audit/{corner_check,swap_uniqueness,tetra_arb,tetra_check}.py` and
  `tilted_chsh_audit/check_exact_prl_v13.log`.

## Added in the prl_v14 bundle
* `rigorous_rates/`, `eta_study/`, `rates_v11/` (development folders copied without `__pycache__/`, scratch files
  `vdbg.py`, `vtest.py`, `t1-4.py`, and without the large solver logs of `rates_v11/`). Code edits on copying: absolute
  paths `/workspace/ndwu/...` in `eta_study/*.py` and `rigorous_rates/make_cases.py` replaced by paths relative to the
  script (`../rigorous_rates`, `../rates_v11`, `../rates_numerics`); `tilted_chsh_audit/parse_tex.py` tries
  `prl_v14.tex` first. No certificate or data file was changed.
* `rigorous_rates/PROOFS.md` section 3, Proposition 6, ends with a remark about equality with the chord for the von
  Neumann entropy; it is not used in prl_v14, which states only the proved bound (at most 0.008 bits at pi/8).

## Added in the prl_v17 / long_v1 bundle (24 Oct 2026)
* `socp/`: the benchmark of the exact quantum maximum (Eq. (1) of prl_v17; Sec. "Exact quantum maxima as a computational
  tool" of long_v1). `socp_core.py` (fixed-overlap second-order-cone dual and the outer maximisation over c; two-qubit
  lower bound), `npa.py` (floating-point NPA levels 1+AB, 2, 3 with tightened Clarabel tolerances), `bench.py` (runs the
  benchmark, writes `bench.json`, log `bench.log`), `mktable.py` (writes `bench_table.tex`). Quick check:
  `python3 -c "from socp_core import hQ_socp; import numpy as np; print(hQ_socp(np.array([0,0,0,0,1,1,1,-1.]))['value'])"`
  prints 2.828427124746193 (CHSH, 2*sqrt(2)). All values are floating point.
* `finitesize/`: finite-size comparison of the CHSH score, the optimal single witness and the full statistics under the
  same entropy accumulation (Liu et al., Nat. Phys. 17, 448 (2021), App. A, Thm. 3, with spot checking), and the i.i.d.
  confidence-region illustration. `fs_core.py` (shared pieces; imports `eta_study/core.py`, `rigorous_rates/npa_float.py`,
  `rates_v11/eat.py` and the device data in `lossycone/data/t4/` via paths relative to the bundle), `check_asym.py`
  (asymptotic values, writes `asym.json`, about 20 s), `run_eat.py <eta> <a|b|c>` (EAT rates, writes `out/eat_*.jsonl`,
  hours), `cr.py`, `cr_opt_chsh.py` (confidence regions, `out/cr*.json`), `REPORT.md` (summary; the development paths
  in it refer to the original working directory). All numbers are floating point and are not certified.
  Code edits on copying: absolute `/workspace/ndwu/...` paths replaced by paths relative to the script.
