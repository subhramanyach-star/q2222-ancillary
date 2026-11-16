# MANIFEST: claim -> certificate -> verifier (manuscript prl_v19)

R. Ramanathan and Y. Wu (HKU). Numbering is that of **prl_v19** (Letter: Results 1-5, Eq. (1), Fig. 1; Supplement:
Theorems S1-S6, Lemma S1, Propositions S1-S2, Tables S1-S4). Run every command inside the directory named in its row
(data are loaded by relative path). Software: Python 3.13, python-flint 0.9 (Arb), numpy, sympy; no verifier below
calls an SDP or conic solver unless marked **[solver]**. Integrity: `md5sum -c --quiet HASHES.md5` (prints nothing).

**Repository: placeholder.** Ref. [anc] of prl_v19 still reads "[REPOSITORY URL TO BE ADDED]". *Submission blocker*:
the bundle must be deposited at a fixed, citable location (e.g. Zenodo DOI or a tagged GitHub release) before submission.

## A. Proved by hand (no computer step)

| claim | where proved |
|---|---|
| Theorem S1 (exact dual, Eq. (1)), given the generation theorem of Ref. [RamanathanNDWU] | Supp. Sec. S1, Steps 1-5 |
| Lemma S1 (properties of P_g), Theorem S2 (the best witness is maximised by its own data; necessary condition) | Supp. Sec. S2 |
| Proposition S1 (reduction under loss; pull-back identity M^T g in N_Q(q)) | Supp. Sec. S2 |
| Proposition S2 (i)-(ii) (CHSH variants: relabelling reduction to two classes; local random recordings) | Supp. Sec. S4 |
| EAT "Fact" (full statistics is a single score in EAT) | Supp. Sec. S5 |

## B. Computer-assisted claims (proofs)

| claim (prl_v19) | verifier: directory / command | inputs | expected output | dependencies | time |
|---|---|---|---|---|---|
| **Theorem S3** (Result 3): maximised set of p_theta = C(p_theta), all theta in (0, pi/4]; hexagon/octagon; theta* | see `VERIFY.md` sections 1-5 (prl_v14 numbering: its "Theorem 2" = Theorem S3). Core: `tilted_chsh_fullrange/fullrange`: `python3 comb.py`, `exact_c0_all.py`, `plane_check.py`, `sigma_check.py`, `pi4_check.py`, `v24_conditions.py`, `v24_datacheck.py`, `vgen.py c_*.json`; `tilted_chsh`: `verify_thm3.py V1|V3`; `tilted_chsh_fullrange/smalltheta`: `sgen.py`, `tgen3.py` | certificate JSONs next to each script (`c_*.json`, `certs_v13.json`, `thm3_data.json`) | `ALL OK (open intervals): True`, `ALL c0 IDENTITIES EXACT: True`, `ALL PLANE IDENTITIES EXACT: True`, `sigma image exact: True`, `ALL MATCH O: True`, `ALL OK = True` per strip; negative controls `ALL OK = False` | sympy, python-flint | fast checks ~2 min; full strips up to ~10 min each (listed in `VERIFY.md`) |
| BB octagon as endpoint theta = pi/4 (Result 3) | `tilted_chsh_fullrange/fullrange`: `python3 pi4_check.py`; `tsirelson_octagon`: `exact_certificates.py`, `exact_certificates_odd.py`, `verify_thm2a_even.py`, `step5_equality_sets.py` | in-folder JSON/txt | `ALL MATCH O: True`; `SUMMARY ... True` | sympy, flint | ~5 s |
| **Theorem S4** (Result 4): exact maximised sets M^{-T} C(q_eta) at the best devices found, 6 main + all 112 grid efficiencies (656 vertices) | `lossycone`: `python3 vpt_t4.py` | `data/t4/<eta>.json` (rational device, polygon, per-vertex certificate curves delta(c)) | last line `RESULT ALL OK = True`; reference `vpt_t4.log` | flint (160-bit Arb), `exact.py`, `vpt.py` | ~3 min |
| Theorem S4, negative controls | `lossycone`: `python3 tamper.py` | perturbed certificates of the six main devices (`data/points.json`, `data/certs.json`) | every line `REJECTED` | flint | <1 min |
| **Theorem S5** (Result 4, Table S2 "our witness" and "CHSH bench. <=" columns, Fig. 1(a)) | `lossycone`: `python3 verify_t4_rates.py` | `data/t4/<eta>.json` (witness weights, exact NPA dual certificates per guess, CHSH-cap certificate) | per point `H_lo >= ...  every CHSH <= ...  OK`; e.g. `0.900 ... 0.198626 ... 0.168313`; `ALL OK = True (112 points)`; reference `verify_t4_rates.log` | `../rigorous_rates/verify_all.py` (exact LDL^T), `../eta_study/cert_eta.py` (loss functional), flint | ~10 s |
| **Theorem S6** (Result 4: every eta in [0.75, 0.98], interval gains 5.33/8.52/10.79/12.21/8.04 %) | same run of `verify_t4_rates.py` (block "interval guarantees") | as above | `[0.75, 0.80]: gain >= 5.33%` ... `[0.95, 0.98]: gain >= 8.04%`; `grid-point gains: 8.87% .. 18.01%` | as above + Proposition S2(iii) below | (same run) |
| **Proposition S2(iii)** (CHSH benchmark: second class of sign variants below S_up at all 112 grid points) | `chsh_variants`: `python3 verify_variants.py` | `class2_certs.json` (exact level-3 certificates), `../lossycone/data/t4/<eta>.json` | 112 lines `OK`; `min margin 3.200e-03 at eta = 0.980`; `ALL OK = True`; reference `verify_variants.log` | `rigorous_rates/verify_all.py`, `eta_study/cert_eta.py` | ~2 s |
| Proposition S2(iii), certificate generation | `chsh_variants`: `python3 class2_cert.py` **[solver]** (overwrites `class2_certs.json`) | none | `ALL OK = True 112`; reference `class2_cert.log` | cvxpy/Clarabel, `rigorous_rates/rates_core.py` | ~15 s |
| Table S2, "asym. CHSH, best found" and "polygon, best found" columns (certified lower bounds at the stored devices) | `verify_extra`: `python3 verify_lower.py` (wapcone part) | `../wapcone/certs/eta*_{wap,comb,full}.json` | one `OK` line per file | `rigorous_rates/verify_all.py`, `eta_study/cert_eta.py` (polygon vertices in Arb) | (same run) |
| **Proposition S2(iv)** (CHSH benchmark under source noise: white / dephasing / experiment-matched rows) | `chsh_variants`: `python3 verify_noise_variants.py` | `noise_class2_certs.json` | `ALL OK = True`, min margin `1.960e-02`; reference `noise_class2.log` | `rigorous_rates/verify_all.py` | ~2 s |
| Proposition S2(iv), certificate generation | `chsh_variants`: `python3 noise_class2.py` **[solver]** | `../noise/certs/wd_*.json` for the eight experiment rows | `ALL OK = True n= 23`; writes `noise_class2_certs.json` | cvxpy/Clarabel, `rates_core` | ~5 s |
| **Table S4 / Result 5** (experiment-matched: CHSH cap in the noise model, tilted and polygon lower bounds) | `verify_extra`: `python3 verify_lower.py` (noise part; rows `wd_eta*`) | `../noise/certs/*.json` | e.g. `wd_eta0.8749_q0.9936-1.0.json S_up exact True every CHSH <= 0.111332 tilted >= 0.116178 cone >= 0.121833 OK`; final `ALL OK = True (119 files)`; reference `verify_lower.log` | as above | ~5 s |
| Earlier-witness and white-noise studies cited in the companion paper | `eta_study`: `verify_eta.py`, `neg_tests.py`; `rigorous_rates`: `verify_all.py`, `verify_vn.py` | in-folder certs | see `VERIFY.md` sections 7-8 | flint | fast |

Upper ends of the two-sided intervals (explicit qubit attacks, ball arithmetic) are stored in the `attack` blocks and
verified in ball arithmetic by the producing scripts (`rigorous_rates/rates_core.verify_attack`); the stand-alone verifiers above re-check only the lower bounds and the CHSH caps, which are all that the separations need.

## C. Numerical (not proved; no claim of the manuscript is a theorem about these)

| item | script | output |
|---|---|---|
| Table S1 (benchmark of Eq. (1) against closed forms; NPA columns) | `socp`: `python3 bench.py` **[solver]**, `mktable.py` | `bench.json`, `bench.log`, `bench_table.tex` |
| Eq. (1) evaluator | `socp/socp_core.py` **[solver]** | upper bounds by weak duality for any returned delta; exactness needs the interval check of Theorems S3-S4 |
| Device optimisation (best devices found) | `lossycone`, `wapcone/cert_wc.py` (design script; imports a search module not shipped), `noise/nz.py` | stored rational devices only; optimality is **not** claimed |
| "single score reproduces full statistics to 1e-6 bits at 5 of 6 efficiencies" | floating-point NPA level 3 | `lossycone/REPORT.md` |
| Table S3 and Sec. S5 (finite size, EAT, i.i.d. confidence regions) | `finitesize`: `run_eat.py`, `cr.py`, `cr_opt_chsh.py`, `check_asym.py` **[solver]** | `out/`, `logs/` |
| Table S4 "best CHSH" and "est. gain" columns | `noise/nz.py` **[solver]** | `noise/out/tables.md` |
| von Neumann rates (companion paper only) | `eta_study/vn_opt.py`, `rates_v11/` | exact-PSD dual, but Gauss-Radau nodes/weights in double precision: **numerical**, not certified |

## D. Known gaps
1. Repository URL (blocker, above).
2. `wapcone/cert_wc.py` and `noise/ncert.py` are generation scripts; the stand-alone check is `verify_extra/verify_lower.py`.
