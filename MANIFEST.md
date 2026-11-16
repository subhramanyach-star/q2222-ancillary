# Manifest of computer-assisted checks

This file maps the computer-assisted statements of the accompanying paper to the verifier, the certificate files, and the expected output. Run each command from the root of the repository unless a row says otherwise. A check has succeeded only if the process exits with status 0 and prints the stated marker.

The normal-cone construction, the exact quantum-bound formula, the necessary condition for certified-randomness maxima, the reduction of the loss model to two-qubit strategies, and the canonical CHSH maximum under equal-efficiency loss are proved in the paper and have no computer-assisted step here.

## Detection loss

| Statement in the paper | Command | Inputs | Expected output |
|---|---|---|---|
| Bell expressions maximized by the numerically optimized strategies (112 efficiencies, 656 vertices) | `python3 lossycone/vpt_t4.py` | `lossycone/data/t4/<eta>.json` | `RESULT ALL OK = True`; 112 points, 656 vertices accepted |
| The same construction at the six table efficiencies | `python3 lossycone/vpt.py` | `lossycone/data/points.json` and `lossycone/data/certs.json` | `RESULT ALL OK = True` |
| Explicit Bell expressions certify randomness above the CHSH bound (lower ends of the selected-expression intervals, and the CHSH caps) | `python3 lossycone/verify_t4_rates.py` | the same T4 certificates | `ALL OK = True (112 points)` |
| Upper ends of the selected-expression intervals at the six reference efficiencies \(0.75,0.80,0.85,0.90,0.95,0.98\) | `python3 lossycone/verify_t4_attack.py` | T4 attack data; the six published upper endpoints | certified \(H_{\mathrm{hi}}\) at most \(0.009847\), \(0.037864\), \(0.095105\), \(0.198627\), \(0.385566\), \(0.597368\); `RESULT ALL OK = True` |
| The constructions cover every efficiency in \([0.75,0.98]\) | the interval-guarantee block of `python3 lossycone/verify_t4_rates.py` | the T4 grid | adjacent-band guaranteed lower gains \(5.32\%\), \(8.52\%\), \(10.78\%\), \(12.21\%\), \(8.04\%\); grid-point lower end at least \(8.8\%\); table-efficiency lower end at least \(9.2\%\) |
| The remaining CHSH sign class lies strictly below the certified canonical maximum | `python3 chsh_variants/verify_variants.py` | `chsh_variants/class2_certs.json` and the T4 certificates | `ALL OK = True`; minimum certified gap at least \(3.19\times 10^{-3}\) |
| Rounding of interval endpoints to floating-point values | `python3 rigorous_rates/test_certified_float.py` | none | `certified_float unit check: ALL OK` |
| A tampered cone certificate is rejected | `python3 lossycone/tamper.py` | `lossycone/data/points.json` and `lossycone/data/certs.json` (the shipped files are not modified) | `RESULT TAMPER ALL OK = True`; an unperturbed certificate is accepted and every perturbation is rejected |

The six-decimal selected-expression intervals in the paper, whose upper ends are certified by the attack path and whose lower ends are certified by the rate path, are:

| \(\eta\) | interval |
|---|---|
| \(0.75\) | \([0.009845, 0.009847]\) |
| \(0.80\) | \([0.037859, 0.037864]\) |
| \(0.85\) | \([0.095104, 0.095105]\) |
| \(0.90\) | \([0.198626, 0.198627]\) |
| \(0.95\) | \([0.385565, 0.385566]\) |
| \(0.98\) | \([0.597367, 0.597368]\) |

These are intervals for the selected-expression guessing problem, under the loss model stated in the paper.

## Tilted-CHSH supporting polygons

The theorem on Bell tests maximized by tilted-CHSH sources is checked in three parts: exact identities for the vertex formulae (including the closed-form \(V_2,V_4\) bounds and the comparison with the Tsirelson octagon at maximal entanglement); exact endpoint identities; and a ball-arithmetic cover of the overlap interval. Run these commands in the directory named in the first column.

| Directory | Command | Expected output |
|---|---|---|
| `tsirelson_octagon/` | `python3 exact_certificates.py && python3 exact_certificates_odd.py && python3 verify_thm2a_even.py && python3 step5_equality_sets.py` | each script reports success (`True` / `ALL OK = True`) |
| `tilted_chsh_fullrange/fullrange/` | `python3 comb.py && python3 exact_c0_all.py && python3 plane_check.py && python3 sigma_check.py && python3 pi4_check.py && python3 v24_conditions.py && python3 v24_datacheck.py` | `ALL OK (open intervals): True`, `ALL c0 IDENTITIES EXACT: True`, `ALL PLANE IDENTITIES EXACT: True`, `sigma image exact: True`, and the corresponding success markers for the remaining scripts |
| `tilted_chsh_fullrange/fullrange/` | `python3 vgen.py c_Wa.json --strips 23 24` | `ALL OK = True` |
| `tilted_chsh_fullrange/fullrange/` | `python3 vgen.py c_Wa_negpin_a.json` | the negative control reports `ALL OK = False` |
| `tilted_chsh_fullrange/smalltheta/` | `python3 sgen.py c_V1s.json` | `ALL OK = True` |
| `tilted_chsh/` | `python3 verify_thm3.py V1` and `python3 verify_thm3.py V3` | `ALL OK = True` for each family |

A floating-point check is not sufficient for these certificates, because the margin vanishes exactly at the target overlap and at tight endpoints. Those zeros are handled in exact arithmetic; the remaining interval is covered by enclosures in ball arithmetic.
