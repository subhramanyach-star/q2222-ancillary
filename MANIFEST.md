# Manifest of checks

This file maps the statements of the accompanying paper that are supported by a shipped script to the command, its inputs, and the expected output. Run each command from the root of the repository unless a row says otherwise. A check has succeeded only if the process exits with status 0 and prints the stated marker.

The normal-cone construction, the exact quantum-bound formula, the supporting quadrilateral of the I_δ family, the necessary condition for certified-randomness maxima, the reduction of the loss model to two-qubit strategies, and the canonical CHSH maximum under equal-efficiency loss are proved in the paper. The I_δ script checks the exact identities of that proof; it is not a computer-assisted existence argument. The detection-loss and CHSH-variant statements below are computer-assisted and are checked by the remaining verifiers.

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

## I_δ supporting quadrilateral

The vertices, multipliers, and local margins of the supporting quadrilateral at every behaviour of the I_δ family are derived by exact algebra in the paper. The following command re-derives those closed forms and checks the identities used in the appendix. It is a check of exact algebra, not a computer-assisted existence proof.

| Statement in the paper | Command | Inputs | Expected output |
|---|---|---|---|
| Closed-form vertices, multipliers, and local slacks of the I_δ quadrilateral, for every 0 < δ ≤ π/6 | `python3 idelta/certify_idelta.py` | none | `All I_δ identities and interval sign checks passed.` |
