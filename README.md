# Ancillary files

This repository accompanies the paper

> Yutian Wu, Xie Sicheng, and Ravishankar Ramanathan, *Explicit construction of the normal cone at every nonlocal extreme behaviour of \(\mathcal{Q}_{2222}\) and device-independent applications*.

It contains the stand-alone verifiers and the certificate data needed to check the computer-assisted claims. The analytic proofs are in the paper and do not depend on these files.

`MANIFEST.md` maps each computer-assisted statement to a command, its inputs, and the expected output. `DEPENDENCIES.md` lists the required software. `SHA256SUMS` records the checksum of every shipped file.

## Requirements

The detection-loss and CHSH-variant checks use Python 3 and [python-flint](https://pypi.org/project/python-flint/) 0.9. The tilted-CHSH polygon checks also use sympy; interval coverage uses python-flint for ball arithmetic. No SDP or conic solver is required to rerun the published checks, and no network access is required once the repository is present.

```
pip install python-flint==0.9.0 sympy
```

Run every command from the root of this repository, or from the directory named in `MANIFEST.md`. Data are loaded by paths relative to each script.

## Quick start

```
python3 lossycone/vpt_t4.py
python3 lossycone/vpt.py
python3 lossycone/verify_t4_rates.py
python3 lossycone/verify_t4_attack.py
python3 chsh_variants/verify_variants.py
python3 rigorous_rates/test_certified_float.py
python3 lossycone/tamper.py
```

Each command should exit with status 0. The first five print a final `RESULT ALL OK = True` or `ALL OK = True` line; the rounding test prints `certified_float unit check: ALL OK`; the last command prints `RESULT TAMPER ALL OK = True` after accepting an unperturbed cone certificate and rejecting every perturbation.

The cone check accepts 112 efficiencies and 656 vertices. The rate and attack checks certify the selected-expression intervals in the paper and the guaranteed gains on every subinterval of \([0.75,0.98]\). The variant check confirms that the remaining CHSH sign class lies a certified gap of at least \(3.19\times 10^{-3}\) below the canonical CHSH maximum at every grid point.

## Tilted-CHSH supporting polygons

```
cd tsirelson_octagon
python3 exact_certificates.py && python3 exact_certificates_odd.py \
  && python3 verify_thm2a_even.py && python3 step5_equality_sets.py

cd ../tilted_chsh_fullrange/fullrange
python3 comb.py && python3 exact_c0_all.py && python3 plane_check.py \
  && python3 sigma_check.py && python3 pi4_check.py \
  && python3 v24_conditions.py && python3 v24_datacheck.py
```

Further strip and small-angle checks are listed in `MANIFEST.md`. A vertex certificate is accepted only when the exact identities and a ball-arithmetic cover of the full overlap interval both succeed; a certificate that fails these checks is rejected.

## Layout

| directory | role |
|---|---|
| `lossycone/` | Supporting-cone and randomness-rate certificates under detection loss, at 112 efficiencies |
| `lossycone/data/points.json`, `lossycone/data/certs.json` | The six table-efficiency devices and their cone certificates |
| `lossycone/data/t4/` | One JSON certificate file for each of the 112 efficiencies |
| `chsh_variants/` | Certificates for the remaining CHSH sign class |
| `tilted_chsh/`, `tilted_chsh_fullrange/`, `tsirelson_octagon/` | Supporting-polygon certificates for the tilted-CHSH family, including the Tsirelson octagon at maximal entanglement |
| `rigorous_rates/` | Library modules used by the rate and variant verifiers (interval arithmetic and exact identities) |

The modules in `rigorous_rates/` are included because the published rate and variant checkers import them. They are not a separate collection of results.
