# Notes on this directory

This is a copy of the independent audit of Theorem 2 (tilted CHSH; "Theorem 3" in the draft in which it was written,
21 Oct 2026 HKT). Only input paths differ from the original:
- `parse_tex.py` reads the paper source from `$PRL_TEX`, else `../../prl_v13.tex`, `../../prl_r4.tex` or `../../prl.tex` (v13 bundle: `prl_v13.tex` added first).
- `check_exact.py` part (b) reads `../tilted_chsh/{thm3_data.json, certs_v13.json}`.
- `run_full.sh` changes into its own directory.
The B&B outputs are in `bnb_V1.jsonl`, `bnb_V3.jsonl` (one JSON record per strip); `bnb_V*.err` hold the peak memory.

Logs: `check_exact.log`, `bnb_V*.jsonl`, `run_full.log` come from the original run (check_exact 227 s; B&B V1 162 s and V3
205 s wall, about 65 MB). On 21 Oct 2026 (HKT) `check_exact.py` was rerun against the current source (prl_r4.tex): 303 PASS,
0 FAIL, 229 s, with PASS/FAIL lines identical to `check_exact.log`; `check_hexrange2.py` was rerun with identical output.

Text-only correction (21 Oct 2026 HKT). In `endpoint_lemma.tex`, part (iv) used to conclude Phi_eps(delta) = h_{K_eps}(g) <= 1
under (H1)-(H2). It now concludes Phi_eps(delta) <= 1, and Phi_eps(delta) = h_{K_eps}(g) = 1 when some L in D_eps has g.L = 1.
A counterexample to the old statement is included (g = 0, delta = (1/4, 0)). The corollary names the tight point that V1, V3,
V1', V3' have at each endpoint (from the exact facet saturations of check (a)). No code, data or log changed.

Rerun for the v13 bundle (22 Oct 2026, HKT): `check_exact.py` against prl_v13.tex: 303 PASS, 0 FAIL, 227 s; output in `check_exact_prl_v13.log` (PASS/FAIL lines identical to `check_exact.log`; `check_exact_summary.json` byte-identical).
