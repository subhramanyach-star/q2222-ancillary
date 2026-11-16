**prl_v14 note.** Every number that prl_v13 took from these floating-point runs is certified in prl_v14 by `../rigorous_rates/` (exact certificates + attacks); Fig. 2 of prl_v14 is drawn by `../rigorous_rates/mkfig2_v14.py`. The files below are kept as the search that located the certificates.

# rates_numerics/: NPA numerics for Fig. 2 and Table IV of prl_v13 (floating point; not used in any proof)

These are numerical NPA level-3 lower-bound estimates of the min-entropy (App. I of prl_v13), solved with CVXPY + Clarabel
(fallback COPT). They illustrate the theorems; no theorem, proposition or corollary of the paper depends on them.
Needs numpy, cvxpy (with Clarabel), matplotlib. `npa.py`: guessing-probability programs; `families.py`: targets and Bell expressions.

| command | writes | used in prl_v13 | runtime |
|---|---|---|---|
| `python3 rates.py th1` | `rates_th1.json` (log `log_th1.txt`) | Table IV, theta = pi/8; Fig. 2(a) diamonds ("best found in hexagon"); headline 0.227 / 0.777 / 0.832 / 0.835 bits at v = 0.999 | SDP runs, not re-timed |
| `python3 rates.py th3` | `rates_th3.json` (`log_th3.txt`) | Table IV, theta = 27.76 deg = 0.1542 pi (tilt of Liu et al.) | SDP runs, not re-timed |
| `python3 rates.py th2` | `rates_th2.json` (`log_th2.txt`) | App. I "full results, including theta = 3pi/16" | SDP runs, not re-timed |
| `python3 rates.py PT` | `rates_PT.json` | App. I, P_T numbers (CHSH 0.912 ... 0.307; full 0.912 ... 0.336) | SDP runs, not re-timed |
| `python3 rates.py Id2` | `rates_Id2.json` (`log_Id2.txt`) | Fig. 2(b) squares (half-way to W1), WBC target phi = pi/6 | SDP runs, not re-timed |
| `python3 dense.py` | `dense.json` (`log_dense.txt`) | Fig. 2 solid/dashed curves (14 visibilities) | SDP runs, not re-timed |
| `python3 mkfig2_v13.py` | `fig2_rates_v13.pdf` | Fig. 2 (reads `dense.json`, `rates_th1.json`, `rates_Id2.json`) | 2 s |

Exact thresholds drawn as dotted lines in Fig. 2 (Corollary 3, proved analytically in App. J): v = (C+r)/2 = 0.96593 at
theta = pi/8, and v = 5/(3 sqrt3) = 0.96225 for phi = pi/6.

Removed relative to the earlier bundle (not used in prl_v13): `rates_Id1.json`/`log_Id1.txt` (phi = pi/12),
`dualfun.py`/`dualfun.json`/`log_dual.txt` (full-statistics dual functional), `mkfig3.py` (Fig. 3 of an earlier draft;
superseded by `mkfig2_v13.py`, which plots the same data). `rates.py Id1` still runs if wanted.

Not included: the certified BFF / finite-size EAT rate computation done after v13; it is not part of prl_v13.
