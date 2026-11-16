# rates_v11 STATUS (started 2026-10-21 22:40 HKT)
Goal: rigorous BFF H(A1|X=1,E) lower bounds (dual-certified) for CHSH vs tilted vs best-in-cone vs full stats under Werner noise; EAT finite-size rates (Liu et al. 2021 Nat. Phys. EAT, Thm 3 + spot-checking Lemma 3); figure + REPORT.md.
No MOSEK on box; using Clarabel (cvxpy) with tight tolerances + own dual post-processing (bffcert.py).
## Log
- 22:40 set up; copied families.py (hexagon vertices) from anc_submission_v7/rates_numerics.
- 22:58 bffcert.py works: SDPA (sdpa-python, pip --user) + multi-trial (Z scalings 1,0.5,0.75; SDPA params; adaptive row scaling on 2 smallest-t nodes) + exact-rational dual certification. CHSH cross-check: certified <= analytic everywhere, deficit 1e-4..2e-3 (m=8), 7.8e-3 at S=2.8256 (Radau t=1-node truncation ~1/(2m^2 ln2)); m=12 reduces it to 2.8e-3. Clarabel stalls (AlmostSolved) on these SDPs; COPT size-limited.
- 22:58 eat.py written (Liu et al. 2021 Thm 3 + Lemma 3, statistics-agnostic V, Bernstein completeness).
- 22:59 run_main.py started in background (out/run_main.log, out/bff_main.jsonl; resumable): tilted & full, 4 thetas x 14 v, m=8 (+m=12 for v>=0.97).
- 23:08 queued: chain.sh (run_cone.py after run_main) -> chain2.sh (run_main again for added v=0.994) -> chain3.sh (check_chsh.py, x*=1, m=8,12). Markers out/cone.done, out/main2.done, out/chsh.done. analyze.py = tables + EAT; figure TBD (make_fig.py).
- 23:23 resumed after interruption: background chain intact (run_main at pi/16; 86 records). Findings so far: pi/8 tilted dies at v~0.966=(C+r)/2 as predicted; gains vs CHSH on same state +88% (v=.99), +143% (v=.999); 27.76deg +56% (v=.99).
- 23:35 parent steering: theta-envelope + H_min vs vN. Added run_env.py (theta=0.1pi,0.22pi,pi/4 BFF) and run_hmin.py (NPA lvl-2 H_min, all 7 thetas; numerical only), queued in chain4.sh after chsh.done -> out/env.done.
- 23:58 cone search running: pi/8 v=.99 best-in-cone = tilted CHSH itself (lam->tilted, hL=(C+r)/2).
- 00:22 resumed (2nd interruption). Killed chain2/3/4; new chain5.sh: after cone.done -> run_env.py + run_hmin.py (env.done) -> run_main.py v=0.994 (main2.done) -> check_chsh.py (chsh.done). All scripts checkpoint per point (jsonl/json append). Cone result: best sigma-symmetric in-cone expression = tilted CHSH (hL=(C+r)/2) at pi/8 all v and 27.76deg (no gain over tilted).
- 00:34 cone.done; run_env.py running; analyze.py (EAT) running in background -> out/eat_table.json, out/table_eat.md. EAT preview: 27.76deg v=.99 n=1e10: CHSH .49, tilted .73, full .77.
- 00:46 env.done, hmin done; envelope.py -> out/table_envelope.md, out/table_theta_detail.md: theta=pi/4 (CHSH) wins the white-noise envelope at every v for both vN and Hmin. EAT table done for v in {.999,.99,.985,.98}; v=.994 pending main2.
- 00:50 main2.done (v=0.994 all theta). analyze.py made incremental (EAT keys checkpointed, cone=max(cone,tilted), pi/4 added). Running analyze.py -> out/analyze2.log
- 00:58 chsh.done -> out/table_chsh_check.md; analyze.py EAT rerun in progress; make_fig.py now 2x2 with envelope panel
- 01:20 resumed (3rd). analyze.py EAT rerun finished (01:17; table_eat.md incl. v=.994, pi/4). Task 2 started: bffcert2.py (global H(AB|X=x*,Y=2,E), 2x3, projector letters, exact-rational cert; tested: max-ent v=.999 full 1.9605, v=.98 1.667), targets2.py (B2 chosen to max ideal H(AB)), run_global.py running in bg -> out/global/bff_global.jsonl (marker out/global/global.done), ~1 min/point.
- 01:35 make_fig.py updated (panel d now includes max-ent CHSH EAT curve at v=.994); REPORT_part1.md (Sec 0-6) written; global run ~70 s/point (84 points); npa_hmin2.py/run_hmin2.py (global H_min cross-check) ready, to run after global.done.
- 02:01 global run in progress: 36/84 points
- 02:02 resumed (4th): run_global.py alive (37/84), chain6.sh waiting for global.done -> run_hmin2.py. Other workers' SDP jobs present (run_cases/run_vn/composite_run); not starting extra heavy jobs.
- 02:28 global 58/84
- 02:29 upper_global.py -> out/global/upper_global.json (device-dependent purification attack U = H_Shannon - S(rho_v); rigorous upper bounds on DI global H for tilted sources). 0.22pi U(.99)=1.90 > max-ent cert 1.815 -> no rigorous strict separation there.
- 03:02 global.done (84/84), hmin.done; analyze_global.py -> out/global/table_global.md, fig_global_v11.pdf
- 03:03 DONE: REPORT.md (= part1 + part2, Sec 0-8), fig_rates_v11.pdf (2x2), fig_global_v11.pdf, out/global/table_separation.md. Answer: no tilted advantage in global randomness at any tested v.
