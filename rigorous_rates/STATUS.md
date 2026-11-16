# rigorous_rates STATUS (started 2026-10-22 01:30 HKT)
Goal: replace NPA-level-3 floating-point H_min numbers in prl_v13 (Results 2-3, Fig. 2, Table tab:rates, WBC panel, Liu-tilt, P_T appendix numbers) by rigorous two-sided intervals with proofs. Do not edit prl files.
## Identification of the entropy (from anc_submission_v7/rates_numerics/rates.py, npa.py)
- theta=pi/8, 3pi/16, 27.76deg tables/Fig2(a): H_min(A1|X=1,E) = -log2 P_guess, LOCAL guess of Alice's A1, single-score constraint g.q >= v*(g.P) (white noise p_v = v*P in correlator coords), or full statistics q = v*P. NPA level 3, real moment matrix, Clarabel.
- WBC panel Fig2(b) (Id2, phi=pi/6): H_min(AB|X=0,Y=0,E) (GLOBAL guess of (A0,B0)); H_min(A0|E) also in json.
- P_T appendix numbers: H_min(A0|E) and H_min(AB|00,E).
- 'best in cone' = random Dirichlet search + coordinate refinement over hexagon vertex weights (lower bound on sup).
## Log
- 01:30 setup; copied families.py.
- 01:35 Method (see PROOFS.md when written): Lemma 1/2 duality P_g <= max_e h_Q(Pi_e + mu g) - mu*omega (resp. + f, - f.p); h_Q bounded by exact-rational NPA certificates (cert.py: rational Z, exact LDL PSD check, levels 1+AB/2/3 lowest tight chosen); irrational data in arb ball arithmetic (families_arb.py, 200 bits). Upper bounds on H via explicit qubit attacks from LP column generation (rates_core.colgen), verified in arb; full statistics: mixture -> p_v/(1-t) with local remainder of weight t=1e-8 (Fine's theorem check).
- 01:37 First results theta=pi/8 v=0.999: tilted [0.776634,0.776634] (lvl 1+AB), CHSH [0.226659,0.226659] (closed form 0.2266586), paper-best [0.832439,0.832440] (lvl2), full [0.835087,0.835088] (lvl3), TRUE SUP over hexagon [0.832443,0.832443] (attack fooling all 6 vertices + certified g*). => in-cone gain over tilted CHSH for H_min(A1|E) is REAL (+7.2%).
- 01:39 run_cases.py cases.json (299 cases, priority: pi/8, 27.76deg, Id pi/6 AB, PT, dense, rest) in bg -> out/results.jsonl, certs/*.json, log out/run.log. Resumable (skips done ids).
  (restarted with OMP/BLAS threads=1 to stay light next to rates_v11 run_global.py)
- 01:50 fixed seeds for Id/PT full-statistics (target strategies) and vN lambda underflow; restarted both runners (resume). Added run_vn.py (vN attack upper bounds, exact rational rho via Cayley, arb entropies) -> out/vn_results.jsonl. First vN: pi/8 v=.999 tilted attack 0.970652, CONE attack (fools all 6 vertices) 0.970652 => sup over hexagon of H(A1|E) <= 0.97065 vs BFF tilted >= 0.96433: vN in-cone gain <= 0.0063 bits.
- 02:05 verify_all.py (standalone, flint only) + negative tests (tampered Z, beta, mu, v all caught). make_table.py -> TABLE.md. verify_vn.py. H_min run at 142/299; vN run done for pi/8 and 27.76deg (24 cases). vN: cone attack == tilted chord attack (target + local) at every point => vN in-cone gain <= attack - BFF(tilted) ~0.006-0.016 bits.
- 02:17 resumed after interruption: run_cases.py still alive (212/299, now ancillary th3pi_16/Id A0/Id pi/12); run_vn.py finished (24/24). Writing PROOFS.md sections 3/5.
- 02:27 H_min run at 299/299
- 02:27 ALL 299 H_min cases done and verified (verify_all.py ~2 s, 0 FAIL, max width 1.7e-6; every paper NPA-3 value consistent). PROOFS.md complete (Lemmas 1-6, Thms 1-5, Prop 6 vN, recommendations), TABLE.md, DERIVED.md. DONE.
