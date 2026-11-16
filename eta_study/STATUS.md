# eta_study STATUS (started 2026-10-22 09:25 HKT)
Goal: certified DI randomness per round under symmetric detector efficiency eta (no-click -> +1), tilted-CHSH /
Thm-2 cone witnesses vs CHSH (max-ent and optimally tilted) vs full statistics; H_min (two-sided rigorous) and
von Neumann (BFF certified) + EAT; rigorous separations where they exist. Do NOT edit prl files.
Reuse: ../rigorous_rates (cert.py, rates_core.py, qubit.py), ../rates_v11 (bffcert.py, eat.py).
## Log
- 09:25 setup.
- 09:28 core.py (loss model, parametrised NPA P_guess: PGm rows>= / rows==), cone.py (C_theta polygon from slice cone +
  16 deterministic cuts; reproduces hexagon at pi/8 and BB octagon at pi/4). Families: tilted = constraints {CHSH>=, <A0>=}
  (= best alpha<A0>+CHSH by SDP duality), cone = {V_i(theta')>= all vertices} (= best witness in C_theta'), full.
- 09:33 opt.py split per family (4 light procs, ~150 MB each): out/opt_v1.0_{chsh,tilted,full,cone}.jsonl (+ _old.jsonl
  from first run). CHSH done: chsh_opt S_max(eta) > 2 down to eta=2/3 (Eberhard), max-ent dies at 0.8284.
  First points: eta=.995 H_min CHSH .7386, tilted .7449, cone .7608, full .7892.
- 09:40 cert_eta.py works (exact-rational chain). eta=0.98 v=1 PROVED: every CHSH protocol H_min(A1|E) <= 0.525646;
  tilted (alpha<A0>+CHSH) >= 0.535303; Thm-2 cone witness >= 0.558502 (two-sided widths ~2e-8). certs/eta0.98_v1.0.json
- 09:44 bffexact.py: BFF certification with EXACT rational PSD check (LDL, fmpq) replacing float eigen margin; agrees
  with rates_v11 to 1e-7. bff_eta.py (vN + EAT) uses it. verify_eta.py standalone verifier (H_min) written; passes eta=.98.
- 09:50 opt v=1 done for chsh/tilted/full (cone finishing); opt v=0.99 running (3 procs). cert_eta.py running over
  eta grid (v=1) -> out/cert_results.jsonl, certs/; fixed best_alpha (grid+refine; Brent alone failed at .89/.90).
  bff_eta.py v=1 running (vN, exact-PSD BFF, EAT) -> out/bff_v1.0.jsonl.
  Float landscape v=1 (H_min A1): ratio to best CHSH: tilted ~1.01-1.05, cone ~1.03-1.11, full ~1.07-1.14; vs max-ent
  CHSH unbounded below 0.828.
- 09:57 refine.py (warm starts incl. CHSH-optimal device => family >= CHSH guaranteed) running for v=1; load_opt keeps
  best record per key. Found near-threshold (eta<=0.72) and cone (eta<0.75) optimisation failures in first pass.
  vN: at H_min-optimised devices cone/full vN < CHSH vN (devices wrong for vN) -> killed bff_eta.py, wrote vn_opt.py
  (NM over device+alpha with BFF m=6 objective, final exact-PSD BFF m=8, EAT). Running 2 procs (v=1).
  eta=.9 quick check: CHSH vN .3049 (closed form, best device); tilted at tilted device BFF m6 .3228.
- 10:12 v=1 H_min PROVED separations (tilted AND cone > every CHSH protocol) at all 27 grid points eta in [0.70,0.995]
  (verify_eta.py). Re-certifying eta 0.68-0.72 with refined optima (out/cert_v1_low.log). variants_check.py: all 8 CHSH
  sign variants x 4 no-click binnings: canonical (+++-, no-click->+1) is optimal (certified) -> out/variants_check.json.
  vn_opt.py restarted with m=4 search (m=6 too slow), 2 procs, eta grid .995,.98,.95,.92,.89,.86,.83,.80,.77,.74,.71.
  refine.py v=0.99 started.
- 10:16 RESUMED after interruption. Alive: refine v1 (at cone .90), refine v.99, vn_opt x2 (slow: ~12 min/family/eta).
  Done: cert low-eta (.68-.72 re-certified), opt AB11 (two-party, 7 etas). vN restarted with reduced scope:
  vn_v1_c/d (chsh+tilted, .95,.92,.89,.86,.83,.80,.77,.74), vn_v1_e (full at .98,.90,.83); existing vN: .995,.98 tilted.
  Started cert v=0.99 (16 etas) and cert AB11 (6 etas).
  Two-party float: tilted ~= CHSH (alpha->0 for eta>=.9), cone +2..6%, full +15..30% over best CHSH.
- 10:21 neg_tests.py (8 tamperings; 7 rejected, 1 = alpha change yields a valid much weaker bound -> sound).
  v=1 A1: proved at every grid eta in [0.69, 0.995] (0.68: not proved, margins below cert precision).
  AB11 (two-party, v=1): cone witness PROVED > every CHSH at eta .75,.80,.85,.90,.95,.98; tilted family NOT (alpha->0).
  BUG FOUND+FIXED: chsh_side ignored v (S_up must be c0 + v*beta); v=0.99 certs deleted and rerunning; verifier updated.
  vN v=1 so far: .995 tilted .92799 vs CHSH .91944; .98 .78192 vs .76026; .92 .40955 vs .38635; .89 .28681 vs .26947.
- 10:36 v=0.99 re-certified with fixed S_up=c0+v*beta: tilted AND cone PROVED > every CHSH at all grid eta in [0.75,0.995]
  (0.755 has no opt record -> only CHSH bound). verify_eta.py: 56 records, all claims re-verified (out/verify_eta.log). neg_tests rerun: same 7/8 + sound degradation.
  Remaining: vN runs c,d,e (v=1), then tables, figs, REPORT.md, PRL paragraph.
- 10:56 vN v=1 complete (out/vn_v1.0.jsonl): tilted > every-CHSH vN bound at eta .77-.995 (8 pts, +0.9..+6.4%);
  eta .74: numerical m8 .01298 > .01253 but certified .01241 -> not claimed. Full vN .98 .810, .90 .361, .83 .135
  (EAT n/a: ill-conditioned kappa). EAT n=1e10 tilted > CHSH at .77-.995; n=1e8 tilted < CHSH at .995 and .77.
  make_figs.py: new fig_eta_prl (PRL figure), vN gain panel, legend fixes. REPORT.md assembled by build_report.py
  from report_src/*.md + out/tables.md. prl_paragraph_eta.tex written (prl files untouched).
  Running: vn_opt.py 0.99 0.94,0.9 chsh,tilted (out/vn_v99.log).
