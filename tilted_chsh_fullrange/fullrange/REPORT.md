# Full-range tilted-CHSH normal cones (extension of Theorem 2 of prl_r4.tex): REPORT

Last updated 20:12 HKT (small-theta closure added; previous version kept as REPORT_pre_smalltheta.md); see the "Run status" section at the bottom (times in HKT, UTC+8). No prl file was edited.

Notation: theta = pi*x; s = sin 2theta, C = cos 2theta, r = sqrt(1+s^2); c0 = (1-s^2)/(1+s^2).
theta* = 0.20224059976992... pi, with s*^2 = x* the real root of x^3 - 8x^2 + 24x - 16 (s* = 0.95531250102563...).

## 1. Bottom line (what is certified)

**Theorem 2 extended (computer-assisted, same standard as prl_r4 Thm 2):** N_Q(p_theta) = C_theta for EVERY
theta in (0, pi/4] (open at 0 only: at theta = 0 the point p_0 is local, the midpoint of ++++ and +-++). More precisely:

* **Hexagon regime, theta in (0, theta*]:** the normalised section is the hexagon conv{V1, V2, V1', V3, V4, V3'}.
* **Octagon regime, theta in (theta*, pi/4]:** the normalised section is the octagon
  conv{Wa, Wb, V2, Wa', Wb', V3, V4, V3'} (cyclic order; Wa' := sigma(Wb), Wb' := sigma(Wa), see 3.3).
* At theta = theta* the two descriptions agree (Wa = Wb = V1, Wa' = Wb' = V1').
* At theta = pi/4 the octagon is **exactly** the Tsirelson octagon O of prl_r4 Thm 3 (BB's octagon), vertex by vertex
  (exact check, pi4_check.log): our result is consistent with, and at pi/4 re-derives, Theorem 3(a).
* Exposure statement as in Thm 2 holds on the whole range: open polygon = exactly the exposing functionals; an open edge
  supports a segment [p_theta, L_e]; a vertex supports a triangle. **Caveat at theta = theta* only:** V1 and V1' are tight
  on three deterministic points there (V1: ++++, +--+, -+--; V1': +++-, +-++, ----), so at theta* those two vertices
  support conv{p_theta*, L, L', L''} rather than a triangle. Everything else unchanged.

How (0, pi/16) was closed (section 2a, SMALLTHETA): asymptotic-in-t certificates, t = sin^2 2theta, verified in Arb
uniformly on t in [0, 3/20] (3/20 > sin^2(pi/8) = 0.1464, i.e. theta up to 0.0634 pi, overlapping the strip-certified
[pi/16, pi/4]): V1 by ALICE-overlap slicing (party swap; sgen.py), V3 by Bob-overlap slicing with exact t = 0 cancellation
(tgen3.py). No theta-strip runs are needed below pi/16 any more.

Combinatorics are proved EXACTLY on all of (0, pi/4] (not just on the certified range): hexagon on (0, theta*],
octagon on (theta*, pi/4]. Only the global (all-overlap) certificates of the non-analytic vertices limit the range.

## 2. Structure of the proof (same three parts as prl_r4 Thm 2) and what each part rests on

(0) **N_Q(p_theta) ⊆ C_theta** for all theta in (0, pi/4): unchanged from prl_r4 (slice cone + overlap stationarity +
    deterministic cuts are necessary conditions; analytic). At pi/4 this is the P_T computation of Thm 3.

(1) **Polygon (exact).** comb.py / comb.log. Every condition (4 mu_{a|x} > 0 and 14 non-facet margins 1 - g.L > 0 per
    vertex; 2 facet identities g.L = 1 per vertex) is an element of Q(s)(C, r). Its sign on an open s-interval is decided
    exactly: real roots of its norm over Q(s) and of all denominators are isolated with rational intervals (sympy), the
    element is evaluated in Arb on each isolating interval inside the range (must exclude 0) and at a rational point of
    every gap; a possible pole or zero raises. Results:
    - hexagon (V1, V2, V1', V3, V4, V3'): all 108 strict inequalities hold on the open interval (0, s*); facet identities exact.
    - octagon (Wa, Wb, V2, Wa', Wb', V3, V4, V3'): all 144 strict inequalities hold on (s*, 1); facet identities exact;
      at s = 1 all strict (Arb).
    - at s*: only the transition margins vanish (V1 vs -+--, V1' vs ----; Wa vs +--+, Wb vs ++++, etc.); they are >= 0 by
      continuity. All mu > 0 at s*.
    - plane_check.log (exact): each of the 10 vertices satisfies g = sum mu_{a|x} gamma_{a|x}, mu.m = 0, g.p_theta = 1.
    - **Vertex exhaustion (Lemma 4 generalised):** the proof of prl_r4 Lemma 4 applies verbatim with 8 vertices: the
      octagon's 8 facet points (-+--, +--+, +++-, ----, +-++, -++-, ---+, ++++), consecutive vertices share exactly one,
      all 8 distinct, each vertex tight on exactly its two and strict on the others => 8 distinct edges of the convex hull,
      so S = conv of the 8 vertices. For the hexagon on [pi/16, theta*] Lemma 4 is as in prl_r4, now with exact inequalities
      (prl_r4's Arb proof on [pi/16, 0.2022pi] is superseded/extended by comb.py on (0, theta*]). At theta* the extra
      tightness of V1, V1' does not affect the argument (only tightness among the six L_j is used).

(2) **Sufficiency: global certificates h_{K_c}(g) <= Phi_c(delta(c)) <= 1 for all c in [-1,1] at every vertex.**
    - V2, V4: analytic (prl_r4 App.), valid for all theta in (0, pi/4]. Side conditions (a0 sign, kappa < 1/2, c1 < c0,
      (beta + sgn(a0) kappa)^2 < nu^2, A2 > 0) re-proved exactly on (0,1) and at s = 1 (v24_conditions.log); closed
      forms match the vertex data exactly.
    - V1' , V3', Wa', Wb': by sigma (A1 -> -A1, B0 <-> B1), which fixes p_theta and each K_c:
      sigma(V1) = V1', sigma(V3) = V3', sigma(Wa) = Wb', sigma(Wb) = Wa' (exact, sigma_check.log).
    - V1, V3, Wa, Wb: piecewise-polynomial (per theta-strip, Chebyshev in c) certificates
      delta(c) = d0 + d1 (c - c0) + (c - c0)^2 [sum_j R_j T_j(c) + alpha + beta c], with exact forced data d0, d1 at c0
      (exact_c0.log: the identities F = F_c = 0 at c0 hold exactly in Q(s)(C,r) for V1, V3, Wa, Wb), optional exact PIN at
      c = +1 (Wa, Wb, V1 near theta*: v_{+|1}(1) = 0) or MID-PIN at c = -1 (V3 low: delta(-1) at the midpoint of two kink
      candidates; this only chooses delta, it adds no enclosure code). Verified in Arb (python-flint, 120 bits) by
      **vgen.py** (md5 02c3b44abc5f6bedbe89980db2a13d25, final version 17:17:40 HKT) over [-1,1] x strip, all 4 sign pairs:
      tests T0 (direct enclosure > 0), T1 (mean value), T2 (F_cc > 0 on hull(box, c0), using the exact c0 identities and
      mu > 0), T3/T3h (strict monotonicity to an endpoint where F >= 0); endpoint check at c = +-1 (mixed terms > 0 in Arb,
      the others equal 1 - g.L >= 0 by part (1)); exact area check that accepted boxes tile the strip. Strictness: every test
      gives F > 0 for c not in {-1, c0, 1}, which is what part (3) needs.
      Certified strips (all with the final vgen):

      | vertex | theta/pi range | file | log | result |
      |---|---|---|---|---|
      | Wa | [0.20224, 0.25] (24 strips, pin +1) | c_Wa.json | v_Wa_0_19_final.log + v_Wa_20_23.log | ALL OK |
      | Wb | [0.20224, 0.25] (24 strips, pin +1) | c_Wb.json | v_Wb.log | ALL OK |
      | V3 | [0.1875, 0.25] (24 strips) | c_V3h.json | v_V3h.log | ALL OK |
      | V1 | [0.1875, 0.20225] (8 strips, pin +1) | c_V1h.json | v_V1h.log | ALL OK |
      | V1 | [0.0625, 0.125] (8 strips) | c_V1l.json | v_V1l.log | ALL OK |
      | V3 | [0.0625, 0.125] (16 strips, mid-pin -1) | c_V3l.json | v_V3l.log | ALL OK |
      | V1, V3 | [0.125, 0.1875] = Theta | prl/anc/thm3_tilted (certs_v13.json, verify_thm3.py) | verify_V1_full.log, verify_V3_full.log | ALL OK (prl_r4 Thm 2) |
      | V1 | [0.03, 0.0625] | c_V1ll.json | v_V1ll.log | ALL OK (redundant: covered by 2a) |
      | V3 | [0.03, 0.0625] | c_V3ll.json (run as byte-identical c_V3ll_par.json) | v_V3ll.log | ALL OK (redundant: covered by 2a) |
      | V1 | [0.02, 0.03] | c_V1lll.json | v_V1lll.log | still running at 20:12 HKT; NOT NEEDED (covered by 2a) |
      | V3 | [0.02, 0.03] | c_V3lll.json (run as c_V3lll_par.json) | v_V3lll.log | still running at 20:12 HKT; NOT NEEDED (covered by 2a) |
      | V1 | (0, 0.0634] (t in [0, 3/20]) | smalltheta/c_V1s.json, c_V1s_ext.json | smalltheta/s_V1s_final.log, s_V1s_ext_final.log | ALL OK (section 2a) |
      | V3 | (0, 0.0634] (t in [0, 3/20]) | smalltheta/c_V3t.json, c_V3m.json | smalltheta/t3_V3t.log, t3_V3m.log | ALL OK (section 2a) |

      Coverage per regime: hexagon (0, theta*] needs V1, V3 on it: (2a) on (0, 0.0634 pi] + V1l + Theta + V1h, and
      (2a) + V3l + Theta + V3h (V1ll/V3ll give a second, independent cover of [0.03, 0.0625] pi). Octagon
      (theta*, pi/4] needs Wa, Wb, V3: c_Wa, c_Wb, c_V3h. (Strips may overhang the regime boundary; only the regime's own
      range is claimed. Wa/Wb strips start at 0.20224 < theta*/pi = 0.2022406.)

(2a) **Small theta (SMALLTHETA; files in /workspace/ndwu/smalltheta/, hashes HASHES_smalltheta.md5).**
    Variable t = s^2 = sin^2 2theta in [0, 3/20]; all vertex data are analytic in t (C = sqrt(1-t), r = sqrt(1+t)).
    Why strips fail as theta -> 0: the margins 1 - h_{K_c}(g) vanish like t (V3) and t^2 (V1) in the bulk, so one must
    divide out the vanishing order EXACTLY. Both verifiers do this algebraically (no cancellation in Arb) and cover t = 0
    itself (where the inequality is the trivial identity); the claim used is t in (0, 3/20].
    - **V1: Alice-overlap slices (sgen.py, md5 0dc77534..).** P = party swap preserves Q, h_Q(g) = h_Q(Pg), faces map to
      faces (as in Prop 6). Bob slices of g' = P g_V1 at p' = P p_theta have target overlap c0' = 0 for all theta, and
      the t = 0 slice vectors are nonzero axis vectors => single scale (Bob slicing of V1 provably needs an inner layer at
      c = 1 - O(t): at t = 0 strict convexity forces delta(c,0) = const but delta_c(c0) != 0). Exact swapped data
      (swapdata.py: slice generators, mu', forced d0', d1'; g' = P g_V1 and the c0'-identities verified exactly in
      Q(s)(C,r)). Certificate delta = d0(t) + d1(t) c + t c^2 (1/20, 0) (one rational vector!). G = F/t (pairs +-,-+,--),
      G = F/t^2 (pair ++, exact second-order cancellation, L(c,0) == 0 asserted). Data enclosed by exact Taylor
      coefficients at t = 0 (sympy) + Arb Lagrange remainder (K = 10; series-cap bug caught and fixed before the final runs).
      Tests T0, T2 (hull with c0' = 0, needs mu' > 0, checked per t-box), T3+- with the endpoint lemma (comb.py margins).
      Result: ALL OK on [0, 1/63] and [1/63, 3/20] (s_V1s_final.log, s_V1s_ext_final.log, re-run 20:11 HKT with the final
      files). Cross-check test_sgen.py (G vs direct F/t^k, 60 digits). Negative controls c_V1s_neg{a..e}.json all FAIL.
    - **V3: Bob slices (tgen3.py, md5 0d4ae8ff..).** At t = 0, u_{s0}(c,0) = 0 and w00 in {(1,0),(0,1)}; G = F/t written
      via exact hats X^ = (X(t)-X(0))/t in (t, gam, rho) (prep.py, hats_V3_V1.json; sqrt(1-t) = 1 - t gam, sqrt(1+t) =
      1 + t rho). Certificate delta = d0(0) + t[d0^ + d1^ (c-c0) + (c-c0)^2 (R(c) + mid-pin at c = -1)], R = R0 + (1+c)
      sum_j Rt_j T_j, two rational certificates: c_V3t.json on t in [0, 1/63], c_V3m.json on [1/63, 3/20] (same R0,
      Rt refit by fit3hi.py). Structure found: G_{-+}(-1,t) == 0 for all t (pin) with a c-layer of width ~t^2; at
      c = +1 the pair (+-) is tight (deterministic facet); corner (c,t) = (-1,0) degenerate (P = M = 0).
      Tests: T0/T0d/T1/T2/T3/T3h; T4 (G_ccc < 0 on [c0,1], G(1) >= 0, G(c0) = G_c(c0) = 0 => G >= 0); E- (G_c > 0 on
      [-1,cb], bound on d/dc of the slice norm: Hadamard at the corner for t-boxes from 0, and for ta > 0 the mean-value
      bound P^2/M <= 2P(-1)^2/M(-1) + 2(1+cb)P_c^2/min M_c, tgen3); endpoint check at c = +-1 where the mixed-sign terms
      are EXACTLY O(t^2) (c = +1) or O(t) (c = -1): certified by exact rational Taylor coefficients at t = 0 (fmpq_series)
      + Arb remainder over [0, tb] (tgen2/3; test_tgen2_ep.py: exact coefficients contained in the Arb series).
      Result: c_V3t ALL OK (5 strips, 171160 boxes, t3_V3t.log), c_V3m ALL OK (7 strips, 731954 boxes, t3_V3m.log), one
      code version (tgen3). Negative controls: c_V3t_negep.json (endpoint, mixed term ~ -0.1 t^2) FAILS; c_V3t_negF{a,h}.json
      (G < 0 in the edge region c in (-0.9995, -0.8)) FAIL, and every edge box containing float-negative G is rejected
      (negE3.py, incl. the t = 0 corner box). G-identity cross-check vs direct mpmath F/t (test_tgen3.py).
    - Strictness for part (3): all accepting tests give G > 0 for c not in {-1, c0, 1} and t > 0 (T4, E- strictly
      monotone/convex), so exposure follows as below; for V1 the argument is run for P g_V1 in the swapped labelling and
      transported back by P (P maps Q onto Q and faces onto faces).
    - Swap trick for V3 does NOT help: P(V3) has mu'(0) = (0,0,1/2,1/2), u00_+ = w00_+ = 0 and its Alice-slice margin is
      ~t^3 with an O(t) layer at c' = 0 (marg_sw.py), worse than Bob slicing.

(3) **Exposure.** As in prl_r4: for g in the open polygon all 16 deterministic margins are strict (part 1), so
    h_{K_{+-1}}(g) < 1; by sublinearity and the strict certificate inequalities, C*(g) = {c0}; all four mu_{a|x} > 0 and the
    complementary-slackness computation (theta-generic algebra) gives the unique maximiser p_theta. Edges/vertices: faces
    from the tight deterministic points (caveat at theta* above).

## 3. Soundness checks performed (this session)

3.1 pinned_term (new code path near pi/4, where the pinned vector's G_1-norm -> 0): formula re-derived by hand (correct:
    N^2 = (1+c')/2 P^2 + (1-c')/2 M^2, X = (1-c')/N in [0, sqrt(2(1-c'))/|M|], |Y| <= sqrt((1+c')/2)).
    - Negative controls c_Wa_negpin_{a,b,c,d}.json (Wa strip 22, R perturbed by k (c-1)(1+c)^6/64 so the pin is unchanged;
      float F < 0 on c in [0.98..0.9995, 1), down to -4e-8, on the pinned pairs): vgen FAILS on all four.
    - contain_pin.py: 4000 point-value and finite-difference (300-bit, h = 1e-40) containment tests of the F and F_c
      enclosures on boxes near c = 1 (strips 20-23, 12; manual enclosure picked 764 times): 0 violations.
3.2 Mid-pin certificates: negative controls c_V3l_negmid_{0,1}.json (float F ~ -1e-2 on c in [-0.98,-0.6]): vgen FAILS
    inside the negative region.
3.3 Earlier controls (previous worker, after the 16:48 soundness fix): perturbed Wa certificate fails; correct strip passes.
3.4 Version discipline: every result in the table was produced by vgen.py md5 02c3b44a... (runs started after 17:17:40);
    all pre-16:48 results are invalid and kept only in invalid_runs/. HASHES.md5 lists hashes of all certificate/data files.
3.5 Exact symbolic checks: sigma_check.log, plane_check.log, exact_c0.log, comb.log, v24_conditions.log, pi4_check.log.

## 4. Endpoints

* **theta = pi/4 (CHSH):** p_{pi/4} = P_T, c0 = 0. The 8 octagon vertices are g(Lambda = 1/2, d) with |d| = rho = 1 - 1/sqrt2
  at angles -90, -135, 180, 135, 90, 45, 0, -45 degrees (Wa, Wb, V2, Wa', Wb', V3, V4, V3'): exactly the regular Tsirelson
  octagon O of prl_r4 Thm 3, in the same cyclic order. So the octagon family deforms continuously from the Tsirelson
  octagon (pi/4) and collapses to the hexagon at theta* as the edges [Wa, Wb] and [Wa', Wb'] shrink to points.
  (Our pi/4 certificates for Wa/Wb are pinned at c = +1; Thm 3's are different but both valid.)
* **theta -> 0:** p_theta -> (a0, a1, b0, b1, e00, e01, e10, e11) = (1, 0, 1, 1, 1, 1, 0, 0), the midpoint of the local
  deterministic points ++++ and +-++ (not extremal), so the theorem is stated on (0, pi/4]. The small-theta certificates
  (2a) are uniform in t = sin^2 2theta down to t = 0, so no strip near 0 is needed. At t = 0, g(V3) = (0,0,1/2,1/2,0,0,1/2,
  -1/2) is a local functional (h_{K_c} = 1 for all c) and g(V1) -> (2/5,1/5,0,1/5,1/5,1/5,1/5,-2/5).

## 5. Open gaps

* None for the range (0, pi/4] at the computer-assisted standard of prl_r4 Thm 2. theta = 0 is excluded (local point).
* The certificates are computer-assisted (Arb branch-and-bound), like prl_r4 Thm 2; V2, V4 are analytic.
* The open-ended statement "N_Q(p_theta) ⊆ C_theta" is the prl_r4 analytic argument; not re-audited here.
* V1lll/V3lll ([0.02, 0.03] pi) Arb runs are no longer needed; V1ll/V3ll (passed) are a redundant second cover.

## 6. Files (all in /workspace/ndwu/fullrange/)

oct_data.json (exact Wa, Wb, Wa', Wb' data), prl/anc/thm3_tilted/thm3_data.json (V1..V4 data); c_*.json (certificates);
vgen.py (verifier); v_*.log (verification logs); comb.py (+kfield.py, kfun.py) / comb.log; exact_c0_all.py / exact_c0.log;
sigma_check.py/.log; plane_check.py/.log; v24_conditions.py/.log; pi4_check.py/.log; contain_pin.py; neg_pin_design.py,
neg_mid_design.py, c_Wa_negpin_*.json, c_V3l_negmid_*.json, v_negpin_*.log, v_negmid_*.log; design2.py, midpin.py (design
tools; not needed for verification except midpin.cands, which vgen uses to define the mid-pinned delta); STATUS.md
(chronological log); HASHES.md5.

## Run status
* 17:59 HKT: V1ll running (2/16 strips done, all OK so far); then V3ll, V1lll, V3lll (runall.sh).
* 20:12 HKT: V1ll ALL OK, V3ll ALL OK. Small-theta: V1 (sgen) and V3 (tgen3) ALL OK on t in [0, 3/20]. Final range
  (0, pi/4]. V1lll / V3lll still running (not killed), not needed.
