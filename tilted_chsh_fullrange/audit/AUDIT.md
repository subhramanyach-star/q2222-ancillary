# Independent audit of the full-range Theorem 2 (tilted CHSH, theta in (0, pi/4])

Auditor: executor subagent, 21 Oct 2026, 20:15-20:45 HKT. Audit artefacts: `fullrange/audit/` (scripts, logs, md5 snapshot).
No certificate, verifier or data file was modified (md5 before/after identical, `audit/md5_before.txt`).

## Verdict
**PASS, with fixable caveats about the write-up (no mathematical gap found).** The claim
N_Q(p_theta) = C_theta for every theta in (0, pi/4] (hexagon on (0, theta*], octagon on (theta*, pi/4], Tsirelson octagon at
pi/4), with the stated face structure (open polygon = exposing set; open edge -> segment; vertex -> triangle; tetrahedron for
V1, V1' at theta* only), is supported at the same computer-assisted standard as prl_r4 Thm 2.

Caveats (all are presentation issues; none changes the statement):
1. **Exposure argument must not use Bob-slice equality sets of V1 for small theta.** For theta < pi/16 the V1 certificate is in
   Alice-overlap slices, so C*(V1) in Bob slices is not controlled. The face statements still follow (section (a) below), but
   the paper's proof should say: (i) h_{K_c}(V1) <= h_Q(V1) <= 1 for all Bob overlaps c (K_c is a subset of Q), so Lemma
   tools(b) still applies to every combination containing V1; interior and edge functionals get C* from the *other* vertex
   (V2, V4 have {-1, c0}; V3, V3' have a subset of {-1, c0, 1}); (ii) the vertex faces of V1, V1' come from the swapped problem.
2. **tgen3 corner hypothesis is asserted only in a docstring.** The Hadamard branch of `edge_dc` (t-boxes starting at t = 0)
   needs u^_-(c = -1, t = 0) = 0. `CORNER_ZERO_CHECKED` is defined but never used. I verified it exactly
   (`audit/corner_check.py`: P(-1,t), M(-1,t) have zero constant term, both certificates). Ship that script with the ancillary files.
3. **Missing negative control for tgen3's T4/T2 path** (c in (c0, 1), t small). Added: `audit/c_V3t_negT4_{0,1}.json`
   (R -> R + 16(1-c)v, float G < 0 at c ~ 0.9993, t = 1/1100, endpoint values unchanged): tgen3 rejects both.
4. REPORT section (3) words exposure as "by sublinearity ... C*(g) = {c0}". That is correct for the open polygon (it uses V2/V4),
   but the vertex/edge statements for V1, V1' at small theta need caveat 1.

## (a) Party-swap face transport for V1 (sgen.py)
* Swap P: (a0,a1,b0,b1,E00,E01,E10,E11) -> (b0,b1,a0,a1,E00,E10,E01,E11) is correct, maps Q onto Q, so
  h_Q(g) = h_Q(Pg) = max_c h_{K_c}(Pg) and F_Q(Pg) = P F_Q(g). `swapdata.py` recomputes the slice generators at p' = P p_theta
  (c0' = 0, Gram = I), solves for mu' and asserts g' = P g_V1 exactly in Q(s)(C,r); the c0'-identities
  (||u|| = 2 mu'_{s0|0}, ||w|| = 2 mu'_{s1|1}, F = 0, F_c = 0) are checked exactly. mu' > 0 is checked per t-box in sgen (Taylor + remainder).
* G-formulas re-derived by hand: G = F/t for (+-),(-+),(--) (axis vectors u00, w00, ||u00||_c = |u00|); for (++) the second-order
  formula Lhh - ||u^||^2/(||u||+a) + q_u<u00,u^>/(a(||u||+a)^2) + (w terms) equals F/t^2 given L(c,0) == 0 (asserted). Correct.
  Lagrange remainders (series coefficient K over the ball [0, tb]; Taylor form in c with C(K,j) B_K) correct. cap = 64 assertion present.
* Endpoint handling in sgen: at c = +-1, F = T - |u.f| - |w.f| exactly; when the two signs are certified equal, F = 1 - g'.L' with
  L' deterministic; 1 - g'.L' = 1 - g.(P^-1 L') >= 0 by comb.py (all 16 margins of V1). Sound (no H2 needed).
* In-slice uniqueness for the swapped problem (Lemma tools(d) applied to (Pg, Pp)): CS forces z_{a|x} = (1 + a q_x) v_{a|x}; after
  eliminating Bob' marginals, q is unique iff det(d_0, d_1) != 0, d_x = (v_{+|x} - v_{-|x})/2. Exact: det = -(1+s^2)/2 != 0
  (`audit/swap_uniqueness.py`). So argmax_{K'_0} Pg_V1 = {Pp_theta}.
* Strictness: T0, T2 (with G = G_c = 0 at c' = 0 and G_cc > 0), T3+- give G > 0 for c' not in {-1, 0, 1}. Maximisers at c' = +-1 are
  local (Alice's observables coincide up to sign) and lie in conv of the tight deterministic points with A0 = +-A1. Face formula
  (applied to the swapped problem) gives F(Pg_V1) = conv{Pp, P(++++), P(+--+)}; transported by P: F(V1) = conv{p, ++++, +--+}. V1' by sigma.
* Interior / edges: F(sum lambda_k U_k) = intersection of F(U_k) (any convex set, h(U_k) = 1 at common p); each L is tight on
  exactly two adjacent vertices (comb.py), so intersections are {p} (interior) and [p, L_e] (edges). Transfer is genuine.
* theta* caveat: V1 tight on ++++, +--+, -+--; the middle point (-+--) is the Wa/Wb shared point and is tight on no neighbour, so
  edges stay segments; Gram det of (L_i - p_theta*) = 213.26 > 0 (`audit/tetra_check.py`), so the vertex face is a genuine tetrahedron.

## (b) Coverage (exact rational strip ends checked; no internal gaps in any certificate file)
| vertex | needed on | covered by |
|---|---|---|
| V1 | (0, theta*] | sgen t in [0, 3/20] (theta <= 0.06329 pi) ; c_V1l [1/16, 1/8] ; prl_r4 Theta [1/8, 3/16] ; c_V1h [3/16, 809/4000 = 0.20225] > theta*/pi = 0.2022406 |
| V3 | (0, pi/4] | tgen3 c_V3t [0, 1/63] + c_V3m [1/63, 3/20] ; c_V3l [1/16, 1/8] ; Theta ; c_V3h [3/16, 1/4] |
| Wa, Wb | (theta*, pi/4] | c_Wa, c_Wb [632/3125 = 0.20224, 1/4] (starts below theta*) |
| V2, V4 | (0, pi/4] | analytic; side conditions exact on (0,1] (v24_conditions.log) |
| V1', V3', Wa', Wb' | as above | sigma (exact: sigma(Wa) = Wb', sigma(Wb) = Wa') |
sin^2(2 * pi/16) = 0.1464 < 3/20, so the small-theta pieces overlap [pi/16, ...]. pi/4 endpoint: Wa (strip 23 ends at 1/4), Wb,
V3h end at exactly 1/4; comb.py covers s = 1 by Arb (all strict). Overhanging strips are harmless: vgen does not use H1/mu > 0
except through comb.py on the claimed regime. Polygon combinatorics: comb.log, exact on (0, s*) and (s*, 1), Arb at s*, s = 1.
pi4_check.log: octagon at s = 1 equals the Tsirelson octagon vertex by vertex. Redundant covers V1ll/V3ll [0.03, 0.0625] ALL OK.

## (c) Verifiers, reruns, hashes, controls
Fresh-process reruns from the archived certificates (`audit/run_audit.sh`, finished 20:22 HKT) reproduce the archived logs
box-for-box:
| run | result | archived |
|---|---|---|
| vgen c_Wa strip 12 | OK, 2456 boxes | 2456 |
| vgen c_Wa strip 23 (pi/4, pinned_term path) | OK, 1667 | 1667 |
| vgen c_V1h strip 7 (to theta*) | OK, 3543 | 3543 |
| vgen c_V3h strip 23 (pi/4) | OK, 3558 | 3558 |
| vgen c_V3l strip 0 | OK, 95491 | 95491 |
| sgen c_V1s (all 3 t-strips) | OK, identical counts | s_V1s_final.log |
| sgen c_V1s_ext strip 3 [1/10, 3/20] | OK, 2048 c-boxes | identical |
| tgen3 c_V3t strip 0 [0, 1/1024] (corner) | OK, 33758 | 33758 |
Negative controls (must FAIL; all do): vgen generic c_Wa_neg (archived), pinned c_Wa_negpin_a (rerun) + b-d (archived), mid-pin
c_V3l_negmid_0 (rerun) + _1; sgen c_V1s_nega (rerun) + b-e (archived; c -> 1 and T2 region); tgen3 endpoint c_V3t_negep (rerun),
edge E- c_V3t_negF{a,h} (archived) + negE3.py; tgen3 T4/T2 region c_V3t_negT4_{0,1} (new, this audit). Each code path has a control.
Hashes: vgen 02c3b44a..., sgen 0dc77534..., tgen3 0d4ae8ff..., all certificate/data hashes in HASHES.md5 and
HASHES_smalltheta.md5 verify (c_V3ll/c_V3lll were renamed *.queued_elsewhere/_par.json, byte-identical).
vgen code read in full: F = T - ||u|| - ||w|| (Theorem 1 dual); TP exact centred expansion; mean-value theta enclosures; T0/T1/T2/
T3/T3h logic; endpoint check (non-mixed terms delta-free = 1 - g.L, mixed > 0 in Arb; pinned pairs all deterministic since the pinned
vector vanishes exactly); pinned_term bounds re-derived; exact area check. Endpoint lemma use at c = +-1 is sound in all three
verifiers (vgen via H1 + strict H2; sgen/tgen3 via certified sign equality, which makes F(+-1) a deterministic margin exactly).
tgen3 exact Taylor coefficients: fmpq_series from the exact hats (prep.py checks the hats against the original expressions at 50 digits);
`audit/test_tgen3_ep.log`: 0 mismatches against the Arb series at t = 0, endpoint passes from t = 0. test_sgen: |G - F/t^k| <= 2.7e-18.
T4 logic (G_ccc < 0 on [lo, 1], G(c0) = G_c(c0) = 0, G(1) >= 0 => G > 0 off c0) re-derived: correct for both sides of c0.
E- bounds (direct, mean-value, Hadamard) re-derived: correct given caveat 2.

## Not re-audited here
N_Q(p_theta) subset of C_theta (prl_r4 analytic argument), Theorem 1, Lemma 4 logic (only extended to 8 vertices; the argument is the
same), the prl_r4 Theta certificates (audited in earlier rounds).
