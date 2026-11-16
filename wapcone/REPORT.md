# WAP normal cones, combined family, and the per-round optimum under loss
Scope: H_min(A1|X=1,E) per round, symmetric detector efficiency η (no-click → +1), v = 1. Each protocol uses the same detectors, each with its own best state and settings.
Labels: **[cert]** = two-sided interval from cert_eta.winner (exact-rational NPA duals + arb, plus an explicit qubit attack), width ≤ 2e-7. **[num]** = float NPA level 3 (level 4 agrees to 5e-8 at η = 0.90). Devices were optimised numerically (Nelder–Mead), so every row is a lower bound on what its family can do.
Files: fam.py, ceil4.py, recover98.py, cont.py, wp98.py (optimisation) → fam.jsonl. cert_wc.py → certs/, cert.jsonl. wcones.py → cones.json. Log: STATUS.md.

## Table (H_min per round)
| η | every CHSH protocol ≤ [cert] | WAP asym-CHSH [cert] | our polygon C_θ′ [cert] | combined: polygon ∪ WAP-cone ∪ tilted [cert] | best witness found [cert] | gain of best over CHSH | numerical optimum over devices (L3) [num] |
|---|---|---|---|---|---|---|---|
| 0.98 | 0.525646 | 0.580212 | 0.558502 | 0.587046 | **0.597368** | +13.6% | 0.597369 |
| 0.95 | 0.330204 | 0.366937 | 0.361305 | 0.374280 | **0.385336** | +16.7% | 0.389044 (not certified) |
| 0.90 | 0.168312 | 0.182071 | 0.186183 | 0.186750 | **0.197767** | +17.5% | 0.200222 (not certified) |
| 0.85 | 0.081930 | 0.085560 | 0.089927 | 0.090009 | **0.095105** | +16.1% | 0.095105 |
| 0.80 | 0.033444 | 0.033995 | 0.036046 | 0.036053 | **0.037864** | +13.2% | 0.037864 |
| 0.75 | 0.009013 | 0.009041 | 0.009465 | 0.009522 | **0.009847** | +9.3% | 0.009847 |

Notes on the table:
- The WAP values are higher than in litcheck because the device was re-optimised at NPA level 2. Level 1+AB is not tight for WAP witnesses.
- How the "best witness" column was obtained:
  - η = 0.98, 0.85, 0.80, 0.75: the witness comes from the level-3 full-statistics dual at the level-3-optimal device.
  - η = 0.95 and 0.90: the witness comes from the level-2-optimal device. There the level-3 dual is ill-conditioned: it is dominated by −4P(−1,−1|1,1), a Hardy-type near-zero, so the certificate came out loose (H_lo of 0.016 and 0.011). Those two points are therefore numerical only for now.
- At η = 0.98 the optimum is a near-maximally-entangled branch (θ = 0.778 rad), found from a random start. The partially entangled branch gives 0.594950. At η = 0.95, continuing from the near-max-ent branch gives only 0.365114, so the partially entangled branch wins there.

## Cones (step 1)
- **Ideal WAP targets.** The target p_a is: |Φ+⟩, A1 = Z, A0 = X, B_y = cos φ Z ± sin φ X, tan φ = 1/a, with a(η) = 0.73–1.00.
  - C(p_a) is an octagon at every η, and all four slice cones are active.
  - The WAP functional lies in C(p_a) (residual ≤ 3e-10).
  - N_Q(p_a) ⊆ C(p_a) is proved (paper Lemma 3).
  - Equality is numerical: max_c h_{K_c}(V) − 1 ≤ 0 to within 2e-9 for every vertex. The correlator-only part of the cone is exact (general/REPORT Thm A).
- **Lossy WAP devices.** The cones have 5 vertices (η ≤ 0.80), 7 (0.85–0.90) and 8 (0.95–0.98), each numerical to 2e-10.
- **Best witness in the WAP cone** ("wapcone", target a′ free) adds only +0.0001–0.0006 over WAP itself.
- **Combined family** (conic hull of polygon(θ′) ∪ wapcone(a′) ∪ {CHSH, ±A0}, device re-optimised): beats both WAP and the polygon at every η, but only by +0.0001–0.0068. It stays well below the best witness.

## The normal cone of the observed point (main finding)
"wappush" is the best witness in N_{Q_η}(p_obs), the normal cone of the lossy set at the observed behaviour itself, with the device optimised. Numerically it reproduces the full-statistics optimum at the same device: to 1e-6 at η = 0.98, 0.90, 0.85, 0.80 and 0.75, and to 1.3e-5 at η = 0.95. This is explained by the following statement.

**Proposition (elementary; proved, given its hypotheses).** Notation:
- L_η is the (affine, invertible) loss map, and Q_η = L_η(Q) is the set of honest behaviours at efficiency η.
- For a witness g and a score ω, P(g,ω) = sup{P_guess(q) : q ∈ Q, g·q ≥ ω} is the adversary's best guessing probability (no detector assumptions). It is non-increasing in ω.
- F(p) = min_g P(g, g·p) is the best single-witness bound at p. By SDP duality this equals the full-statistics bound (Nieto-Silleras–Pironio–Silman; Bancal–Sheridan–Scarani).

Suppose p* ∈ Q_η minimises F over Q_η (that is, maximises H_min over honest devices), g* attains F(p*), and ω ↦ P(g*,ω) is strictly decreasing at ω* = g*·p*. Then g*·p* = max_{p∈Q_η} g*·p, i.e. g* ∈ N_{Q_η}(p*). Equivalently, the pull-back g*_η lies in N_Q(L_η^{-1}p*), which is the ideal normal cone our characterisation computes.

*Proof.* Suppose some p′ ∈ Q_η had g*·p′ > ω*. Then F(p′) ≤ P(g*, g*·p′) < P(g*, ω*) = F(p*), which contradicts the minimality of F at p*. ∎

Assumptions and caveats:
- (i) The proposition needs the global optimum over honest devices. We have only numerical optima over real qubit x–z devices, and we found several branches (η = 0.98). The optimum column is therefore a lower bound on the true best, not a proved ceiling.
- (ii) Strict monotonicity of the score-to-guessing curve is assumed. It holds numerically at every point.
- (iii) The cone itself, N_Q = C(p), is proved only where Thm 2 or Thm A apply. Elsewhere it is numerical, to ~1e-10.

## Exact vs numerical
- **Proved:**
  - every CHSH-protocol upper bound;
  - every certified entry in the table (each is a rigorous lower bound on H_min for that explicit rational witness and honest device, two-sided tight);
  - the proposition above, under its stated hypotheses;
  - N_Q ⊆ C(p).
- **Numerical:**
  - the device optimisation, so the "optimum" column is not a proven maximum;
  - the L3 values at η = 0.95 and 0.90;
  - N_Q = C(p) at the WAP points and at the optimal devices;
  - the agreement wappush = full.

## Suggested headline (2–3 sentences)
"With symmetric lossy detectors, we compute the optimal per-round Bell witness for device-independent randomness: it lies in the normal cone of its own observed behaviour, which our characterisation reduces to a small polygon. The resulting explicit witnesses certify 9–18% more min-entropy per round than any CHSH protocol over η = 0.75–0.98 (+16–18% at η = 0.85–0.95), with exact-rational certificates. That is a further 3–6% over the better of the asymmetric-CHSH (Woodhead–Acín–Pironio) and polygon witnesses at every η, each protocol on the same detectors with its own best state and settings."
