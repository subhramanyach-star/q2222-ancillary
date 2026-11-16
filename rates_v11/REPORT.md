# Certified device-independent randomness rates for the tilted-CHSH family (rates_v11)

*Status: final. One-sided results are in Sec. 1–6 and global two-party results in Sec. 7. All numbers below are produced by the
scripts in this directory; raw records are in `out/` (JSONL, one record per point).*

## 0. Summary

| question | answer | rigour |
|---|---|---|
| Same partially entangled state: tilted witness vs CHSH score | tilted ≫ CHSH (e.g. 27.76°, v = 0.99: 0.838 vs 0.537 bits asymptotic; EAT n = 10¹⁰: 0.790 vs 0.522) | certified (BFF dual + exact-rational residual; CHSH curve analytic) |
| Best σ-symmetric witness in the exact cone vs tilted CHSH | no gain: the search returns tilted CHSH (π/8 all v; 27.76° all v; 3π/16 v = 0.99) | each evaluated bound certified; the *search* over the cone is heuristic (not a proof of optimality within the cone) |
| White-noise envelope over θ, one-sided H(A₁\|X=1,E) | won by θ = π/4 (max-entangled CHSH) at **every** v ∈ [0.86, 0.999], for von Neumann (BFF) and H_min (NPA) | vN: certified lower bounds for θ < π/4 vs the analytic (tight) CHSH curve at π/4, so the comparison is conservative only for the tilted side — see Sec. 6 |
| Global two-party randomness H(AB\|X=x*,Y=2,E) | max-ent CHSH wins at every v ∈ [0.95, 0.999]; best tilted is 0.22π (v = 0.99: 1.753 vs 1.815 bits) | both sides certified lower bounds; strict loss rigorously proved (device-dependent upper bound) for θ ≤ 27.76° at v ≥ 0.98, and for 3π/16 and π/16 at v ≥ 0.994 |

**Bottom line for the paper (one-sided):** the large "gain over CHSH" is a *same-state* statement (witness matched to the
state vs the CHSH score evaluated on that state). It is not an advantage over the best CHSH experiment: under white noise a
maximally entangled CHSH experiment certifies more one-sided randomness per round than any tilted target with any witness
(including full statistics) at every visibility we tested, asymptotically and at finite size.

## 1. Quantities and protocol

* Scenario 2222, correlators z = (⟨A₀⟩,⟨A₁⟩,⟨B₀⟩,⟨B₁⟩,⟨A₀B₀⟩,⟨A₀B₁⟩,⟨A₁B₀⟩,⟨A₁B₁⟩).
* Target p_θ: |ψ_θ⟩ = cos θ|00⟩ + sin θ|11⟩, A₀ = σ_z, A₁ = σ_x, B_{0,1} = cos μ σ_z ± sin μ σ_x, cos μ = 1/r,
  r = √(1+sin²2θ) (the tilted-CHSH-optimal strategy; θ = π/4 is the CHSH Tsirelson point).
* Noise: white noise ρ → vρ + (1−v)𝟙/4, i.e. z → v·z.
* Randomness: one-sided H(A₁|X=1,E) (A₁ = σ_x has a uniform marginal for every θ, so the ideal value is 1 bit for all θ).
* Constraints ("protocols"): (i) **CHSH score** S ≥ 2rv on the tilted state; (ii) **tilted-CHSH score** g_θ·z ≥ v g_θ·p_θ;
  (iii) **best exposing witness**: σ-symmetric convex combinations of the exact hexagon vertices of the normal cone C_θ
  (every such g is maximised by p_θ), searched on the simplex; (iv) **full statistics** z = v p_θ (reference upper envelope).

## 2. Method

**Asymptotic rate (von Neumann).** Brown–Fawzi–Fawzi (BFF), arXiv:2106.13692, Lemma 2.3: Gauss–Radau quadrature with m
nodes (t_m = 1 dropped, which only lowers the bound), one non-commutative polynomial SDP per node ("speed-up 3": separate
infimum per node), operator-norm bound ‖Z_a‖ ≤ α_i = (3/2)max(1/t_i, 1/(1−t_i)). Each node is relaxed to a real NPA moment
matrix (level 2 plus {A,B,AB}×{Z,Z*}); the real part of a complex moment matrix is feasible, so the relaxation is valid.
m = 8 everywhere, m = 12 additionally for v ≥ 0.97 (best of the two reported).

**Solver and certification (`bffcert.py`).** SDPA (sdpa-python) with several trials per node (Z-scalings α^{1, 0.5, 0.75},
alternative SDPA parameters, adaptive row scaling on the two hardest small-t nodes). The solver's dual PSD matrix S is
turned into a **rigorous** bound: S is PSD-repaired in the scaled frame (negative eigenvalues clipped, a floating-point
eigenvalue margin subtracted), the residual r = c − ⟨S,F_k⟩ − Aᵀμ is computed in **exact rational arithmetic**, and every
non-zero residual entry is charged against the a-priori bound |⟨w⟩| ≤ α^{#Z(w)}. The result is an affine function
K₀ + κ·z valid for **all** quantum strategies — a lower bound on H at the observed point and, at the same time, an EAT
min-tradeoff function. The solver's primal value is never used as a bound. Clarabel stalls on these SDPs; MOSEK is not
available; COPT is size-limited.

**H_min (cross-check).** NPA level-2 guessing-probability SDPs (Clarabel), numerical only (not dual-certified).

**Finite size (`eat.py`).** Liu et al., Nat. Phys. 17, 448 (2021) [arXiv:1912.11159], App. A Thm 3 + spot-checking Lemma 3,
generalised from the CHSH game to an arbitrary affine min-tradeoff function. Conservative simplifications: Δ(f,p) ≥ 0
dropped; the variance term replaced by its supremum over all no-signalling test distributions. ε_h = 1.545·10⁻¹²,
ε_EAT = 3.09·10⁻¹², completeness ε_c = 10⁻⁶ (Bernstein), test probability γ, Rényi α and the affine offset optimised per n.
Min-tradeoff functions are the certified (K₀, κ) from BFF at design visibilities within 0.03 below the operating point.
The CHSH curve uses tangents of the analytic bound.

**Analytic CHSH cross-check** (`check_chsh.py`, `out/table_chsh_check.md`). Against the tight bound
H ≥ 1 − h(½ + ½√(S²/4 − 1)) (Pironio et al. NJP 11, 045021 (2009)): certified ≤ analytic at every S; deficit
10⁻⁴–2·10⁻³ for S ≤ 2.8, 7.7·10⁻³ at S = 2.8256 for m = 8 (3.7·10⁻³ for m = 12). The deficit is the Radau truncation
(t = 1 node dropped, ≈ 1/(2m² ln 2)) plus relaxation slack, so the pipeline is conservative and tight to ~10⁻³ bits.

## 3. Same-state gains (asymptotic, certified; `out/table_asym.md`)

| θ | v | CHSH score (analytic) | tilted CHSH | best in cone | full stats | max-ent CHSH (π/4) |
|---|---|---|---|---|---|---|
| π/8 | 0.999 | 0.396 | 0.964 | 0.964 | 0.965 | 0.989 |
| π/8 | 0.99 | 0.373 | 0.699 | 0.699 | 0.762 | 0.919 |
| π/8 | 0.98 | 0.347 | 0.408 | 0.408 | 0.616 | 0.857 |
| π/8 | 0.97 | 0.321 | 0.118 | 0.118 | 0.513 | 0.803 |
| 27.76° | 0.999 | 0.567 | 0.978 | 0.978 | 0.978 | 0.989 |
| 27.76° | 0.99 | 0.537 | 0.839 | 0.839 | 0.855 | 0.919 |
| 27.76° | 0.98 | 0.505 | 0.692 | 0.692 | 0.750 | 0.857 |
| 27.76° | 0.95 | 0.413 | 0.266 | 0.266 | 0.538 | 0.708 |
| 3π/16 | 0.99 | 0.722 | 0.895 | 0.896 | 0.901 | 0.919 |

* The tilted witness dies at v = h_L(g_θ)/h_Q(g_θ) = (cos 2θ + r)/2 (0.966 at π/8, 0.931 at 27.76°); the CHSH score on
  the same state dies at 1/r (0.816, 0.772). Below v ≈ 0.977 (π/8) / 0.963 (27.76°) the CHSH score wins on the
  same state; the full-statistics bound is above both everywhere.
* **Cone search** (`run_cone.py`, `out/bff_cone.jsonl`): over σ-symmetric combinations of the hexagon vertices, the best
  certified witness coincides with tilted CHSH to within 10⁻³ bits at every (θ, v) tested — no gain over tilted CHSH.
  The full-statistics bound is the ceiling for *any* witness; the gap tilted→full (e.g. 0.70 → 0.76 at π/8, v = 0.99) is
  only reachable with witnesses outside the σ-symmetric family / outside the ideal cone.

## 4. Finite-size (EAT) rates (`out/table_eat.md`, `out/eat_table.json`)

H_min^ε/n per round (bits):

| θ | v | n | CHSH score | tilted (= best in cone) | full-stat. dual expr. | max-ent CHSH (π/4) |
|---|---|---|---|---|---|---|
| 27.76° | 0.994 | 10⁸ | 0.480 | 0.684 | 0.720 | 0.825 |
| 27.76° | 0.994 | 10¹⁰ | 0.535 | 0.848 | 0.860 | 0.919 |
| 27.76° | 0.99 | 10⁸ | 0.468 | 0.631 | 0.684 | 0.803 |
| 27.76° | 0.99 | 10¹⁰ | 0.522 | 0.790 | 0.814 | 0.893 |
| π/8 | 0.994 | 10¹⁰ | 0.371 | 0.737 | 0.777 | 0.919 |
| π/8 | 0.99 | 10¹⁰ | 0.361 | 0.625 | 0.708 | 0.893 |
| π/8 | 0.999 | 10¹⁰ | 0.384 | 0.878 | 0.882 | 0.955 |

Same-state gains survive finite size (27.76°, v = 0.99, n = 10¹⁰: +51% tilted/CHSH, +56% full/CHSH), but the
maximally entangled CHSH experiment is higher at every n and v in the table.

## 5. White-noise envelope over θ (`out/table_envelope.md`, `out/table_theta_detail.md`)

θ ∈ {π/16, 0.1π, π/8, 27.76°, 3π/16, 0.22π, π/4}, v ∈ {0.86 … 0.999}. At every v, the best θ for the tilted witness,
for full statistics and for H_min is θ = π/4. Examples (vN, certified): v = 0.99 — π/4 0.919, 0.22π 0.915, 3π/16 0.895,
27.76° 0.839, π/8 0.699; v = 0.95 — π/4 0.708, 0.22π 0.687, 3π/16 0.574, 27.76° 0.266, π/8 0. At v = 0.999 the
next best (0.22π, tilted 0.985) is within the m = 8 truncation error of π/4 (0.989), so "π/4 wins" is only a statement
of no advantage there, not of a strict loss. H_min (NPA, numerical) gives the same ordering at every v.

## 6. Why maximally entangled CHSH wins the one-sided white-noise envelope

1. **No room for a tilted advantage in this figure of merit.** A₁ = σ_x has a uniform marginal on every |ψ_θ⟩ and the
   target self-tests, so the ideal one-sided randomness is exactly 1 bit for all θ. The known benefits of partial
   entanglement — up to 2 bits of *global* randomness from weakly entangled states (Acín–Massar–Pironio 2012) and higher
   tolerance to *detection loss* (Eberhard 1993) — are not rewarded by H(A₁|X=1,E) under white noise. The only thing left
   to compare is robustness to white noise, and that is a contest the maximally entangled state wins.
2. **Nonlocal margin.** White noise scales z → v z, i.e. it moves every point straight toward the maximally mixed
   (origin) behaviour. For a witness g normalised to g·p = 1, certification needs v > h_L(g). For the exposing witnesses
   h_L(g_θ) = (cos 2θ + r)/2: 0.707 (π/4), 0.795 (0.22π), 0.872 (3π/16), 0.931 (27.76°), 0.966 (π/8), 0.997 (π/16).
   Less entangled targets sit closer to a face of the local polytope, so a given amount of noise consumes a larger
   fraction of their margin: at v = 0.99, (1−v)/(1−h_L) = 3.4% for π/4 but 14.5% for 27.76° and 29% for π/8.
   Even full statistics cannot recover this — the noisy partially entangled state is simply less nonlocal than the noisy
   singlet at the same visibility (CHSH-violation thresholds 1/r: 0.707 at π/4 vs 0.816 at π/8).
3. **Steep self-testing slope.** Eve's best attack mixes the quantum strategy with deterministic (fully predictable)
   strategies. For the tilted witness, the deterministic strategies that saturate h_L are only (1 − h_L) below the
   quantum value in score, so they cost Eve very little; the certified entropy therefore falls steeply from 1 bit
   (π/8: 0.964 → 0.699 between v = 0.999 and 0.99, slope ≈ 29 bits per unit v) while the CHSH curve is shallow
   (0.989 → 0.919, slope ≈ 7.7).
4. **What the same-state gain does mean.** On a *fixed* partially entangled state the CHSH score is a badly matched
   witness (S = 2rv instead of 2√2), so matching the witness to the state is worth a lot (Sec. 3). That is a statement about
   choosing the Bell expression given the source, not about which source to build.
5. **Scope of the conclusion.** White noise is the natural model for imperfect visibility, but it is not the dominant
   imperfection in photonic DI-QRNG experiments, which are loss-limited; there partially entangled states are used
   precisely because of Eberhard's detection-efficiency advantage (e.g. Liu et al. 2021). The envelope result therefore does
   not transfer to loss-dominated experiments; that comparison needs a detection-efficiency noise model (not done here).

## 7. Global (two-party) randomness H(AB|X=x*, Y=y*, E)

**Setup** (`bffcert2.py`, `targets2.py`, `run_global.py`). This is a 2×3 scenario in which Bob gets an extra key setting B₂, as in Acín–Massar–Pironio 2012 and Brown–Fawzi–Fawzi. Alice keeps the tilted-CHSH settings A₀ = σ_z, A₁ = σ_x, and Bob keeps B₀ and B₁. B₂ = cos φ σ_z + sin φ σ_x, and the pair (x*, φ) is chosen to maximise the ideal Shannon entropy of p(ab|x*,2). For θ = π/4 this is the maximally entangled CHSH (Tsirelson) protocol with A₀ = σ_z, B₂ = σ_x, giving an ideal 2 bits. For tilted θ, x* = 0 and φ = π/2, except at 0.1π and π/16 where x* = 1 and φ = π/2.

White noise is applied as z → v z.

There are two protocols:
- **witness**: the exposing witness score ≥ v·g·p on the 2222 block, plus the three B₂-row correlators ⟨B₂⟩, ⟨A₀B₂⟩, ⟨A₁B₂⟩. The witness is tilted CHSH for θ < π/4 and CHSH at π/4.
- **full**: all 11 correlators.

The method is BFF with m = 8 and projector letters. Each node is a real NPA relaxation with monomials {1,A,B,AB}×{Z_ab, Z_ab*} for the 4 outcomes, giving a 116×116 moment matrix. The SDPs are solved with SDPA and certified with the same exact-rational dual certification as Sec. 2. Ideal global randomness is below 2 bits for every θ < π/4 when all measurements are real qubit measurements. With uniform outcomes this would need σ_x⊗σ_y-type complex settings, and certifying those needs more settings; that was not explored.

### 7.1 Certified H(AB|X=x*,Y=2,E), bits per round (witness / full) — `out/global/table_global.md`

| v | pi/4 witness / full | pi/8 witness / full | 27.76deg witness / full | 3pi/16 witness / full | 0.22pi witness / full | 0.1pi witness / full | pi/16 witness / full |
|---|---|---|---|---|---|---|---|
| 0.999 | 1.9632 / 1.9632 | 1.3718 / 1.4139 | 1.6352 / 1.6672 | 1.8204 / 1.8400 | 1.9257 / 1.9268 | 1.2825 / 1.3011 | 0.5868 / 0.7089 |
| 0.994 | 1.8787 / 1.8790 | 0.7741 / 0.9074 | 1.3410 / 1.3922 | 1.6676 / 1.6897 | 1.8265 / 1.8348 | 0.5771 / 0.7031 | 0.0000 / 0.2683 |
| 0.990 | 1.8149 / 1.8152 | 0.4778 / 0.7079 | 1.1510 / 1.2311 | 1.5556 / 1.5860 | 1.7529 / 1.7647 | 0.3038 / 0.5506 | 0.0000 / 0.1720 |
| 0.980 | 1.6667 / 1.6672 | 0.1408 / 0.4888 | 0.7497 / 0.9552 | 1.3174 / 1.3737 | 1.5858 / 1.5999 | 0.0000 / 0.3581 | 0.0000 / 0.0920 |
| 0.970 | 1.5309 / 1.5316 | 0.0282 / 0.3943 | 0.4571 / 0.7779 | 1.1009 / 1.2040 | 1.4332 / 1.4559 | 0.0000 / 0.2672 | 0.0000 / 0.0547 |
| 0.950 | 1.2832 / 1.2844 | 0.0000 / 0.2841 | 0.1192 / 0.5616 | 0.7250 / 0.9405 | 1.1443 / 1.1993 | 0.0000 / 0.1709 | 0.0000 / 0.0146 |

### 7.2 Envelope: best tilted θ < π/4 vs maximally entangled CHSH


| v | max-ent CHSH witness+B2 | max-ent full | best tilted witness (theta) | best tilted full (theta) | Hmin max-ent witness | Hmin best tilted witness (theta) | Hmin max-ent full | Hmin best tilted full (theta) |
|---|---|---|---|---|---|---|---|---|
| 0.999 | 1.9632 | 1.9632 | 1.9257 (0.22pi) | 1.9268 (0.22pi) | 1.7756 | 1.4973 (0.22pi) | 1.7756 | 1.5452 (0.22pi) |
| 0.994 | 1.8787 | 1.8790 | 1.8265 (0.22pi) | 1.8348 (0.22pi) | 1.4876 | 1.2425 (0.22pi) | 1.4877 | 1.2907 (0.22pi) |
| 0.990 | 1.8149 | 1.8152 | 1.7529 (0.22pi) | 1.7647 (0.22pi) | 1.3596 | 1.1234 (0.22pi) | 1.3597 | 1.1832 (0.22pi) |
| 0.980 | 1.6667 | 1.6672 | 1.5858 (0.22pi) | 1.5999 (0.22pi) | 1.1442 | 0.9167 (0.22pi) | 1.1444 | 1.0114 (0.22pi) |
| 0.970 | 1.5309 | 1.5316 | 1.4332 (0.22pi) | 1.4559 (0.22pi) | 0.9942 | 0.7700 (0.22pi) | 0.9944 | 0.8886 (0.22pi) |
| 0.950 | 1.2832 | 1.2844 | 1.1443 (0.22pi) | 1.1993 (0.22pi) | 0.7802 | 0.5602 (0.22pi) | 0.7803 | 0.7045 (0.22pi) |

(The H_min columns are NPA level-2 guessing-probability SDPs from `npa_hmin2.py` / `run_hmin2.py`. They are numerical and not dual-certified.)

### 7.3 Rigorous strict separation (`upper_global.py`, `out/global/table_separation.md`)

The certified numbers above are lower bounds on both sides. To show that max-ent CHSH is *strictly* better, a rigorous upper bound on the tilted side is needed. One admissible attack is for Eve to hold a purification of ρ_v = vψ_θ + (1−v)𝟙/4. For rank-1 projective measurements this gives H(AB|E) = H_Shannon(p_ab) − S(ρ_v), which is an upper bound on the DI quantity. We maximise it over all real (xz-plane) measurement pairs for Alice and Bob and call it U_real. "<" below means U_real < the certified max-ent lower bound, i.e. a rigorous strict loss for that source with any real key pair.

| theta | v=0.999 | v=0.994 | v=0.990 | v=0.980 | v=0.970 | v=0.950 |
|---|---|---|---|---|---|---|
| 0.22pi | 1.965 ≥ | 1.926 ≥ | 1.900 ≥ | 1.839 ≥ | 1.785 ≥ | 1.687 ≥ |
| 3pi/16 | 1.883 < | 1.846 < | 1.820 ≥ | 1.761 ≥ | 1.709 ≥ | 1.614 ≥ |
| 27.76deg | 1.755 < | 1.719 < | 1.694 < | 1.639 < | 1.589 ≥ | 1.500 ≥ |
| pi/8 | 1.625 < | 1.592 < | 1.569 < | 1.519 < | 1.474 < | 1.394 ≥ |
| 0.1pi | 1.724 < | 1.689 < | 1.664 < | 1.609 < | 1.560 ≥ | 1.471 ≥ |
| pi/16 | 1.882 < | 1.844 < | 1.818 ≥ | 1.760 ≥ | 1.707 ≥ | 1.612 ≥ |
| max-ent CHSH certified L | 1.963 | 1.879 | 1.815 | 1.667 | 1.531 | 1.283 |

### 7.4 Answer for global randomness

* **No tilted target beats the maximally entangled CHSH protocol at any tested visibility v ∈ {0.95, 0.97, 0.98, 0.99, 0.994, 0.999}.** This holds with its exposing witness and also with full 2×3 statistics. The best tilted θ is always the least tilted one tested, 0.22π:

  | v | 0.22π, witness | max-ent CHSH, witness | gap |
  |---|---|---|---|
  | 0.999 | 1.926 | 1.963 | −0.037 |
  | 0.99 | 1.753 | 1.815 | −0.062 |
  | 0.98 | 1.586 | 1.667 | −0.081 |
  | 0.95 | 1.144 | 1.283 | −0.139 |

  The gap grows as v falls, and the ordering is monotone in θ: more tilt gives less certified global randomness. The NPA H_min check gives the same ordering, with larger gaps (v = 0.99: 1.360 vs 1.123).
* **Rigour.** The certified lower bounds for both sides are rigorous, up to the stated BFF truncation and relaxation, which are conservative. The *strict* statement "tilted is worse" is proved rigorously only where U_real < L_max-ent:
  - θ ≤ 27.76° (π/8, 0.1π, 27.76°): for v ≥ 0.98, and for π/8 down to 0.97.
  - 3π/16 and π/16: for v ≥ 0.994.

  For 0.22π, and at lower v, the device-dependent upper bound is too loose to separate. There the conclusion rests on comparing certified lower bounds, which is numerically clear but not a proof.
* **Why.** At v = 1, max-ent CHSH already reaches the 2-bit ceiling for two binary outputs. A tilted source with real measurements cannot have a uniform 4-outcome distribution: uniform outcomes need ⟨A⟩ = ⟨B⟩ = ⟨AB⟩ = 0, which forces A = B = σ_x and then ⟨AB⟩ = sin 2θ ≠ 0. So its ideal value is strictly below 2 bits: 1.60 at π/8, 1.75 at 27.76°, 1.97 at 0.22π. It also loses robustness faster, for the same nonlocal-margin reason as Sec. 6. The well-known AMP advantage of partial entanglement is about *randomness per unit of entanglement*, approaching 2 bits as θ → 0 with vanishing noise tolerance. It is not about randomness per round at fixed white-noise visibility.
* **Not covered.** This study does not cover complex-measurement constructions with extra settings for both parties, more than two outputs, or noise models dominated by detection loss. Detection loss is where partially entangled states are known to help (Eberhard).

## 8. Files

* Code:
  - One-sided: `bffcert.py`, `run_main.py`, `run_env.py`, `run_cone.py`, `run_hmin.py`, `npa_hmin.py`, `eat.py`, `analyze.py`, `envelope.py`, `check_chsh.py`, `make_fig.py`.
  - Global: `bffcert2.py`, `targets2.py`, `run_global.py`, `run_global_cfg.py`, `npa_hmin2.py`, `run_hmin2.py`, `upper_global.py`, `analyze_global.py`.
* One-sided data:
  - Raw: `out/bff_main.jsonl`, `out/bff_cone.jsonl`, `out/hmin_env.json`, `out/eat_table.json`, `out/chsh_check.json`.
  - Tables: `out/table_asym.md`, `out/table_eat.md`, `out/table_envelope.md`, `out/table_theta_detail.md`, `out/table_chsh_check.md`.
* Global data:
  - Raw: `out/global/bff_global.jsonl`, `out/global/hmin_global.jsonl`, `out/global/upper_global.json`, `out/global/envelope_global.json`.
  - Tables: `out/global/table_global.md`, `out/global/table_separation.md`.
* Figures:
  - `fig_rates_v11.pdf/.png`: 2×2 one-sided figure: (a) 27.76°, (b) π/8, (c) θ-envelope, (d) EAT vs max-ent.
  - `fig_global_v11.pdf/.png`: global figure, witness and full statistics.
