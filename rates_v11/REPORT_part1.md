# Certified device-independent randomness rates for the tilted-CHSH family (rates_v11)

*Status: one-sided results (Sec. 1–6) final; global two-party results in Sec. 7. All numbers below are produced by the
scripts in this directory; raw records are in `out/` (JSONL, one record per point).*

## 0. Summary

| question | answer | rigour |
|---|---|---|
| Same partially entangled state: tilted witness vs CHSH score | tilted ≫ CHSH (e.g. 27.76°, v = 0.99: 0.838 vs 0.537 bits asymptotic; EAT n = 10¹⁰: 0.790 vs 0.522) | certified (BFF dual + exact-rational residual; CHSH curve analytic) |
| Best σ-symmetric witness in the exact cone vs tilted CHSH | no gain: the search returns tilted CHSH (π/8 all v; 27.76° all v; 3π/16 v = 0.99) | each evaluated bound certified; the *search* over the cone is heuristic (not a proof of optimality within the cone) |
| White-noise envelope over θ, one-sided H(A₁\|X=1,E) | won by θ = π/4 (max-entangled CHSH) at **every** v ∈ [0.86, 0.999], for von Neumann (BFF) and H_min (NPA) | vN: certified lower bounds for θ < π/4 vs the analytic (tight) CHSH curve at π/4, so the comparison is conservative only for the tilted side — see Sec. 6 |
| Global two-party randomness | see Sec. 7 | |

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

