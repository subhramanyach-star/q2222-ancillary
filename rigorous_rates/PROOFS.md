# Rigorous min-entropy statements for Results 2–3, Fig. 2 and Table `tab:rates` of prl_v13

*Ravishankar Ramanathan — working notes, `/workspace/ndwu/rigorous_rates/` (22 Oct 2026). Nothing here edits the prl files.*

## 0. What the paper's numbers are (identified from the code that produced them)

`anc_submission_v7/rates_numerics/{rates.py,npa.py,dense.py}` (the source of every number in Results 2–3, Fig. 2 and
Table `tab:rates`) computes

* **Tilted targets** (θ = π/8, 27.76°, 3π/16; Fig. 2(a), Table, Liu-tilt numbers):
  **H_min(A₁|X=1,E) = −log₂ P_g**, the guessing probability of Alice's output **A₁** (a *local*, one-output guess),
  under either one single-score constraint **g·z ≥ v (g·p_θ)** (white noise, z ↦ v z) or the **full statistics z = v p_θ**.
* **WBC panel** (Fig. 2(b), I_φ with φ = π/6): **H_min(A₀B₀|X=0,Y=0,E)** (a *global* guess of the pair).
  H_min(A₀|E) and φ = π/12 are in the ancillary json.
* **P_T paragraph** (Appendix "Searches"): H_min(A₀|E) and H_min(A₀B₀|00,E) at the Tsirelson point.
* All values were floating-point NPA level-3 (Clarabel) optima; "best in cone" was a 120-sample Dirichlet search +
  coordinate refinement over the six hexagon vertices, i.e. only a lower bound on the supremum over the hexagon.

Coordinates: z = (⟨A₀⟩,⟨A₁⟩,⟨B₀⟩,⟨B₁⟩,⟨A₀B₀⟩,⟨A₀B₁⟩,⟨A₁B₀⟩,⟨A₁B₁⟩); ±1-valued observables.
Q = the set of quantum behaviours (commuting-operator model; all upper bounds below hold in that generality, and
all lower bounds are realised by explicit two-qubit strategies, so the results hold for the tensor-product model too).

## 1. The guessing problems and two duality lemmas

Eve's guess e ranges over the outputs being guessed (e = a for A_x; e = (a,b) for (A_x,B_y)). Let Π_e be the
corresponding affine functional of z: Π_a = (1 + a⟨A_x⟩)/2, Π_{ab} = (1 + a⟨A_x⟩ + b⟨B_y⟩ + ab⟨A_xB_y⟩)/4.
As is standard [Pironio et al. 2010; Nieto-Silleras et al. 2014; Bancal et al. 2014], conditioning on Eve's
measurement outcome decomposes any attack into sub-normalised quantum behaviours z̃_e = t_e z_e (t_e ≥ 0, Σt_e = 1,
z_e ∈ Q), and

* single score: P_g(g,ω) = sup { Σ_e t_e Π_e(z_e) : Σ_e t_e g·z_e ≥ ω },
* full statistics: P_g(p) = sup { Σ_e t_e Π_e(z_e) : Σ_e t_e z_e = p }.

h_Q(c) := sup_{z∈Q} c(z) for an affine functional c.

**Lemma 1 (single score).** For every μ ≥ 0, P_g(g,ω) ≤ max_e h_Q(Π_e + μg) − μω.
*Proof.* Σ_e t_e Π_e(z_e) ≤ Σ_e t_e[Π_e + μg](z_e) − μω ≤ Σ_e t_e h_Q(Π_e + μg) − μω ≤ max_e h_Q(Π_e+μg) − μω,
using Σ_e t_e g·z_e ≥ ω and Σ t_e = 1. ∎ (Strong duality holds, so the bound is tight for the optimal μ; we never need
this — tightness is shown by explicit attacks, Lemma 5.) The function ω ↦ P_g(g,ω) is non-increasing and concave, and
each certificate (μ, bound on h_Q) gives an affine upper bound valid **for all ω simultaneously**.

**Lemma 2 (full statistics).** For every f ∈ ℝ⁸, P_g(p) ≤ max_e h_Q(Π_e + f) − f·p.
*Proof.* Σ_e t_e Π_e(z_e) = Σ_e t_e [Π_e + f](z_e) − f·Σ_e t_e z_e ≤ max_e h_Q(Π_e + f) − f·p. ∎

So every number reduces to **upper bounds on quantum maxima of 2222 Bell expressions with marginal terms**, which
can be certified exactly:

**Lemma 3 (exact NPA / sum-of-squares certificate).** Let Γ be the real moment matrix of the words of an NPA level
(1+AB, 2 or 3; 9, 13, 25 words) of any quantum realisation; Γ ⪰ 0, Γ₁₁ = 1, and entries related by the algebra
(A_x² = B_y² = 1, [A_x,B_y] = 0, w ~ w† for the real part) are equal. If a rational symmetric matrix Z satisfies
(i) Z ⪰ 0, (ii) Σ_{(i,j)∈class k} Z_ij = −c_k for the eight classes of the moments in z, and = 0 for every other
non-identity class, then c·z ≤ β := tr Z for every quantum behaviour.
*Proof.* ⟨Z,Γ⟩ = Σ_k y_k Σ_{class k} Z_ij = β − c·z, and ⟨Z,Γ⟩ ≥ 0 because both are PSD. (The real part of a complex
PSD moment matrix is real PSD with the same real moments, so a real relaxation is valid.) ∎
We obtain Z from an SDP solver, add a margin ε·𝟙 (ε = 10⁻⁹…10⁻⁶), round to dyadic rationals, restore (ii) **exactly**
by spreading each class residual uniformly over the class, and check (i) **exactly** by rational LDLᵀ with symmetric
pivoting (python-flint `fmpq`). No floating-point number enters the bound.

**Lemma 4 (irrational data).** |h_Q(c) − h_Q(c̃)| ≤ ‖c − c̃‖₁ (all |z_j| ≤ 1). The targets, the hexagon vertices
(solutions of the 4×4 linear systems of `families.py`), g_θ, WBC functionals, thresholds ω = v g·p are evaluated
in 200-bit ball arithmetic (`families_arb.py`, arb), certificates are issued for a dyadic c̃, and the ‖·‖₁ error
(≈10⁻¹⁶) is added rigorously.

**Lemma 5 (explicit attacks ⇒ lower bounds on P_g).** If z_k ∈ Q (explicit qubit strategies or deterministic points),
π_k ≥ 0, Σπ_k = 1 and Σπ_k g·z_k ≥ ω, then P_g(g,ω) ≥ Σ_k π_k max_e Π_e(z_k) (Eve holds the classical label k and
guesses the likeliest output). For full statistics we use Σπ_k z_k ≈ p/(1−t) and write p = (1−t)Σπ_k z_k + t w with
t = 10⁻⁸; Fine's theorem (16 positivity + 8 CHSH inequalities, checked in ball arithmetic) shows w is local, hence a
mixture of deterministic points that Eve guesses perfectly, so P_g(p) ≥ (1−t)Σπ_k max_e Π_e(z_k) + t.
Strategies are stored with **dyadic-rational angles and state amplitudes**; their behaviours are recomputed in ball
arithmetic, so the constraints are verified rigorously (no tolerance).

**Lemma 6 (supremum over the whole normal cone).** Let V₁,…,V₆ be the hexagon vertices (V_i·p_θ = 1, so every
g ∈ C_θ = conv{V_i} has g·p_θ = 1 and threshold ω = v). If an attack satisfies V_i·z̄ ≥ v for all six i, then it
satisfies g·z̄ ≥ v for every g ∈ C_θ, so its guessing probability lower-bounds P_g(g, v) **for every g in the
hexagon at once**: sup_{g∈C_θ} H_min(g) ≤ −log₂ P_attack. Conversely Lemma 1 at one explicit g* ∈ C_θ (rational
convex weights, checked) gives sup ≥ H_min(g*). The two bounds are computed independently.

Explicit attacks are found by column generation (LP over a dictionary of qubit strategies; new columns by maximising
the LP-dual functional over qubit strategies — by Jordan's lemma qubit blocks suffice; one-sided trivial blocks give
local behaviours and are covered by the 16 deterministic points). The LP dual also supplies μ (Lemma 1), the full-
statistics functional f (Lemma 2) and g* (Lemma 6). None of this search needs to be trusted: only the stored
certificates and attacks are checked.

**Exact closed forms used as cross-checks.**
(a) CHSH score S: P_g ≤ ½ + ½√(2 − S²/4), tight [Pironio et al., Nature 464, 1021 (2010)]. On the tilted target,
CHSH(p_θ) = 2r (r = √(1+sin²2θ)), so with white noise **H_min^CHSH(θ,v) = −log₂[½(1 + √(2 − r²v²))]**
(exact). Our certified intervals reproduce it to 10⁻⁸.
(b) **Exact zeros.** If v ≤ h_L(g) (local bound of g normalised to g·p = 1), a deterministic local point attains the
score, so P_g = 1 and H_min = 0 exactly. For g_θ, h_L = (cos2θ + r)/2 = 0.96593 (π/8), 0.93104 (27.76°); every hexagon
vertex has h_L = 1, and by Corollary `cor:window` every in-cone witness certifies 0 for v below the exact window
threshold. These are the zeros in Table `tab:rates`.

## 2. Results (two-sided, rigorous)

Each entry is an interval [H_lo, H_hi] containing the exact single-round min-entropy (in bits):
H_lo from Lemmas 1–4 (exact-rational NPA certificate at the level shown), H_hi from an explicit qubit attack
(Lemma 5). "NPA-3 (paper)" is the floating-point value printed in prl_v13. See `TABLE.md` (generated by
`make_table.py` from `out/verified.json`, i.e. from the standalone verifier, not from the solver run).

Each cell: exact value lies in [H_lo, H_hi]; a single number means both ends agree to 6 decimals (interval width < 10⁻⁶).
H_lo: exact-rational NPA certificate (level in brackets). H_hi: explicit qubit attack. Paper value (NPA-3 float) below each row.

### θ = π/8: H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.776634 | 0.383682 | 0.185275 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.777 | 0.384 | 0.185 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | [0.832439, 0.832440] | 0.438340 | 0.206499 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.832 | 0.438 | 0.206 | 0.000 | 0.000 |
| *(NPA level used)* | 2 | 2 | 2 | 1+AB | 1+AB |
| **sup over the whole hexagon C_θ** | 0.832443 | 0.438352 | 0.206530 | 0 (exact) | 0 (exact) |
| *(NPA level used)* | 2,3 | 2 | 2 | – | – |
| CHSH | 0.226659 | 0.210974 | 0.194352 | 0.148883 | 0.084669 |
| *(paper NPA-3)* | 0.227 | 0.211 | 0.194 | 0.149 | 0.085 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.835087, 0.835088] | 0.514289 | 0.376273 | 0.212696 | 0.103973 |
| *(paper NPA-3)* | 0.835 | 0.514 | 0.376 | 0.213 | 0.104 |
| *(NPA level used)* | 3 | 1+AB | 1+AB | 2 | 2 |

### θ = 27.76° (Liu et al.): H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.842291 | 0.548525 | 0.393120 | 0.118590 | 0 (exact) |
| *(paper NPA-3)* | 0.842 | 0.549 | 0.393 | 0.119 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | 0.871352 | 0.606007 | 0.439178 | 0.126896 | 0 (exact) |
| *(paper NPA-3)* | 0.871 | 0.606 | 0.439 | 0.127 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 2,3 | 1+AB |
| **sup over the whole hexagon C_θ** | 0.871353 | 0.606026 | 0.439191 | 0.126918 | 0 (exact) |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 3 | – |
| CHSH | 0.350085 | 0.326529 | 0.302144 | 0.237970 | 0.152205 |
| *(paper NPA-3)* | 0.350 | 0.327 | 0.302 | 0.238 | 0.152 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | 0.871943 | [0.630677, 0.630678] | 0.505821 | 0.322627 | 0.186098 |
| *(paper NPA-3)* | 0.872 | 0.631 | 0.506 | 0.323 | 0.186 |
| *(NPA level used)* | 1+AB | 2 | 1+AB | 2 | 2 |

### θ = 3π/16 (ancillary): H_min(A₁|X=1,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| tilted CHSH g_θ | 0.884187 | 0.660479 | 0.538610 | 0.318483 | 0.094460 |
| *(paper NPA-3)* | 0.884 | 0.660 | 0.539 | 0.318 | 0.094 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| paper "best in cone" (its weights) | 0.895283 | 0.692396 | 0.578470 | 0.346602 | 0.098757 |
| *(paper NPA-3)* | 0.895 | 0.692 | 0.578 | 0.347 | 0.099 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 1+AB |
| **sup over the whole hexagon C_θ** | 0.895283 | 0.692398 | 0.578473 | 0.346607 | 0.098841 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 1+AB |
| CHSH | 0.527518 | 0.485829 | 0.445285 | 0.347412 | 0.229271 |
| *(paper NPA-3)* | 0.528 | 0.486 | 0.445 | 0.347 | 0.229 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.895491, 0.895492] | 0.697997 | 0.596080 | 0.421675 | 0.267516 |
| *(paper NPA-3)* | 0.895 | 0.698 | 0.596 | 0.422 | 0.268 |
| *(NPA level used)* | 1+AB | 1+AB | 2 | 1+AB | 2 |

### WBC φ = π/6: H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 1.729981 | 1.069929 | 0.588463 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.730 | 1.070 | 0.588 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | [1.693991, 1.693992] | 0.606132 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.694 | 0.606 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB,2 | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.439524 | 0.409375 | 0.378186 | 0.296221 | 0.187197 |
| *(paper NPA-3)* | 0.440 | 0.409 | 0.378 | 0.296 | 0.187 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB,2 |
| full statistics | [1.759378, 1.759379] | [1.163943, 1.163944] | 0.846721 | 0.475418 | 0.248684 |
| *(paper NPA-3)* | 1.759 | 1.164 | 0.847 | 0.475 | 0.249 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### WBC φ = π/6: H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 0.891970 | 0.610278 | 0.359617 | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.892 | 0.610 | 0.360 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 0.884285 | 0.372425 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.884 | 0.372 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,3 | 1+AB,3 | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.356581 | 0.332529 | 0.307668 | 0.242401 | 0.155460 |
| *(paper NPA-3)* | 0.357 | 0.333 | 0.308 | 0.242 | 0.155 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | 0.892805 | 0.669448 | 0.541128 | 0.351942 | 0.200250 |
| *(paper NPA-3)* | 0.893 | 0.669 | 0.541 | 0.352 | 0.200 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 2 | 2 |

### WBC φ = π/12 (ancillary): H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 1.692537 | 0.644473 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.693 | 0.644 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 1.573891 | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 1.574 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.260635 | 0.241782 | 0.221783 | 0.167054 | 0.090246 |
| *(paper NPA-3)* | 0.261 | 0.242 | 0.222 | 0.167 | 0.090 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB,2 | 1+AB,2 |
| full statistics | [1.707286, 1.707288] | 0.902843 | 0.623863 | 0.303884 | 0.113527 |
| *(paper NPA-3)* | 1.707 | 0.903 | 0.624 | 0.304 | 0.114 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### WBC φ = π/12 (ancillary): H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| I_φ (centre) | 0.900733 | 0.400217 | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.901 | 0.400 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to W₁ | 0.848375 | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.848 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB,2 | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex W₁ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| CHSH | 0.214065 | 0.199042 | 0.183093 | 0.139324 | 0.077203 |
| *(paper NPA-3)* | 0.214 | 0.199 | 0.183 | 0.139 | 0.077 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.902934, 0.902935] | 0.598917 | 0.450536 | 0.237619 | 0.097246 |
| *(paper NPA-3)* | 0.903 | 0.599 | 0.451 | 0.238 | 0.097 |
| *(NPA level used)* | 1+AB | 2 | 3 | 3 | 3 |

### Tsirelson point P_T: H_min(A₀|X=0,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| CHSH (centre) | 0.911547 | 0.737568 | 0.642251 | 0.472341 | 0.307179 |
| *(paper NPA-3)* | 0.912 | 0.738 | 0.642 | 0.472 | 0.307 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to g₀ | 0.908267 | 0.726020 | 0.624817 | 0.437441 | 0.194229 |
| *(paper NPA-3)* | 0.908 | 0.726 | 0.625 | 0.437 | 0.194 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex g₀ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [0.911632, 0.911633] | 0.739922 | 0.648357 | 0.492569 | 0.336009 |
| *(paper NPA-3)* | 0.912 | 0.740 | 0.648 | 0.493 | 0.336 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 2 |

### Tsirelson point P_T: H_min(A₀B₀|00,E)

| witness | v = 0.999 | v = 0.99 | v = 0.98 | v = 0.95 | v = 0.9 |
|---|---|---|---|---|---|
| CHSH (centre) | 1.122306 | 0.911258 | 0.794408 | 0.584108 | 0.377572 |
| *(paper NPA-3)* | 1.122 | 0.911 | 0.794 | 0.584 | 0.378 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| half-way to g₀ | 1.118612 | 0.899307 | 0.777202 | 0.544032 | 0.236652 |
| *(paper NPA-3)* | 1.119 | 0.899 | 0.777 | 0.544 | 0.237 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| vertex g₀ | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) | 0 (exact) |
| *(paper NPA-3)* | 0.000 | 0.000 | 0.000 | 0.000 | 0.000 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |
| full statistics | [1.133620, 1.133621] | 0.954577 | 0.857586 | 0.646259 | 0.422413 |
| *(paper NPA-3)* | 1.134 | 0.955 | 0.858 | 0.646 | 0.422 |
| *(NPA level used)* | 1+AB | 1+AB | 1+AB | 1+AB | 1+AB |

### Fig. 2 dense curves (certified points)

| v | (a) tilted | (a) CHSH | (a) full | (b) I_φ centre | (b) CHSH | (b) full |
|---|---|---|---|---|---|---|
| 0.999 | 0.776634 | 0.226659 | [0.835087, 0.835088] | 1.729981 | 0.439524 | [1.759378, 1.759379] |
| 0.995 | 0.540022 | 0.219598 | 0.640180 | 1.388400 | 0.425862 | [1.448192, 1.448193] |
| 0.99 | 0.383682 | 0.210974 | 0.514289 | 1.069929 | 0.409375 | [1.163943, 1.163944] |
| 0.985 | 0.273148 | [0.202561, 0.202562] | 0.433678 | 0.809204 | 0.393498 | [0.977841, 0.977842] |
| 0.98 | 0.185275 | 0.194352 | 0.376273 | 0.588463 | 0.378186 | 0.846721 |
| 0.97 | 0.047068 | 0.178506 | 0.298916 | 0.228099 | 0.349099 | 0.671741 |
| 0.96 | 0 (exact) | 0.163369 | 0.248581 | 0 (exact) | 0.321849 | 0.557727 |
| 0.95 | 0 (exact) | 0.148883 | 0.212696 | 0 (exact) | 0.296221 | 0.475418 |
| 0.94 | 0 (exact) | 0.134999 | 0.185295 | 0 (exact) | 0.272038 | [0.411815, 0.411816] |
| 0.92 | 0 (exact) | 0.108865 | 0.140783 | 0 (exact) | 0.227455 | 0.317337 |
| 0.9 | 0 (exact) | 0.084669 | 0.103973 | 0 (exact) | 0.187197 | 0.248684 |
| 0.88 | 0 (exact) | 0.062172 | 0.072777 | 0 (exact) | 0.150614 | 0.195376 |
| 0.86 | 0 (exact) | 0.041183 | 0.045999 | 0 (exact) | 0.117230 | [0.149035, 0.149036] |
| 0.84 | 0 (exact) | 0.021540 | 0.022925 | 0 (exact) | 0.086696 | 0.107753 |

### von Neumann H(A₁|X=1,E): BFF certified lower bound (../rates_v11) and explicit-attack upper bound (this work)

| θ | v | tilted [BFF, attack] | sup over hexagon [BFF best found, attack for all g ∈ C_θ] | full [BFF, attack] |
|---|---|---|---|---|
| pi/8 | 0.999 | [0.9643, 0.9707] | [0.9643, 0.9707] | [0.9653, 0.9725] |
| pi/8 | 0.99 | [0.6988, 0.7065] | [0.6985, 0.7065] | [0.7619, 0.7729] |
| pi/8 | 0.98 | [0.4075, 0.4130] | [0.4075, 0.4130] | [0.6155, 0.6295] |
| pi/8 | 0.95 | [0.0000, 0.0000] | [–, 0.0000] | [0.3784, 0.3960] |
| 27.76deg | 0.999 | [0.9777, 0.9855] | [0.9775, 0.9855] | [0.9765, 0.9858] |
| 27.76deg | 0.99 | [0.8385, 0.8550] | [0.8385, 0.8550] | [0.8551, 0.8729] |
| 27.76deg | 0.98 | [0.6923, 0.7100] | [0.6922, 0.7100] | [0.7502, 0.7729] |
| 27.76deg | 0.95 | [0.2664, 0.2750] | [0.2660, 0.2750] | [0.5384, 0.5611] |


### Derived numbers quoted in the text

| quantity | paper | rigorous interval |
|---|---|---|
| pi/8 v=0.999: sup_hexagon / CHSH | 3.7 | [3.6727, 3.6727] |
| pi/8 v=0.999: paper-best / CHSH | 3.7 | [3.6727, 3.6727] |
| pi/8 v=0.99: sup_hexagon / CHSH | 2.1 | [2.0778, 2.0778] |
| pi/8 v=0.99: paper-best / CHSH | 2.1 | [2.0777, 2.0777] |
| pi/8 v=0.98: sup_hexagon / CHSH | 1.06 ("6%") | [1.0627, 1.0627] |
| pi/8 v=0.98: paper-best / CHSH | 1.06 ("6%") | [1.0625, 1.0625] |
| pi/8 v=0.999: sup_hexagon / tilted | 1.072 | [1.0719, 1.0719] |
| pi/8 v=0.999: paper-best / tilted | 1.072 | [1.0719, 1.0719] |
| pi/8 v=0.99: sup_hexagon / tilted | 1.142 | [1.1425, 1.1425] |
| pi/8 v=0.99: paper-best / tilted | 1.142 | [1.1425, 1.1425] |
| pi/8 v=0.98: sup_hexagon / tilted | 1.115 | [1.1147, 1.1147] |
| pi/8 v=0.98: paper-best / tilted | 1.115 | [1.1145, 1.1146] |
| pi/8 v=0.999: sup_hexagon / full | 0.997 | [0.99683, 0.99683] |
| pi/8 v=0.999: tilted - CHSH (bits) | 0.55 | [0.5500, 0.5500] |
| pi/8 v=0.999: sup_hexagon - tilted (bits) | 0.06 | [0.0558, 0.0558] |
| pi/8 v=0.999: sup_hexagon - CHSH (bits) | 0.61 | [0.6058, 0.6058] |
| pi/8 v=0.999: sup_hexagon - paper search (bits) | – | [3.18e-06, 3.26e-06] |
| 27.76deg v=0.999: sup_hexagon / CHSH | 2.5 | [2.4890, 2.4890] |
| 27.76deg v=0.99: sup_hexagon / CHSH | 1.9 | [1.8560, 1.8560] |
| 27.76deg v=0.98: sup_hexagon / CHSH | 1.5 | [1.4536, 1.4536] |
| 27.76deg v=0.999: tilted / CHSH | 2.4 | [2.4060, 2.4060] |
| WBC pi/6 AB v=0.999: centre - half-way (bits, >0 required by Cor.1(b)) | >0 | [0.0360, 0.0360] |
| WBC pi/6 AB v=0.99: centre - half-way (bits, >0 required by Cor.1(b)) | >0 | [0.4638, 0.4638] |
| P_T A0 v=0.999: CHSH - half-way (bits, >0 by Cor. chshopt) | >0 | [0.0033, 0.0033] |
| P_T A0 v=0.99: CHSH - half-way (bits, >0 by Cor. chshopt) | >0 | [0.0115, 0.0115] |
| P_T A0 v=0.98: CHSH - half-way (bits, >0 by Cor. chshopt) | >0 | [0.0174, 0.0174] |
| P_T A0 v=0.95: CHSH - half-way (bits, >0 by Cor. chshopt) | >0 | [0.0349, 0.0349] |
| P_T A0 v=0.9: CHSH - half-way (bits, >0 by Cor. chshopt) | >0 | [0.1129, 0.1129] |


## 3. Answers to the specific questions

Notation: target p_θ (|ψ_θ⟩ = cosθ|00⟩ + sinθ|11⟩, A₀ = Z, A₁ = X, B_y = cosμ Z ± sinμ X, cosμ = 1/r,
r = √(1+sin²2θ)); white noise p_v = v p_θ; every witness normalised to g·p_θ = 1, so the score threshold is ω = v.
"Certified" = two-sided interval from Lemmas 1–5 with exact-rational certificates (widths ≤ 10⁻⁶ unless shown);
"exact" = closed form or exact zero.

**Theorem 1 (CHSH, exact).** For every quantum behaviour with CHSH value ≥ S,
P_g(A_x|E) ≤ ½[1 + √(2 − S²/4)], tight. On p_θ, CHSH = 2r, so the CHSH row of Table `tab:rates` is
H_min = −log₂ ½[1 + √(2 − r²v²)] exactly: 0.226659, 0.210974, 0.194352, 0.148883, 0.084669 (π/8);
0.350085, 0.326529, 0.302144, 0.237970, 0.152205 (27.76°), for v = 0.999, 0.99, 0.98, 0.95, 0.90.
*Proof.* Pironio et al. 2010 (closed form and qubit attack); our certificates reproduce it to 10⁻⁸ (level 1+AB). ∎

**Theorem 2 (tilted CHSH, certified).** Let g_θ = (𝔞A₀ + CHSH)/√(8+2𝔞²), 𝔞 = 2cos2θ/r. For every quantum behaviour
(any dimension) with g_θ·z ≥ v, H_min(A₁|X=1,E) ≥ H_θ(v), and an explicit two-qubit attack achieves H_θ(v) + 10⁻⁶, with
H_{π/8}(v) = 0.776634, 0.383682, 0.185275 at v = 0.999, 0.99, 0.98 and H_{27.76°}(v) = 0.842291, 0.548525, 0.393120,
0.118590 at v = 0.999, 0.99, 0.98, 0.95. H_θ(v) = 0 exactly for v ≤ h_L(g_θ) = (cos2θ + r)/2 (0.96593 at π/8, 0.93104 at
27.76°). *Proof.* Lemma 1 at the stated μ (stored in certs/), with h_Q(Π_e + μ g̃_θ) bounded by a **level-1+AB**
(9×9) rational certificate (Lemma 3, 4); upper side Lemma 5 with the stored attack; zeros by a deterministic point. ∎
(No closed form is claimed; the level-1+AB certificate is tight to 10⁻⁸ at every point, which suggests one exists.)

**Theorem 3 (the best witness in the exact normal cone, certified — replaces the "best found" search).**
Let C_θ be the hexagon of Theorem `thm:tilted` (θ = π/8, 27.76°). Then
  sup_{g∈C_θ} H_min(A₁|X=1,E; g, v) ∈ [0.832443, 0.832443] (π/8, v = 0.999), 0.438352 (0.99), 0.206530 (0.98),
  and 0.871353, 0.606026, 0.439191, 0.126918 for 27.76° at v = 0.999, 0.99, 0.98, 0.95;
  it is 0 exactly below the window threshold of Corollary `cor:window`.
The supremum is attained (to 10⁻⁶) by an explicit σ-symmetric g* ∈ C_θ, which coincides (to 10⁻³ in ℓ₂) with the
paper's search result; the paper's search was within 3.3·10⁻⁶ bits of the true supremum.
*Proof.* Lower side: Lemma 1 at g* (rational convex weights of the six arb-exact vertices; level-2 certificates, one
level-3). Upper side: Lemma 6 — one stored qubit/deterministic mixture satisfies V_i·z̄ ≥ v for all six vertices, so
it attacks every g ∈ C_θ simultaneously. ∎
**Corollary 3.1 (the in-cone gain over tilted CHSH is real for H_min).** At θ = π/8 the best in-cone witness certifies
**7.19%, 14.25%, 11.47%** more min-entropy than g_θ at v = 0.999, 0.99, 0.98 (strict, rigorous: e.g. 0.832443 ≥ 1.0719 ×
0.776634); at 27.76°: 3.45%, 10.5%, 11.7%. Against CHSH: 3.6727×, 2.0778×, 1.0627× (π/8) and 2.489×, 1.856×, 1.454×
(27.76°). The sup reaches **99.683%** of the true full-statistics min-entropy (0.835087–0.835088) at π/8, v = 0.999 — a
statement about the exact quantities, no longer about relaxations.

**Theorem 4 (full statistics, certified).** H_min(A₁|E) given z = v p_θ: π/8: [0.835087, 0.835088], 0.514289,
0.376273, 0.212696, 0.103973; 27.76°: 0.871943, [0.630677, 0.630678], 0.505821, 0.322627, 0.186098 (v = 0.999…0.90).
*Proof.* Lemma 2 with the LP-dual functional f (rational), certificates at levels 1+AB/2/3 (level 3 only needed at
π/8, v = 0.999); attack = nine-component qubit mixture reproducing p_v exactly up to a local remainder of weight 10⁻⁸. ∎

**Theorem 5 (WBC panel and P_T, certified).** WBC φ = π/6, H_min(A₀B₀|00,E): centre I_φ 1.729981, 1.069929, 0.588463,
0, 0; half-way to W₁ [1.693991, 1.693992], 0.606132, 0, 0, 0; vertex W₁ 0; CHSH 0.439524, 0.409375, 0.378186, 0.296221,
0.187197; full [1.759378, 1.759379], [1.163943, 1.163944], 0.846721, 0.475418, 0.248684. The centre exceeds the
half-way point by ≥ 0.0360 bits (v = 0.999) and ≥ 0.4638 (v = 0.99), as Corollary 1(b) requires. P_T, H_min(A₀|E):
CHSH 0.911547, 0.737568, 0.642251, 0.472341, 0.307179 > half-way 0.908267, 0.726020, 0.624817, 0.437441, 0.194229
(strictly, at every v), full 0.911632–0.911633, 0.739922, 0.648357, 0.492569, 0.336009 — rigorous confirmation of
Corollary `cor:chshopt`. All certificates are level 1+AB (one level 2). ∎

**Proposition 6 (von Neumann rates; the in-cone gain is absent).** For H(A₁|X=1,E) (the entropy entering EAT rates),
explicit mixed-state qubit attacks (exact rational ρ via a Cayley transform; Eve purifies; entropies in ball arithmetic)
give, for **every** g ∈ C_θ simultaneously, H ≤ 0.970652, 0.706523, 0.413045, 0 (π/8; v = 0.999, 0.99, 0.98, 0.95) and
H ≤ 0.985497, 0.854991, 0.709982, 0.274955 (27.76°). These equal the "chord" (v − h_L)/(1 − h_L) of tilted CHSH, i.e.
the attack mixes the ideal target with local points. Since the certified BFF lower bounds for g_θ (../rates_v11) are
0.96433, 0.69880, 0.40752 (π/8) and 0.97772, 0.83850, 0.69232, 0.26643 (27.76°), **no witness in the cone beats tilted CHSH
by more than 0.006–0.008 bits (π/8) / 0.008–0.018 bits (27.76°)** in von Neumann entropy, versus 0.056 bits (7.2%) in
min-entropy. The remaining gap is the BFF slack (Gauss–Radau truncation); we make no claim about its size.

**Remark (scope).** All statements are single-round (asymptotic, i.i.d.-attack) entropies under white noise. They are
not finite-size rates; for those see ../rates_v11 (EAT). The certified one-sided white-noise envelope there shows that
the maximally entangled CHSH experiment (θ = π/4) certifies more H(A₁|E) per round than any tilted target with any
witness, including full statistics, at every v tested (0.86–0.999), for both H_min and von Neumann; the large gains
above are **same-state** gains (witness matched to a given partially entangled source vs the CHSH score on that source).



## 4. How to verify

```
pip install python-flint          # only dependency of the verifier (no SDP/LP solver needed)
cd /workspace/ndwu/rigorous_rates
python3 verify_all.py             # re-checks every certificate in certs/ exactly; writes out/verified.json
python3 make_table.py             # rebuilds TABLE.md from out/verified.json
python3 verify_vn.py              # re-checks the von Neumann attack upper bounds in certs_vn/
```
Files: `npa_struct.py` (moment-matrix combinatorics), `exact_psd.py` (rational LDLᵀ), `families_arb.py` (targets and
functionals in ball arithmetic), `verify_all.py` (standalone verifier). Generation (not needed for verification):
`cert.py`, `rates_core.py`, `qubit.py`, `run_cases.py`, `vn_attack.py`, `run_vn.py`.

Runtime: `verify_all.py` re-checks all 299 cases (315 exact-rational certificates: 262 at level 1+AB, 39 at level 2, 14 at
level 3, plus 299 attacks) in about 2 s on the box; `verify_vn.py` takes under 1 s. Negative tests (a tampered Z entry, β, μ
or v) are all rejected. Every NPA level-3 value printed in prl_v13 lies inside (or within rounding of) its interval;
the largest interval width is 1.7·10⁻⁶ bits.

## 5. Recommendation for Results 2–3 and Fig. 2

1. **Say which entropy, and drop "numerical estimates".** "…certifies 0.832 bits of min-entropy of Alice's outcome A₁
   per round (Eve's guessing probability, white noise of visibility v), against 0.227 bits for CHSH. These values are
   rigorous: each is bracketed by an exact-rational sum-of-squares certificate valid in every dimension and an explicit
   qubit attack, which agree to 10⁻⁶ (Appendix, ancillary verifier)."
2. **Replace "the best expression we found" by the certified optimum.** "The best expression in the hexagon certifies
   0.832 bits — exactly the supremum over all target-maximising expressions, certified by one attack that defeats every
   expression in the hexagon simultaneously." The 7.2/14.2/11.5% and 3.7/2.1/1.06× numbers stand as printed.
3. **99.7%:** "…reaches 99.7% of the min-entropy certified by the full statistics" (now a statement about exact values;
   the earlier 'same relaxation level' qualifier is no longer needed).
4. **Qualify the in-cone gain as a min-entropy effect** (one sentence): "For the von Neumann entropy that enters
   entropy-accumulation rates, tilted CHSH is optimal in the hexagon to within 0.01 bits."
5. **Keep the same-state framing explicit** (one sentence, Results 2 or Connection to experiments): "These are gains
   for a given partially entangled source. Under white noise alone a maximally entangled CHSH experiment certifies more
   one-sided randomness per round at every visibility; the case for tilted sources rests on detection efficiency and
   two-party randomness, not on visibility." This pre-empts the obvious referee objection.
6. **Fig. 2 caption:** replace "numerical NPA level-3 lower-bound estimates … not used in any proof" with "rigorous:
   at each plotted point the exact single-round min-entropy lies within 10⁻⁶ of the curve (exact-rational SOS
   certificate + explicit qubit attack, ancillary `verify_all.py`); between points the curves are interpolations."
   Diamonds become "best expression in the hexagon (certified supremum)". The CHSH curve is the closed form of Thm 1.
7. **Appendix "Programs/Cross-checks"**: replace by a short "Certification" paragraph (Lemmas 1–3, 5, 6 in five
   lines) and keep the Clarabel description only as how certificates were found. Remove "Tolerances of 10⁻¹⁰ failed…"
   and "Level 2 … is not tight there and we do not use it": the lowest certificate levels that are tight are 1+AB
   (tilted, CHSH, WBC, P_T), 2 (in-cone optima) and 3 (only full statistics at π/8, v = 0.999).
