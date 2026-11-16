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
