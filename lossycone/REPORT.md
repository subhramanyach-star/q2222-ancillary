# Exact optimisation of Bell witnesses under loss — rigorous results (23 Oct 2026, 18:20 HKT)
Full statements/proofs: theorem.tex / theorem.pdf (5 pp). Data: data/t4/<eta>.json, t4_summary.json. Verifier: vpt.py (Arb), tamper.py (negative tests).

## Proved
**Pen and paper (T3).**
- Lemma 1: P_g(g,ω) = sup{P_g(q): g·q ≥ ω}; concavity; once a device certifies anything, P_g(g,·) drops at least linearly with slope κ = (1−P_g)/(ω*−h_L(g)) > 0.
- Lemma 2: at every interior point of Q there is a tight witness (one Bell value certifies as much as the full statistics).
- **Theorem 1:** if p* locally minimises P_g on a convex honest set K (e.g. Q_η) with P_g(p*) < 1, every tight witness g lies in N_K(p*). No assumptions on monotonicity. The same holds for any convex entropy, including von Neumann H(A|E).
- Corollary (improvement test) and Prop. (ε-tight version: g·p ≤ g·p* + ε/κ).
- Prop. (reduction): the optimum over Q_η sits at L_η(q*) with q* a nonlocal extreme (two-qubit) point, N_{Q_η}(L_η q) = M_η^{−T} N_Q(q), so the pull-back of the optimal witness lies in C(q*).

**Computer-assisted (T1, T2) at η = .75 .80 .85 .90 .95 .98:** N_Q(q_η) = C(q_η) exactly (5/5/7/7/6/8 vertices, all μ > 0, each vertex tight on exactly 2 deterministic points), and the cone witness certifies:

| η | H_min (bits) | CHSH bound | gain |
|---|---|---|---|
| .75 | [.009845,.009847] | .009013 | +9.2% |
| .80 | [.037859,.037864] | .033445 | +13.2% |
| .85 | [.095104,.095105] | .081930 | +16.0% |
| .90 | [.198626,.198627] | .168313 | +18.0% |
| .95 | [.385565,.385566] | .330204 | +16.7% |
| .98 | [.597367,.597368] | .525646 | +13.6% |

**T4 (continuum in η, rigorous):** all 110 grid points in [0.75,0.98] certified (643 vertices; 68.8k Arb boxes). Loss emulation (L_a∘L_b = L_ab ⇒ Q_η ⊆ Q_η' for η ≤ η') gives a gap-free statement: on [η_i,η_{i+1}] run the η_i protocol with clicks discarded, and bound CHSH by S_up(η_{i+1}). Guaranteed gain for every η: ≥5.3% [.75,.80], ≥8.5% [.80,.85], ≥10.7% [.85,.90], ≥12.2% [.90,.95], ≥8.0% [.95,.98]. At the grid points the gain is 8.8–17.6%. Between grid points the witness is not in the normal cone; it is only a valid lower bound.

## Still numerical
Only the device optimisation. It is not proved that our devices are optimal, so the rates are certified lower bounds on the best achievable rate. The cone statements hold pointwise at the rational devices. The fact that the cone-best witness matches the full-statistics optimum (≤1e-6 at 5/6 points) is a numerical observation.
Scope: symmetric η, independent losses, no-click → +1, H_min(A_1|X=1,E), comparison against CHSH-value-only certification.
Caveat: at η=.95 the L3 optimiser drifts to a Hardy-type device (NPA lower bound 0.047 vs attack 0.389). We use the interior device warm-started from η=.949 (0.3856 bits).
