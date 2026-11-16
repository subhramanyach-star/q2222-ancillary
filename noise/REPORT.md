# Does the randomness gain survive realistic source noise? (23 Oct 2026)
Files: nz.py (optimiser), ncert.py (certificates), mk.py -> out/tables.md (full tables), certs/*.json, STATUS.md (log, experiment parameters).

Model: symmetric efficiency eta (no-click -> +1). Source noise on cos th|00>+sin th|11>: (a) white noise, visibility v;
(b) dephasing, coherence p (mixture with the classically correlated state); (c) both (experiment-matched).
Each protocol has its own best state and settings at each (eta, noise): best CHSH protocol, alpha<A0>+CHSH, best expression in C_theta'.
Quantity: H_min(A1|X=1,E).

What is proved and what is numerical:
- White noise: the CHSH value is certified for any CHSH protocol whose source is v*(any state)+(1-v)*I/4, and it is tight (a qubit attack matches it).
  The tilted and C_theta' lower bounds are certified with exact NPA duals plus interval arithmetic. "Proved" = certified lower bound > certified CHSH upper bound.
- Dephasing: the tilted and C_theta' lower bounds are certified. On the CHSH side, the best dephased CHSH protocol is a float optimum (multi-start).
  The only certified CHSH bound is a loose one, valid when the source is p*(any state)+(1-p)*(any separable state). Where we beat that bound, the separation is proved; elsewhere the gain is numerical.

Findings:
- White noise reduces the gain steadily, and it goes to zero at the noise level where CHSH itself stops certifying.
  Approximate thresholds: v ~ 0.982 at eta=.76, 0.956 at .80, 0.92 at .84, 0.87 at .88.
  Example, eta=.88: C_theta' gives +9.7% at v=.995, +5.8% at .96 and +1.5% at .90.
- Dephasing (the dominant SPDC imperfection: H/V visibility > D/A visibility) barely reduces the gain.
  Example: eta=.88 goes from +10.3% to +9.4% at p=.96.
- At the reported parameters of the seven photonic experiments, the gain is +2.4..+4.6% (tilted) and +3.6..+9.7% (C_theta').
  It is proved at 7 of the 8 rows (Liu 2018 is numerical only).
