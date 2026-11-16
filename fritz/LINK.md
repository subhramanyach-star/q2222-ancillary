# Fritz and our normal cones: what is proved where (23 Oct 2026, HKT)

## 1. What Fritz proved (JMP 53, 072202 (2012), arXiv:1202.0141; read in full)
- **Main theorem (Thm 4):** in every (n,2,2) scenario, the no-signalling cone is linearly isomorphic to the dual of the
  local (Bell) cone, via F^{⊗n}, where F is the substitution A0→A0+A1, A1→A0−A1. **Marginals are included**: he
  homogenises with a "no measurement" setting. So tight Bell inequalities ↔ extremal NS boxes, and CHSH ↔ PR box.
  This is the polytope duality Local ↔ NS. It is not a description of the quantum set.
- Prop 7: no such duality for (k,l)≠(2,2). Sec. 3: n→n+1 extension (Mermin–Klyshko), Werner–Wolf/Żukowski–Brukner as a special case.
- **Quantum set:** only Question 14 (is Q self-dual under the same map?). Prop 15: a symmetric self-dual cone has
  CHSH max 2√2. Prop 16: false for n=3. n=2: "holds on the level of correlations", citing an unpublished private communication from
  Palazuelos–Villanueva. Marginals are left open.
- Other Fritz work checked: temporal CHSH (NJP 12, 083055 (2010)): temporal = spatial quantum region **for correlators only**.
  Local orthogonality (Nat. Commun. 4, 2263 (2013)) and Acín–Fritz–Leverrier–Sainz (CMP 334, 533 (2015)) give principles
  and relaxations, not exact descriptions or normal cones. **No Fritz paper characterises the full Q_2222 by duality.**
- Le–Meroni–Sturmfels–Werner–Ziegler (Quantum 7, 947 (2023)): the first published proof of the correlator self-duality
  Q°=½HQ (the 4-dim zero-marginal body), the support function (Tsirelson/Landau/Masanes), the "normal cycle" (= correlator
  normal cones) and the unique exposing functional on the exposed patches. They state that their techniques do not
  carry over to the full-marginal 222 case. (Their abstract mentions "four marginals and four correlations", but the body they study
  is the 4-dim correlator body.)

## 2. The link, made precise
Homogenise z ↦ x=(1,z) ∈ R³⊗R³. A Bell expression g with bound h is the covector f = h·e₀₀ − g, and N_Q(p) ≅ {f∈K_Q*: f(x_p)=0}
(the face of the cone of Tsirelson inequalities exposed by p). Our polygons are its h=1 sections.
- (i) **Deterministic cuts = Fritz.** K_Q* ⊆ K_L* = Φ(K_NS). So every normalised g∈N_Q(p) satisfies g·L≤1, i.e.
  Φ⁻¹(e₀₀−g) is a valid NS box. Tight cuts correspond to zero entries of that box. Check: the images of the BB octagon
  vertices g0, g1 are NS boxes with exactly 2 zero probabilities, matching their faces conv{P_T,L,L'} (fritzmap2.py).
  So Fritz gives the *outer polytope* that contains every normalised quantum normal cone. He does not say which of
  those functionals are quantum-maximised at p. In our work the slice cone and overlap stationarity decide that,
  with sufficiency from the NDWU dual.
- (ii) **Correlators: partly implied by prior work.** Q_corr is the projection of Q. So a correlator-only g is maximised at p over Q
  iff E(p) maximises it over Q_corr, and that is decided by Tsirelson's support function (proved in Le et al. Thm 2).
  Our Thm "correlator functionals" (general/Thm A) is therefore implied as a *bound*. What is ours is the local
  criterion (stationarity in Bob's overlap suffices) and that it holds at points with nonzero marginals. At P_T:
  N_Q(P_T)∩{a=b=0} = CHSH ray, the centre of the octagon, consistent with Le et al.
- (iii) **Marginals: Fritz's map is not a self-duality (Question 14 negative for n=2 with marginals).**
  **Priority:** Barizien–Bancal (PRL 133, 010201 (2024)) already make the extremality observation (quote below). The exact
  separating certificate below is ours and independent. See Sec. 5.

## 2b. Proposition and proofs
**Proposition.** Let n=2, with two binary observables per party and marginals included. The image of the cone of Tsirelson inequalities K_Q* under
Fritz's map Φ: g_{s1s2} ↦ F^{s1t1}F^{s2t2}g_{t1t2}, with F=[[1,0,1],[0,1,0],[1,0,−1]] in the index order (A0, no measurement, A1), is not K_Q.
(Relabelling the two observables, or the outcomes, conjugates Φ by a symmetry of Q, so the statement does not depend on the labelling.)

*Proof 1 (analytic, extremality; this is Barizien–Bancal's remark).* Φ is linear and invertible, so Φ(K_Q*)=K_Q would make Φ
map extreme rays to extreme rays. Take the Tsirelson inequality f = 2√2·e00 − CHSH ∈ K_Q*. Its image is the ray of the Tsirelson
point, up to flipping Bob's outcomes (Fritz Prop. 15, rechecked here). That point is an extreme point of Q. But f is not extreme in K_Q*.
Let g± = CHSH/(2√2) ± ρ((A0+A1)/√2 − B0) with ρ = 1 − 1/√2. These are two opposite vertices of the Barizien–Bancal octagon, and each has quantum maximum 1
[BB 2024; our Thm octagon]. Then f/(2√2) = ½[(e00 − g+) + (e00 − g−)], with both summands in K_Q* and not proportional to f. ∎

*Proof 2 (explicit, exact, independent of the octagon): a Tsirelson inequality whose image is not quantum.*
- g: 14/17·A0 + CHSH ≤ 52/17. This is valid on Q: the quantum maximum of αA0+CHSH is √(8+2α²) (Acín–Massar–Pironio, PRL 108, 100402 (2012)),
  and here 8+2(14/17)² = (52/17)² exactly. It is also valid on L, since the local maximum is 2+α = 48/17.
- The image is y⁰⁰ = 52/17 > 0. Normalised: z* = (α0,α1,β0,β1,E00,E01,E10,E11) = (−7/26, −7/26, 0, 0, −17/26, −17/26, −17/26, 17/26).
  It is a valid NS box (smallest probability 1/52), as it must be by Fritz's Thm 4.
- The separating functional is w·z = −(α0+α1) − β0 − (E00+E01+E10−E11). Equivalently, after flipping Alice's outputs:
  I = CHSH + A0 + A1 − B0, whose local maximum is 3.
- The certificate is an exact rational 9×9 matrix Y ≻ 0 on the NPA level-1+AB monomials {1,A0,A1,B0,B1,AxBy} such that
  61/20 − w·z = ⟨Y, Γ⟩ for every moment matrix Γ, identically in all 17 moment classes. Since Re Γ ⪰ 0 for every quantum
  realisation (any dimension; Re⟨W⟩ = Re⟨W†⟩), we get w·q ≤ 61/20 on Q.
- But w·z* = 41/13 > 61/20, a gap of 27/260 ≈ 0.104. The NPA optimum is ≈ 3.0179, so the true quantum maximum is ≤ 3.0179 < 41/13.
- Verification uses exact rationals only (Python `fractions`): `python3 verify_q14.py` checks the bound identity, the image,
  the NS check, the a=0 sanity check (CHSH ↦ Tsirelson point) and an exact LDL^T of Y with all pivots > 0. The data are in `cert_q14.json`
  (sha256 ad459d81…), and the builder is `make_cert.py`.
- NPA level 1 alone does not separate this point (visibility 1.038). Level 1+AB is needed, so I found no simpler pen-and-paper
  version of Proof 2. The one-sided Cauchy–Schwarz bound CHSH ≤ 2√2·√(1−m²) is satisfied by the image. Proof 1 is the clean analytic argument.
- Earlier float results (fritzmap_results.json): the NPA visibility of the image of the tilted-CHSH inequality is 0.973, 0.955, 0.948 and 0.967 at α = ¼, ½, 1 and 3/2,
  and 0.9465 at α = 14/17. For the octagon vertices it is 0.950.
- Still **open**: whether Q_2222 (with marginals) is self-dual under *some other* linear map. BB state this as open. Our polygonal cones would
  have to correspond to polygonal 2-faces of Q under any such map. That might give a proof, but I have not pursued it.

## 3. Overlap summary
- Implied by earlier work: the deterministic-cut constraint (Fritz, standard). Correlator-only optimality bounds and the zero-marginal
  normal cones (Tsirelson; Le et al.). The CHSH-point octagon (Barizien–Bancal).
- New in ours: exact normal cones in the full 8-dim set with marginals, at partially entangled points (Thm 2, all
  θ∈(0,π/4]). The local sufficiency criterion. Its use for choosing DI-randomness witnesses. The (numerical)
  exact rational certificate that Fritz's map fails with marginals (the failure itself was already observed by Barizien–Bancal 2024).

## 4. Drafts
See link.tex: an intro/discussion paragraph and an appendix remark. It compiles standalone; it references \ref{thm:tilted}, \ref{thm:corr} and \chsh from the PRL.

## Verified citations
- T. Fritz, Polyhedral duality in Bell scenarios with two binary observables, J. Math. Phys. 53, 072202 (2012), doi:10.1063/1.4734586, arXiv:1202.0141.
- T. Fritz, Quantum correlations in the temporal CHSH scenario, New J. Phys. 12, 083055 (2010), arXiv:1005.3421.
- T. P. Le, C. Meroni, B. Sturmfels, R. F. Werner, T. Ziegler, Quantum correlations in the minimal scenario, Quantum 7, 947 (2023), arXiv:2111.06270.
- K. T. Goh et al., Geometry of the set of quantum correlations, Phys. Rev. A 97, 022104 (2018), arXiv:1710.05892.
- T. Fritz et al., Local orthogonality as a multipartite principle for quantum correlations, Nat. Commun. 4, 2263 (2013).
- A. Acín, T. Fritz, A. Leverrier, A. B. Sainz, A combinatorial approach to nonlocality and contextuality, Commun. Math. Phys. 334, 533 (2015).
- Already in refs_v15: Cirel'son LMP 4, 93 (1980); Barizien–Bancal PRL 133, 010201 (2024); Masanes quant-ph/0512100.
- Not re-verified this run: Landau, Found. Phys. 18, 449 (1988); Masanes, quant-ph/0309137 (2003).

## 5. Literature verdict on Question 14 (n=2, marginals)
Verdict: **implicitly settled, in a remark.** For Fritz's specific map, Barizien–Bancal (2024) already note that it does not send Q to Q*.
Nobody has stated it as an answer to "Question 14", and no explicit non-quantum image was given before ours. Self-duality under an arbitrary map is **open**.
I checked all 23 papers citing arXiv:1202.0141 in Semantic Scholar (s2cites.json) and searched the texts of the relevant ones (lit/).
- Barizien & Bancal, PRL 133, 010201 (2024), arXiv:2401.12791: "It was recently shown that this subset is self-dual, i.e. Qc ≅ Q*c [18, 19].
  This striking property sets the quantum set apart from both the local and the no-signaling sets." And, in the conclusion: "In [18], a map was introduced to prove that
  L* ≅ NS and it was proven that this map also sends the subset Qc with uniformly random marginals statistics to its dual. However, since this
  map sends the extremal point PT onto βCHSH, which admits the decomposition Eq. (24), it cannot be used to map Q to Q*. While other maps could
  be found, analysis such as ours might help proving that the quantum set Q is not self-dual, making Q* a set of physical interest with a
  possibly very different geometrical structure." ([18] = Fritz 2012, [19] = Le et al.)
- Le et al., Quantum 7, 947 (2023): "We thought that self-duality (Thm. 3) was entirely our own idea, until we found it mentioned in [27,
  Sect.4] with only a private communication by Palazuelos and Villanueva (2011) as a reference. We could not trace a published proof." And: "There
  it was stated that the self-duality of Q was also known to C. Palazuelos and I. Villanueva (in a private communication), and that self-duality
  fails for the analogously defined three-partite body, but leaving the possibility for N > 3." (Their Q is the zero-marginal correlator body.)
- Fritz 2012, Sec. 4: "for n = 2, which is the first nontrivial case, this does turn out to be correct on the level of full (two-party)
  correlations, i.e. upon not taking into account the marginal probabilities [19]". Prop. 16: "The image of (Q_3)* under (23) is not equal to Q_3."
- Liu, Chung & Ramanathan, arXiv:2309.06304: "the duality of the quantum correlation set, which is convex but is in general not a polytope unlike
  L and NS, is less known and has been studied for some specific scenarios".
- Ramanathan, Tuziemski, M. & P. Horodecki, PRL 117, 050401 (2016): this answers Fritz's *other* open question (Sec. 4.2, whether any extremal nonlocal NS box is quantum)
  negatively in full generality. Combined with Fritz's own argument (Prop. 16), this rules out his self-duality in every scenario where a tight Bell inequality without
  quantum violation maps to a nonlocal extremal box. For n=2 there is none: every tight inequality is a positivity inequality or CHSH.
- Goh et al., PRA 97, 022104 (2018), and Mikos-Nuszkiewicz & Kaniewski (arXiv:2302.10658) do not cite Fritz 2012 and do not discuss self-duality.
  Ishizaka was not re-checked: my arXiv ID was wrong and the paper is not among the citing papers.
