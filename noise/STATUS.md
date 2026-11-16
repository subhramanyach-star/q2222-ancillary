# noise STATUS (started 2026-10-23 14:50 HKT)
Question (RR): does the certified H_min(A1|X=1,E) gain of tilted / C_theta' witnesses over the best CHSH protocol survive
realistic source noise at lab efficiencies? Grid eta in {0.80,0.84,0.88,0.92}; noise (a) white noise visibility v,
(b) dephasing (coherence factor p on |00><11|, = mixture of the tilted state with its dephased version).
Each protocol re-optimised (own best state and settings) at each (eta, noise). Do NOT edit prl or files outside noise/.
## Log
- 14:50 setup. Reusing ../core.py, ../cert_eta.py, ../cone.py, ../../rigorous_rates.
- 14:58 nz.py (optimiser, float NPA 1+AB) validated: white eta=.88 v=.99 reproduces earlier certs (chsh .115701, tilted .120388, cone .12616 vs cert .12617). Launched 8 procs: white v in .995..0.93, deph p in .995..0.85, eta .80,.84,.88,.92 -> out/opt_{white,deph}.jsonl, logs out/log_*.txt
- 14:50 ncert.py (certifier; CE.winner with model-specific honest_arb). Tested: white .88 v=.995: CHSH<=.121673 (S_up = honest
  optimum to 1e-8, tight), tilted>=.127024, cone>=.133499. deph .88 p=.995: certified CHSH bound (p*any state+(1-p)*separable;
  loose) .126971 vs honest dephased CHSH optimum .125139 (float); tilted>=.131046, cone>=.138035 (both > loose bound => proved).
- 14:52 added eta=0.76 (covers Shalm/Giustina/Bierhorst ~0.75-0.79): white .995..0.98, deph .995..0.96.
  Experiment parameters collected (see EXPERIMENTS section at end, filled at write-up).
- 15:12 BUG: polygon cache cleared after insert -> KeyError in all nz.py procs (~14:58-15:05). Fixed (clear before insert); relaunched all 10 (resumable).
- 15:25 Float grid mostly done (see out/tables.md via mk.py). White noise erodes the gain (eta .80: +4.1/+6.7% at v=.995 ->
  +1.2/+2.1% at v=.97; eta .76: gain -> 0 together with CHSH at v~.98). Dephasing barely changes it (+4-5% / +7-10%).
  Dephasing CHSH side: no tight DI certificate (an Alice-side-involution NPA model is useless: D=1 allowed); keep the
  loose certified bound (p*any + (1-p)*separable) + float dephased optimum.
- 15:27 model 'wd' (white v AND dephasing p; HV visibility = v, max-ent DA visibility = v*p) for experiment-matched points:
  Shalm15 (.7515,.999,.997), Giustina15 (.774,.99,1), Bierhorst18 (.755; vis assumed = Shalm15), Liu18 (.7865,.995,.990),
  Liu21 NatPhys (.841,.994,.991), Shalm21 NatPhys (.763; vis assumed), Liu22 PRL (.8749; F=.9952 only -> brackets
  white v=.9936 and deph p=.9768 at its theta=20deg). Running nz.py wd (7 procs). ncert.py patched for wd.
- 15:37 white + wd optimisation done; certification launched (ncert.py white per eta; wd 2 procs). deph opt still running (q .90/.85).
- 15:34 RESUMED after interruption (note: earlier log times ran ~5 min ahead of box clock). Alive: ncert white .80-.92, nz white .88 extension. Done: opt deph, opt wd, cert wd (8/8), white .80 extension (v=.958). Launched ncert deph (5 procs).
- 15:47 ALL DONE. opt: white 115, deph 123, wd 24 records; certs: white 40, deph 41, wd 8 (no certification errors).
  Tables: out/tables.md (mk.py). Summary write-up: REPORT.md.
## EXPERIMENTS (photonic loophole-free, all SPDC + SNSPD, all with non-maximally entangled states)
| experiment | system eff. (A/B) | state | visibility / fidelity | source |
|---|---|---|---|---|
| Shalm et al. PRL 115, 250402 (2015) | 74.7% / 75.6% | 0.961|HH>+0.276|VV> | max-ent vis. 0.999 (H/V), 0.996 (D/A) | https://arxiv.org/abs/1511.03189 |
| Giustina et al. PRL 115, 250401 (2015) | 78.6% / 76.2% | Eberhard r ~ -2.9 | vis. >99% in H/V and diagonal | https://arxiv.org/abs/1511.03190 |
| Bierhorst et al. Nature 556, 223 (2018) | 75.5% symmetric | 0.982|HH>+0.191|VV> | not reported in main text (NIST source as 2015) | https://arxiv.org/abs/1803.06219 |
| Liu et al. Nature 562, 548 (2018) | 78.8% / 78.5% | cos22.05|HV>+sin22.05|VH> | vis. 99.5% (H/V), 98.5% (D/A); F = 99.02% | https://arxiv.org/abs/1807.09611 |
| Liu et al. Nat. Phys. 17, 448 (2021) | 83.40% / 84.80% | cos27.76|HV>+sin27.76|VH> | vis. 99.4% (H/V), 98.5% (D/A); F = 99.06% | https://arxiv.org/abs/1912.11159 |
| Shalm et al. Nat. Phys. 17, 452 (2021) | 76.3% symmetric | cos14.8|HH>+sin14.8|VV> | not reported (NIST source) | https://arxiv.org/abs/1912.11158 |
| Liu et al. PRL 129, 050502 (2022) | 87.16% / 87.82% | cos20|HV>+sin20|VH> | F = 99.52 +- 0.15% | https://arxiv.org/abs/2110.01480 |
Model mapping: v = H/V visibility, p = (D/A)/(H/V); Liu22: F only -> bracket white v=(4F-1)/3=.9936 / dephasing p=.9768.
Not modelled: background/dark counts (NIST threshold 72.5%), multi-pair emission, asymmetric efficiencies (symmetric mean used).
