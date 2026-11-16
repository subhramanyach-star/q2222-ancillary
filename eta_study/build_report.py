"""Assemble REPORT.md from report_src/*.md and the auto-generated out/tables.md (run make_tables.py first)."""
import re
S=lambda f: open('report_src/'+f).read()
tab=open('out/tables.md').read().replace('## Table','### Table')
figs="""
![PRL figure](figs/fig_eta_prl.png)

*Proposed PRL figure (`figs/fig_eta_prl.pdf`).* (a) H_min(A₁|X=1,E) vs η at v = 1. Lines show numerical optima; bars are proved upper bounds for every CHSH protocol; circles and squares are proved lower bounds for the tilted and C_θ′ witnesses. (b) Proved gain over every CHSH protocol at v = 1 and v = 0.99.

![v=1](figs/fig_eta_hmin_v1.png)
![v=0.99](figs/fig_eta_hmin_v099.png)
![vN](figs/fig_eta_vn_v1.png)
![AB](figs/fig_eta_hmin_AB_v1.png)
"""
notes="""
Column notes:
- "every CHSH protocol ≤" is the certified upper bound H(S_up(η)), valid for any CHSH protocol, any state, any measurements and any dimension.
- Bracketed entries are two-sided intervals wider than 5×10⁻⁷.
- Gains are relative to that upper bound.
- θ′ is the target of the best expression found in C_θ′.
- The full-statistics columns are numerical (NPA 1+AB).
- n/a† marks an ill-conditioned full-statistics dual (|κ| ~ 10³), for which the EAT second-order term dominates.
"""
out=S('head.md')+'\n'+S('ranges.md')+'\n## 4. Tables\n'+notes+'\n'+tab+'\n\n### Figures\n'+figs+S('tail.md')
open('REPORT.md','w').write(out); print(len(out))
