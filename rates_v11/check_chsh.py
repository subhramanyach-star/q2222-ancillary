"""Cross-check of the certified BFF pipeline against the tight analytic CHSH bound (Pironio et al. NJP 11, 045021 (2009)):
H(A|X=x,E) >= 1 - h(1/2 + 1/2 sqrt(S^2/4 - 1)). Writes out/chsh_check.json."""
import numpy as np, json, warnings
warnings.filterwarnings('ignore')
from bffcert import Rate
from eat import chsh_analytic
g = np.array([0, 0, 0, 0, 1, 1, 1, -1.0]); out = []
R = {8: Rate(m=8, xstar=1), 12: Rate(m=12, xstar=1)}
for S in [2.8256, 2.80, 2.75, 2.70, 2.60, 2.50, 2.449, 2.40, 2.30, 2.20, 2.10, 2.05]:
    rec = dict(S=S, analytic=float(chsh_analytic(S)))
    for m in (8, 12):
        r = R[m].bound('score', g=g, omega=S); rec['cert_m%d' % m] = r['H_cert']; rec['num_m%d' % m] = r['H_num']
    out.append(rec); print(rec, flush=True)
    json.dump(out, open('out/chsh_check.json', 'w'), indent=1)
