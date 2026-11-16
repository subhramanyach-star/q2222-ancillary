"""Rigorous UPPER bounds on the DI quantity H(AB|X=x*,Y=y*,E) for a white-noise tilted source: the honest implementation
with Eve holding a purification of rho_v = v psi_theta + (1-v) 1/4 is one admissible attack, and for rank-1 projective
measurements on AB Eve's conditional states are pure, so H(AB|E) = H_Shannon(p_ab) - S(rho_v).
U_key: the key pair actually used (A_{x*} in {sz, sx}, B2 optimised in the xz plane) = targets2.best_key.
U_real: max over ALL real (xz-plane) qubit measurement pairs (alpha, phi).  Writes out/global/upper_global.json."""
import numpy as np, json
from targets2 import best_key
from run_global_cfg import TH, VS
def S_rho(v): ev = np.array([v + (1 - v) / 4] + [(1 - v) / 4] * 3); ev = ev[ev > 0]; return float(-(ev * np.log2(ev)).sum())
def Hsh(p): p = p[p > 1e-15]; return float(-(p * np.log2(p)).sum())
def U_real(th, v):
    C, s = np.cos(2 * th), np.sin(2 * th); best = 0
    for al in np.linspace(0, np.pi, 361):
        for ph in np.linspace(0, np.pi, 361):
            a, b, e = C * np.cos(al), C * np.cos(ph), np.cos(al) * np.cos(ph) + s * np.sin(al) * np.sin(ph)
            p = v * np.array([(1 + sa * a + sb * b + sa * sb * e) / 4 for sa in (1, -1) for sb in (1, -1)]) + (1 - v) / 4
            best = max(best, Hsh(p))
    return best - S_rho(v)
out = []
for lab, th in TH:
    for v in VS:
        h, xs, phi = best_key(th, v)
        out.append(dict(label=lab, v=v, U_key=h - S_rho(v), U_real=U_real(th, v), S_rho=S_rho(v)))
        print(out[-1], flush=True)
json.dump(out, open('out/global/upper_global.json', 'w'), indent=1)
