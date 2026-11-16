"""Visibility-optimised envelope over theta (vN from certified BFF, H_min from NPA level 2) vs CHSH at theta=pi/4."""
import numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
import analyze as A
from families import ptheta
from eat import chsh_analytic
VS = [0.999, 0.994, 0.99, 0.98, 0.97, 0.95, 0.92, 0.90, 0.86]
TH = ['pi/16', '0.1pi', 'pi/8', '27.76deg', '3pi/16', '0.22pi', 'pi/4']
def hmin_chsh(S): return -np.log2(0.5 + 0.5*np.sqrt(max(2 - S*S/4, 0)))
T = A.asym_table(); H = json.load(open('out/hmin_env.json')) if os.path.exists('out/hmin_env.json') else []
out = []; L = ['| v | CHSH max-ent vN | best tilted vN (theta) | best full vN (theta) | CHSH max-ent Hmin | best tilted Hmin (theta) | best full Hmin (theta) |', '|' + '---|'*7]
for v in VS:
    rec = dict(v=v, vN_chsh=float(chsh_analytic(2*np.sqrt(2)*v)), Hmin_chsh=float(hmin_chsh(2*np.sqrt(2)*v)), vN={}, Hmin={})
    for lab in TH:
        d = [x for x in T if x['label'] == lab and abs(x['v'] - v) < 1e-12]
        if d: rec['vN'][lab] = dict(tilted=d[0]['tilted'] if lab != 'pi/4' else rec['vN_chsh'], full=d[0]['full'], CHSH_same_state=d[0]['CHSH'])
        h = [x for x in H if x['label'] == lab and abs(x['v'] - v) < 1e-12]
        if h: rec['Hmin'][lab] = dict(tilted=h[0]['tilted'], full=h[0]['full'], CHSH_same_state=h[0]['CHSH'])
    def best(dic, key):
        c = [(val[key], lab) for lab, val in dic.items() if val.get(key) is not None]
        return max(c) if c else (None, None)
    rec['vN_best_tilted'] = best(rec['vN'], 'tilted'); rec['vN_best_full'] = best(rec['vN'], 'full')
    rec['Hmin_best_tilted'] = best(rec['Hmin'], 'tilted'); rec['Hmin_best_full'] = best(rec['Hmin'], 'full')
    out.append(rec)
    f = lambda t: '-' if t[0] is None else '%.4f (%s)' % t
    L.append('| %.3f | %.4f | %s | %s | %.4f | %s | %s |' % (v, rec['vN_chsh'], f(rec['vN_best_tilted']), f(rec['vN_best_full']), rec['Hmin_chsh'], f(rec['Hmin_best_tilted']), f(rec['Hmin_best_full'])))
json.dump(out, open('out/envelope.json', 'w'), indent=1); open('out/table_envelope.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
# per-theta detail table
L = ['| theta | v | vN: CHSH same state | vN: tilted | vN: full | Hmin: CHSH same state | Hmin: tilted | Hmin: full |', '|' + '---|'*8]
g = lambda x: '-' if x is None else '%.4f' % x
for lab in TH:
    for r in out:
        a = r['vN'].get(lab, {}); b = r['Hmin'].get(lab, {})
        if a or b: L.append('| %s | %.3f | %s | %s | %s | %s | %s | %s |' % (lab, r['v'], g(a.get('CHSH_same_state')), g(a.get('tilted')), g(a.get('full')), g(b.get('CHSH_same_state')), g(b.get('tilted')), g(b.get('full'))))
open('out/table_theta_detail.md', 'w').write('\n'.join(L) + '\n')
