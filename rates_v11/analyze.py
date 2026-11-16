"""Tables, EAT finite-size rates and figure from out/bff_main.jsonl, out/bff_cone.jsonl."""
import numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
from families import ptheta, gtilt, hL
from eat import eat_rate, chsh_analytic, chsh_tangent
def load(fn):
    return [json.loads(l) for l in open(fn)] if os.path.exists(fn) else []
MAIN = load('out/bff_main.jsonl'); CONE = load('out/bff_cone.jsonl')
LABELS = {'pi/16': np.pi/16, '0.1pi': 0.1*np.pi, 'pi/8': np.pi/8, '27.76deg': 27.76/180*np.pi, '3pi/16': 3*np.pi/16, '0.22pi': 0.22*np.pi, 'pi/4': np.pi/4}
gch = np.array([0, 0, 0, 0, 1, 1, 1, -1.0])
def best_main(lab, v, expr):
    rows = [d for d in MAIN if d['label'] == lab and abs(d['v'] - v) < 1e-12 and d['expr'] == expr and d['status'] == 'ok']
    if not rows: return None
    return max(rows, key=lambda d: d['H_cert'])
def best_cone(lab, v):
    rows = [d for d in CONE if d['label'] == lab and abs(d['v'] - v) < 1e-12]
    if not rows: return None
    d = rows[-1]; ms = [m for m in (8, 12) if d.get('H_cert_m%d' % m) is not None]
    m = max(ms, key=lambda m: d['H_cert_m%d' % m]); return dict(H_cert=d['H_cert_m%d' % m], K0=d['K0_m%d' % m], kappa=d['kappa_m%d' % m], g=d['g'], hL=d['hL'], lam=d['lam'], m=m)
def certs(lab, kind):
    """All affine certificates (K0, kappa, design v) available for protocol kind at theta label."""
    out = []
    if kind in ('tilted', 'full'):
        for d in MAIN:
            if d['label'] == lab and d['expr'] == kind and d['status'] == 'ok' and d['K0'] is not None: out.append((d['K0'], np.array(d['kappa']), d['v']))
    elif kind == 'cone':
        for d in CONE:
            if d['label'] != lab: continue
            for m in (8, 12):
                if d.get('K0_m%d' % m) is not None: out.append((d['K0_m%d' % m], np.array(d['kappa_m%d' % m]), d['v']))
    return out
def asym_table():
    rows = []
    for lab, th in LABELS.items():
        p, c0, (s, C, r) = ptheta(th)
        vs = sorted(set(d['v'] for d in MAIN if d['label'] == lab))
        for v in vs:
            S = v*2*r
            t = best_main(lab, v, 'tilted'); f = best_main(lab, v, 'full'); c = best_cone(lab, v)
            Ht = None if t is None else max(t['H_cert'], 0)
            # domination (rigorous): tilted g lies in C_theta, so best-in-cone >= tilted; any score-constrained bound is also a
            # bound for the full-statistics problem (smaller feasible set)
            Hc = None if c is None else max(c['H_cert'], Ht or 0)
            Hf = None if f is None else max(f['H_cert'], Ht or 0, Hc or 0)
            rows.append(dict(label=lab, theta_deg=th*180/np.pi, v=v, CHSH=float(chsh_analytic(S)) if S > 2 else 0.0,
                             tilted=Ht, cone=Hc, full=Hf, full_raw=None if f is None else f['H_cert'], cone_raw=None if c is None else c['H_cert'],
                             cone_lam=None if c is None else c['lam'], cone_hL=None if c is None else c['hL'],
                             maxent_CHSH=float(chsh_analytic(2*np.sqrt(2)*v)),
                             hL_tilted=float(hL(gtilt(th))*1.0), v_CHSH_thr=1/r, v_tilted_thr=(C + r)/2))
    return rows
def eat_curves(lab, v, ns, kinds=('CHSH', 'tilted', 'cone', 'full'), window=0.03):
    th = LABELS[lab]; p, c0, (s, C, r) = ptheta(th); z = v*p; out = {}
    for kind in kinds:
        if kind == 'CHSH':
            S = v*2*r
            if S <= 2: out[kind] = [0.0]*len(ns); continue
            cl = [chsh_tangent(S0) + (S0,) for S0 in np.linspace(max(2.0001, S - 0.25), S, 26)]
            cl = [(a, b, None) for a, b, _ in cl]
        else:
            cl = [c for c in certs(lab, kind) if v - window - 1e-12 <= c[2] <= v + 0.0051]
            if kind == 'cone':
                if not cl: out[kind] = None; continue
                cl += [c for c in certs(lab, 'tilted') if v - window - 1e-12 <= c[2] <= v + 0.0051]   # tilted CHSH lies in C_theta
            if not cl: out[kind] = None; continue
        res = []
        for n in ns:
            best = 0.0
            for K0, kap, _ in cl:
                if K0 + kap @ z <= best: continue
                b = eat_rate(K0, kap, z, n); best = max(best, b['rate'])
            res.append(best)
        out[kind] = res
    return out
if __name__ == '__main__':
    import sys
    T = asym_table(); json.dump(T, open('out/asym_table.json', 'w'), indent=1)
    L = ['| theta | v | CHSH (tight) | tilted CHSH | best exposing | full stats | gain tilted/CHSH | gain best/CHSH | gain full/CHSH | max-ent CHSH |', '|' + '---|'*10]
    f3 = lambda x: '-' if x is None else '%.4f' % x
    def gain(a, b): return '-' if (a is None or b is None or b <= 0) else ('%+.0f%%' % (100*(a/b - 1)))
    for d in T:
        L.append('| %s | %.3f | %s | %s | %s | %s | %s | %s | %s | %s |' % (d['label'], d['v'], f3(d['CHSH']), f3(d['tilted']), f3(d['cone']), f3(d['full']),
                 gain(d['tilted'], d['CHSH']), gain(d['cone'], d['CHSH']), gain(d['full'], d['CHSH']), f3(d['maxent_CHSH'])))
    open('out/table_asym.md', 'w').write('\n'.join(L) + '\n'); print('\n'.join(L))
    if 'noeat' not in sys.argv:
        ns = list(np.logspace(6, 12, 13))
        E = json.load(open('out/eat_table.json')) if os.path.exists('out/eat_table.json') else dict(ns=ns, rates={})
        E['rates'] = {k: v for k, v in E['rates'].items() if not k.endswith('|cone') and v is not None}
        for lab in ('27.76deg', 'pi/8', 'pi/4', '3pi/16', '0.22pi', 'pi/16', '0.1pi'):
            for v in (0.999, 0.994, 0.99, 0.985, 0.98):
                if not any(d['label'] == lab and abs(d['v'] - v) < 1e-12 for d in MAIN): continue
                kinds = [k for k in ('CHSH', 'tilted', 'cone', 'full') if '%s|%.3f|%s' % (lab, v, k) not in E['rates']]
                if lab == 'pi/4': kinds = [k for k in kinds if k in ('CHSH', 'full')]
                if not kinds: continue
                cur = eat_curves(lab, v, ns, kinds=kinds)
                for k, vals in cur.items():
                    if vals is not None: E['rates']['%s|%.3f|%s' % (lab, v, k)] = vals
                json.dump(E, open('out/eat_table.json', 'w'), indent=1)      # checkpoint
                print(lab, v, {k: None if vals is None else ['%.4f' % x for x in vals[::3]] for k, vals in cur.items()}, flush=True)
        L = ['| theta | v | n | CHSH | tilted | best exposing | full-stat dual expr. | best/CHSH |', '|' + '---|'*8]
        for key in sorted(set(k.rsplit('|', 1)[0] for k in E['rates'])):
            lab, v = key.split('|')
            for j, n in enumerate(ns):
                if j % 2: continue
                g = lambda k: E['rates'].get(key + '|' + k, [None]*len(ns))[j]
                ch, ti, co, fu = g('CHSH'), g('tilted'), g('cone'), g('full')
                if all(x is None for x in (ti, co, fu)): continue
                L.append('| %s | %s | %.0e | %s | %s | %s | %s | %s |' % (lab, v, n, f3(ch), f3(ti), f3(co), f3(fu), gain(max(x for x in (ti, co, fu) if x is not None), ch)))
        open('out/table_eat.md', 'w').write('\n'.join(L) + '\n')
