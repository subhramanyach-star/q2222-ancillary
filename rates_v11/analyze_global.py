"""Tables + figure for the global (two-party) runs."""
import numpy as np, json, os, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from run_global_cfg import TH, VS
def load(fn): return [json.loads(l) for l in open(fn)] if os.path.exists(fn) else []
B = load('out/global/bff_global.jsonl'); Hm = load('out/global/hmin_global.jsonl')
def get(lab, v, proto, D=B, key='H_cert'):
    r = [d for d in D if d['label'] == lab and abs(d['v'] - v) < 1e-12 and d['proto'] == proto and d.get(key) is not None]
    return r[-1][key] if r else None
def cert(lab, v, proto):
    w = get(lab, v, 'witness'); f = get(lab, v, 'full')
    if proto == 'witness': return None if w is None else max(w, 0)
    return None if f is None and w is None else max(x for x in (w, f, 0) if x is not None)   # full >= witness (smaller feasible set)
labs = [l for l, _ in TH]; f4 = lambda x: '-' if x is None else '%.4f' % x
L = ['| v | ' + ' | '.join('%s witness / full' % l for l in labs) + ' |', '|' + '---|' * (len(labs) + 1)]
for v in sorted(VS, reverse=True):
    L.append('| %.3f | ' % v + ' | '.join('%s / %s' % (f4(cert(l, v, 'witness')), f4(cert(l, v, 'full'))) for l in labs) + ' |')
L += ['', 'Envelope (certified vN, bits/round): best tilted theta < pi/4 vs maximally entangled CHSH (pi/4)', '',
      '| v | max-ent CHSH witness+B2 | max-ent full | best tilted witness (theta) | best tilted full (theta) | Hmin max-ent witness | Hmin best tilted witness (theta) | Hmin max-ent full | Hmin best tilted full (theta) |', '|' + '---|' * 9]
env = []
for v in sorted(VS, reverse=True):
    def best(proto, D=None):
        c = [((cert(l, v, proto) if D is None else get(l, v, proto, Hm, 'Hmin')), l) for l in labs if l != 'pi/4']
        c = [x for x in c if x[0] is not None]; return max(c) if c else (None, None)
    bw, bf = best('witness'), best('full'); hw, hf = best('witness', 1), best('full', 1)
    rec = dict(v=v, me_w=cert('pi/4', v, 'witness'), me_f=cert('pi/4', v, 'full'), bw=bw, bf=bf,
               hme_w=get('pi/4', v, 'witness', Hm, 'Hmin'), hme_f=get('pi/4', v, 'full', Hm, 'Hmin'), hw=hw, hf=hf); env.append(rec)
    t = lambda x: '-' if x[0] is None else '%.4f (%s)' % x
    L.append('| %.3f | %s | %s | %s | %s | %s | %s | %s | %s |' % (v, f4(rec['me_w']), f4(rec['me_f']), t(bw), t(bf), f4(rec['hme_w']), t(hw), f4(rec['hme_f']), t(hf)))
open('out/global/table_global.md', 'w').write('\n'.join(L) + '\n'); json.dump(env, open('out/global/envelope_global.json', 'w'), indent=1)
print('\n'.join(L))
# figure
plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False})
fig, ax = plt.subplots(1, 2, figsize=(7.0, 2.7), gridspec_kw=dict(wspace=0.28)); cm = plt.cm.viridis
tex = {'pi/16': r'\pi/16', '0.1pi': r'0.1\pi', 'pi/8': r'\pi/8', '27.76deg': r'27.76^\circ', '3pi/16': r'3\pi/16', '0.22pi': r'0.22\pi', 'pi/4': r'\pi/4\ \mathrm{(CHSH)}'}
for k, proto in enumerate(('witness', 'full')):
    a = ax[k]
    order = ['pi/16', '0.1pi', 'pi/8', '27.76deg', '3pi/16', '0.22pi', 'pi/4']
    for i, l in enumerate(order):
        pts = sorted((v, cert(l, v, proto)) for v in VS if cert(l, v, proto) is not None)
        if not pts: continue
        st = dict(color='k', lw=1.8) if l == 'pi/4' else dict(color=cm(i / 6.5), lw=0.9)
        a.plot(*zip(*pts), '-', marker='.', ms=4, label=r'$%s$' % tex[l], **st)
    a.set_xlim(0.948, 1.0005); a.set_ylim(0, 2); a.set_xlabel('visibility $v$')
    a.set_title('exposing witness + $B_2$ row' if proto == 'witness' else 'full statistics (2x3)', fontsize=8.5)
ax[0].set_ylabel(r'certified $H(AB|X{=}x^*,Y{=}2,E)$ [bits]'); ax[0].legend(fontsize=6, loc='upper left', ncol=2, handlelength=1.5)
fig.savefig('fig_global_v11.pdf', bbox_inches='tight'); fig.savefig('fig_global_v11.png', dpi=200, bbox_inches='tight')
