import numpy as np, json, warnings
warnings.filterwarnings('ignore')
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from families import ptheta
from eat import chsh_analytic
import analyze as A
plt.rcParams.update({'font.size': 8.5, 'font.family': 'serif', 'mathtext.fontset': 'cm', 'axes.linewidth': 0.6,
                     'xtick.direction': 'in', 'ytick.direction': 'in', 'legend.frameon': False})
T = json.load(open('out/asym_table.json')); E = json.load(open('out/eat_table.json'))
col = dict(CHSH='#555555', tilted='#1f77b4', cone='#d62728', full='#2ca02c', maxent='#999999')
fig, ax = plt.subplots(2, 2, figsize=(7.0, 5.0), gridspec_kw=dict(wspace=0.28, hspace=0.42)); ax = ax.ravel()
ax = [ax[0], ax[1], ax[3], ax[2]]   # (a),(b) fixed theta; (c) envelope; (d) EAT
for k, lab in enumerate(['27.76deg', 'pi/8']):
    a = ax[k]; th = A.LABELS[lab]; p, c0, (s, C, r) = ptheta(th)
    rows = [d for d in T if d['label'] == lab]; v = np.array([d['v'] for d in rows])
    vv = np.linspace(0.85, 1, 400)
    a.plot(vv, [chsh_analytic(2*r*x) if 2*r*x > 2 else 0 for x in vv], color=col['CHSH'], lw=1.2, label='CHSH score (tight)')
    a.plot(vv, chsh_analytic(2*np.sqrt(2)*vv), color=col['maxent'], lw=0.8, ls=':', label=r'CHSH, max. entangled')
    def series(key):
        x = [d['v'] for d in rows if d[key] is not None]; y = [d[key] for d in rows if d[key] is not None]; return x, y
    x, y = series('tilted'); a.plot(x, y, 'o-', color=col['tilted'], ms=2.5, lw=1, label='tilted-CHSH score')
    x, y = series('cone'); a.plot(x, y, 's', color=col['cone'], ms=3.2, mfc='none', mew=0.9, label='best exposing expr.')
    x, y = series('full'); a.plot(x, y, '^-', color=col['full'], ms=2.5, lw=1, label='full statistics')
    a.axvline((C + r)/2, color=col['tilted'], lw=0.5, ls='--'); a.axvline(1/r, color=col['CHSH'], lw=0.5, ls='--')
    a.set_xlim(0.86, 1.0); a.set_ylim(0, 1.0); a.set_xlabel('visibility $v$')
    a.set_title(r'$\theta=%s$' % ({'27.76deg': '27.76^\\circ', 'pi/8': '\\pi/8'}[lab]), fontsize=9)
    if k == 0: a.set_ylabel(r'certified $H(A_1|X{=}1,E)$ [bits]')
    for vexp in (0.994, 0.985): a.axvline(vexp, color='k', lw=0.3, alpha=0.4)
ax[0].legend(loc='upper left', fontsize=6.3, handlelength=1.8)
a = ax[2]; ns = np.array(E['ns'])
for kind, lsty in [('CHSH', '-'), ('tilted', '-'), ('full', '-')]:
    for vop, ls in [(0.994, '-'), (0.985, '--')]:
        key = '27.76deg|%.3f|%s' % (vop, kind)
        if key in E['rates']: a.semilogx(ns, E['rates'][key], ls=ls, color=col[kind], lw=1, label=('%s' % kind) if vop == 0.994 else None)
key = 'pi/4|0.994|CHSH'
if key in E['rates']: a.semilogx(ns, E['rates'][key], ls='-', color='k', lw=1.6)
a.set_xlabel('rounds $n$'); a.set_ylabel(r'$H_{\min}^{\epsilon}/n$ (EAT) [bits/round]'); a.set_ylim(0, 1.32)
a.set_title(r'EAT finite-size rate, $\theta=27.76^\circ$ vs max-ent.', fontsize=8)
from matplotlib.lines import Line2D
a.legend(handles=[Line2D([], [], color='k', lw=1.6), Line2D([], [], color=col['CHSH'], lw=1), Line2D([], [], color=col['tilted'], lw=1), Line2D([], [], color=col['full'], lw=1),
                  Line2D([], [], color='k', lw=1, ls='-'), Line2D([], [], color='k', lw=1, ls='--')],
         labels=[r'max-ent. CHSH ($\pi/4$), $v=0.994$', 'CHSH score (same state)', r'tilted CHSH (= best in $\mathcal{C}_\theta$)', 'full-stat. dual expr.', '$v=0.994$', '$v=0.985$'], fontsize=5.6, loc='upper left', ncol=2, handlelength=1.6, columnspacing=0.8)
EN = json.load(open('out/envelope.json'))
a = ax[3]; cm = plt.cm.viridis
labs = ['pi/16', '0.1pi', 'pi/8', '27.76deg', '3pi/16', '0.22pi']
tex = {'pi/16': r'\pi/16', '0.1pi': r'0.1\pi', 'pi/8': r'\pi/8', '27.76deg': r'27.76^\circ', '3pi/16': r'3\pi/16', '0.22pi': r'0.22\pi'}
for i, lab in enumerate(labs):
    pts = [(e['v'], e['vN'][lab]['tilted']) for e in EN if lab in e['vN'] and e['vN'][lab].get('tilted') is not None]
    pts.sort(); a.plot(*zip(*pts), '-', marker='.', ms=3, lw=0.9, color=cm(i/len(labs)), label=r'$%s$' % tex[lab])
vv = np.linspace(0.85, 1, 400)
a.plot(vv, chsh_analytic(2*np.sqrt(2)*vv), 'k-', lw=1.6, label=r'$\pi/4$ (CHSH)')
pts = sorted((e['v'], e['Hmin_chsh']) for e in EN); a.plot(*zip(*pts), 'k:', lw=1, label=r'$\pi/4$, $H_{\min}$')
a.set_xlim(0.86, 1.0); a.set_ylim(0, 1.0); a.set_xlabel('visibility $v$'); a.set_ylabel(r'certified $H(A_1|X{=}1,E)$ [bits]')
a.set_title(r'tilted-CHSH score, all $\theta$ (white noise)', fontsize=8)
a.legend(fontsize=6, loc='upper left', ncol=2, handlelength=1.5)
for i, a in enumerate(ax): a.text(0.97 if i < 2 else 0.97, 0.04, '(%s)' % 'abdc'[i], transform=a.transAxes, va='bottom', ha='right', fontsize=9)
fig.savefig('fig_rates_v11.pdf', bbox_inches='tight'); fig.savefig('fig_rates_v11.png', dpi=200, bbox_inches='tight')
print('saved')
