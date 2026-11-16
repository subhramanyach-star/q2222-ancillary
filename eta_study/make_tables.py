"""Markdown tables for REPORT.md from verified certificates (out/verified_eta.json), optimisation records and vN runs."""
import json, os, numpy as np
from fractions import Fraction as F
from core import load_opt, Hvn_chsh
VER=json.load(open('out/verified_eta.json')) if os.path.exists('out/verified_eta.json') else []
def ver(V,guess=('A',1)): return {d['eta']:d for d in VER if abs(d['v']-V)<1e-9 and tuple(d['guess'])==tuple(guess)}
def f6(x): return '%.6f'%x if x is not None else '–'
def iv(p): return '–' if p is None else ('%.6f'%p[0] if p[1] is not None and abs(p[1]-p[0])<5e-7 else '[%.6f, %s]'%(p[0],'%.6f'%p[1] if p[1] is not None else '?'))
def hmin_table(V,guess=('A',1),suf=''):
    R=load_opt(V,suf); C=ver(V,guess)
    lines=['| η | S_max(η) ≤ | CHSH, \\|Φ+⟩ | **every CHSH protocol ≤** | tilted α⟨A₀⟩+CHSH (proved) | best in C_θ′ (proved) | θ′/π | full stats (numerical) | gain tilted | gain cone | separation |',
           '|---|---|---|---|---|---|---|---|---|---|---|']
    for eta in sorted(C,reverse=True):
        c=C[eta]
        if 'tilted' not in c and 'cone' not in c: continue
        mx=R.get((eta,'chsh_maxent'),{}).get('H'); fu=R.get((eta,'full'),{}).get('H')
        cone=R.get((eta,'cone'),{}); tp=cone.get('thetap')
        gt=100*(c['tilted'][0]/c['chsh_H_up']-1) if 'tilted' in c and c['chsh_H_up']>0 else None
        gc=100*(c['cone'][0]/c['chsh_H_up']-1) if 'cone' in c and c['chsh_H_up']>0 else None
        sep=[]
        if 'tilted' in c and c['tilted'][0]>c['chsh_H_up']: sep.append('tilted')
        if 'cone' in c and c['cone'][0]>c['chsh_H_up']: sep.append('cone')
        lines.append('| %.3f | %.6f | %s | %s | %s | %s | %s | %s | %s | %s | %s |'%(eta,c['S_up'],('0' if mx is not None and mx<=1e-9 else f6(mx)),f6(c['chsh_H_up']),
            iv(c.get('tilted')),iv(c.get('cone')),'%.4f'%(tp/np.pi) if tp else '–',f6(fu),'%+.1f%%'%gt if gt is not None else '–','%+.1f%%'%gc if gc is not None else '–',
            ('**proved** ('+', '.join(sep)+')') if sep else 'not proved'))
    return '\n'.join(lines)
def vn_table(V):
    fn='out/vn_v%s.jsonl'%V
    if not os.path.exists(fn): return '(no vN data)'
    B={}
    for l in open(fn): d=json.loads(l); B[(d['eta'],d['fam'])]=d
    C=ver(V)
    lines=['| η | CHSH, best device | **every CHSH protocol ≤** | tilted (BFF, exact-PSD cert.) | α | full stats (BFF cert.) | separation | EAT n=1e8: CHSH / tilted / full | EAT n=1e10: CHSH / tilted / full |','|---|---|---|---|---|---|---|---|---|']
    for eta in sorted(set(e for e,f in B),reverse=True):
        ch=B.get((eta,'chsh')); ti=B.get((eta,'tilted')); fu=B.get((eta,'full'))
        if not ch: continue
        up=Hvn_chsh(C[eta]['S_up']) if eta in C else ch.get('H_upper_every_chsh')
        sep='–'
        if ti and up is not None: sep='**proved**' if ti['H_cert']>up else 'not proved'
        e=lambda r,k: ('%.4f'%r['eat'][k] if r['eat'][k]>0 else 'n/a†') if r else '–'
        lines.append('| %.3f | %.6f | %s | %s | %s | %s | %s | %s / %s / %s | %s / %s / %s |'%(eta,ch['H_cert'],f6(up),f6(ti['H_cert'] if ti else None),'%.3f'%ti['alpha'] if ti else '–',
              f6(fu['H_cert'] if fu else None),sep,e(ch,0),e(ti,0),e(fu,0),e(ch,1),e(ti,1),e(fu,1)))
    return '\n'.join(lines)
if __name__=='__main__':
    out=['## Table 1: H_min(A1|X=1,E), v = 1\n',hmin_table(1.0),'\n\n## Table 2: H_min(A1|X=1,E), v = 0.99\n',hmin_table(0.99),
         '\n\n## Table 3: von Neumann H(A1|X=1,E) and EAT rates, v = 1\n',vn_table(1.0),'\n\n## Table 4: von Neumann, v = 0.99\n',vn_table(0.99),
         '\n\n## Table 5: two-party H_min(A1B1|X=1,Y=1,E), v = 1\n',hmin_table(1.0,('AB',1,1),'_AB11')]
    open('out/tables.md','w').write('\n'.join(out)); print('\n'.join(out)[:3000])
