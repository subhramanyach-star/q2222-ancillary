import json, numpy as np
D = json.load(open('bench.json'))
def fmt_gap(x):
    if x is None: return '--'
    if abs(x) < 2e-6: return '$\\lesssim10^{-6}$'
    m, e = f"{x:.1e}".split('e'); return f"${m}\\times10^{{{int(e)}}}$"
def sci(x):
    if x < 1e-15: return '$<10^{-15}$'
    m, e = f"{x:.0e}".split('e'); return f"$10^{{{int(e)}}}$" if m=='1' else f"${m}\\times10^{{{int(e)}}}$"
fam_short = {'CHSH':'CHSH','tilted CHSH [AMP]':'tilted CHSH~\\cite{AMP}','doubly tilted CHSH [Gigena]':'doubly tilted~\\cite{Gigena2025}',
 'MK double-tilted':'MK double-tilted~\\cite{MikosKaniewski}','generalised Wolfe-Yelin [MK]':'gen.\\ Wolfe--Yelin~\\cite{MikosKaniewski}',
 'WBC I_phi':'$I_\\phi$~\\cite{WBC22}','WBC J_gamma':'$J_\\gamma$~\\cite{WBC22}','extended tilted CHSH [BSB, MPW]':'ext.\\ tilted~\\cite{CustomSOS,MPW}',
 'asymmetric CHSH [WAP]':'asym.\\ CHSH~\\cite{WAP2021}','BB octagon vertex g0':'BB vertex~\\cite{BBExtremal}'}
def lab(l):
    l = l.replace('alpha0','@0').replace('alpha1','@1').replace('alpha','\\alpha').replace('@0','\\alpha_0').replace('@1','\\alpha_1').replace('beta','\\beta').replace('theta','\\theta').replace('phi','\\phi').replace('gamma','\\gamma').replace(',',',\\ ')
    return '$'+l+'$' if l!='CHSH' and l!='g0' else ('--' if l=='CHSH' else '$g_0$')
rows=[]; times=[]; devs=[]
for d in D:
    ref = d['closed'] if d['closed'] is not None else d['qubit_lower']
    dev = abs(d['socp']-ref); devs.append(dev); times.append(d['socp_time_s'])
    src = 'closed' if d['closed'] is not None else 'qubit'
    loc = d['local']; 
    rows.append(f"{fam_short[d['family']]} & {lab(d['label'])} & {loc:.6f} & {d['socp']:.10f} & {sci(dev)}{'' if src=='closed' else ('$^\\ddagger$' if abs(d['socp']-loc)<1e-9 else '$^\\dagger$')} & {fmt_gap(d['npa_1+AB']-ref)} & {fmt_gap(d['npa_2']-ref)} & {fmt_gap(d['npa_3']-ref)} & {d['socp_time_s']:.2f}\\\\")
open('bench_table.tex','w').write('\n'.join(rows)+'\n')
print('\n'.join(rows)); print('max dev',max(devs),'times',min(times),max(times),np.median(times))
for d in D:
    if 'printed_sextic_largest_root' in d: print(d['label'], d['printed_sextic_largest_root'], d['socp'])
