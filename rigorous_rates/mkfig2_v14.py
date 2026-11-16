# v14: Fig. 2 drawn from the certified two-sided intervals of rigorous_rates/out/verified.json (written by the
# standalone verifier verify_all.py). Every plotted point is the midpoint of an interval of width < 2e-6 bits;
# lines join the certified points. Diamonds: certified supremum over the whole hexagon C_theta.
import json, os, numpy as np, matplotlib, re
HERE=os.path.dirname(os.path.abspath(__file__))
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
d={r['id']:r for r in json.load(open(os.path.join(HERE,'out','verified.json')))}
def series(key):
    pts=[]
    for k,r in d.items():
        m=re.fullmatch(re.escape(key)+r'_v([0-9.]+)',k)
        if m:
            assert r['attack_ok'] and r['H_hi']-r['H_lo']<2e-6, k
            pts.append((float(m.group(1)),0.5*(r['H_lo']+r['H_hi'])))
    pts.sort(reverse=True); v=np.array([p[0] for p in pts]); return 1-v,np.array([p[1] for p in pts])
fig,ax=plt.subplots(2,1,figsize=(3.4,3.6),sharex=True)
a=ax[0]
x,y=series('dense_thpi_8_full'); a.plot(x,y,'k-',lw=1.2,label='full statistics')
x,y=series('dense_thpi_8_tilted'); a.plot(x,y,'b-',lw=1.2,label='tilted CHSH')
x,y=series('dense_thpi_8_chsh'); a.plot(x,y,'r--',lw=1.2,label='CHSH')
x,y=series('thpi_8_conesup'); a.plot(x,y,'gD',ms=3.5,label='best in hexagon (supremum)')
a.axvline(1-0.96593,color='b',ls=':',lw=0.8)
a.set_xscale('log'); a.set_ylabel('min-entropy [bits/round]')
a.set_title(r'(a) tilted target $\theta=\pi/8$, outcome $A_1$',fontsize=8); a.legend(fontsize=6.5,frameon=False,loc='upper right'); a.set_ylim(0,1.0)
a=ax[1]
x,y=series('dense_Idpi_6_full'); a.plot(x,y,'k-',lw=1.2,label='full statistics')
x,y=series('dense_Idpi_6_centre'); a.plot(x,y,'b-',lw=1.2,label=r'WBC expression $I_\phi$ (centre)')
x,y=series('dense_Idpi_6_chsh'); a.plot(x,y,'r--',lw=1.2,label='CHSH')
x,y=series('Idpi_6_AB_half'); a.plot(x,y,'gs',ms=3.5,mfc='none',label=r'half-way to vertex $W_1$')
a.axvline(1-0.96225,color='b',ls=':',lw=0.8)
a.set_xscale('log'); a.set_xlabel(r'noise $1-v$'); a.set_ylabel('min-entropy [bits/round]')
a.set_title(r'(b) WBC target $\phi=\pi/6$, outcomes $(A_0,B_0)$',fontsize=8); a.legend(fontsize=6.5,frameon=False,loc='upper right'); a.set_ylim(0,2.0)
fig.tight_layout(pad=0.3); fig.savefig(os.path.join(HERE,'fig2_rates_v14.pdf')); print('ok')
