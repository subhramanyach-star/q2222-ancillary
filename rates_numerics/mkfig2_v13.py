# Fig. 2 of prl_v13 (NPA level-3 numerics; floating point, not a proof). Reads dense.json, rates_th1.json, rates_Id2.json; writes fig2_rates_v13.pdf.
# v13: legend 'best test found in hexagon' -> 'best found in hexagon'. v11: same data as v9 fig3_rates_v9.pdf (NPA level-3 numerics), stacked single-column layout, plain labels.
import json, numpy as np, matplotlib, os
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
D=''
d=json.load(open(D+'dense.json')); v=np.array(d['v']); x=1-v
t1=json.load(open(D+'rates_th1.json')); I2=json.load(open(D+'rates_Id2.json'))
vs=np.array(t1['v']); xs=1-vs
fig,ax=plt.subplots(2,1,figsize=(3.4,3.6),sharex=True)
a=ax[0]; r=d['th1']
a.plot(x,r['full'],'k-',lw=1.2,label='full statistics')
a.plot(x,r['g_theta'],'b-',lw=1.2,label='tilted CHSH')
a.plot(x,r['CHSH'],'r--',lw=1.2,label='CHSH')
a.plot(xs,t1['H_min(A1|E)']['best in hexagon (search)'],'gD',ms=3.5,label='best found in hexagon')
a.axvline(1-0.96593,color='b',ls=':',lw=0.8)
a.set_xscale('log'); a.set_ylabel('min-entropy [bits/round]')
a.set_title(r'(a) tilted target $\theta=\pi/8$, outcome $A_1$',fontsize=8); a.legend(fontsize=6.5,frameon=False,loc='upper right'); a.set_ylim(0,1.0)
a=ax[1]; r=d['Id2']
a.plot(x,r['full'],'k-',lw=1.2,label='full statistics')
a.plot(x,r['I_delta'],'b-',lw=1.2,label=r'WBC expression $I_\phi$ (centre)')
a.plot(x,r['CHSH'],'r--',lw=1.2,label='CHSH')
a.plot(xs,I2['H_min(AB|00,E)']['half-way to W1'],'gs',ms=3.5,mfc='none',label=r'half-way to vertex $W_1$')
a.axvline(1-0.96225,color='b',ls=':',lw=0.8)
a.set_xscale('log'); a.set_xlabel(r'noise $1-v$'); a.set_ylabel('min-entropy [bits/round]')
a.set_title(r'(b) WBC target $\phi=\pi/6$, outcomes $(A_0,B_0)$',fontsize=8); a.legend(fontsize=6.5,frameon=False,loc='upper right'); a.set_ylim(0,2.0)
fig.tight_layout(pad=0.3); fig.savefig('fig2_rates_v13.pdf'); print('ok')
