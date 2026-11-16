# v14: single-column version of eta_study/figs/fig_eta_prl.pdf (same data: eta_study/out/*.jsonl, out/verified_eta.json).
import os, sys, numpy as np
os.chdir(os.path.dirname(os.path.abspath(__file__))); sys.path.insert(0,'.')
import make_figs as M, matplotlib.pyplot as plt
plt.rcParams.update({'font.size':8,'font.family':'serif','mathtext.fontset':'cm'})
R=M.load(1.0); C1=M.certs(1.0); C2=M.certs(0.99); STY=M.STY
fig,ax=plt.subplots(2,1,figsize=(3.4,4.0),sharex=True)
for fam in ('chsh_maxent','chsh_opt','tilted','cone','full'):
    e,h=M.curve(R,fam); k=h>0
    ax[0].plot(e[k],h[k],STY[fam][2],color=STY[fam][1],label=STY[fam][0],lw=1.0)
es=sorted(C1); ax[0].plot(es,[C1[x]['chsh_H_up'] for x in es],'k_',ms=5,mew=0.9,label='upper bound, every CHSH protocol')
for fam,mk in (('tilted','o'),('cone','s')):
    pts=[(x,C1[x][fam][0]) for x in es if fam in C1[x]]
    ax[0].plot(*zip(*pts),mk,mfc='none',color=STY[fam][1],ms=2.6,mew=0.7)
for a in ax:
    for x0 in (2/3,2/(1+np.sqrt(2))): a.axvline(x0,color='0.75',lw=0.6)
ax[0].set_xlim(0.66,1.0); ax[0].set_yscale('log'); ax[0].set_ylim(1e-4,1.05)
ax[0].set_ylabel('$H_{\\min}(A_1|X{=}1,E)$ [bits]')
ax[0].legend(loc='lower right',frameon=False,fontsize=5.8,handlelength=2.0)
for C,ls,lab in ((C1,'-','$v=1$'),(C2,'--','$v=0.99$')):
    for fam,mk in (('tilted','o'),('cone','s')):
        pts=[(x,100*(C[x][fam][0]/C[x]['chsh_H_up']-1)) for x in sorted(C) if fam in C[x] and C[x][fam][0]>C[x]['chsh_H_up']]
        ax[1].plot(*zip(*pts),ls,marker=mk,mfc='none',ms=2.6,mew=0.7,lw=0.9,color=STY[fam][1],label=('tilted CHSH' if fam=='tilted' else 'best in $C_{\\theta}$')+', '+lab)
ax[1].axhline(0,color='k',lw=0.6)
ax[1].set_xlabel('detector efficiency $\\eta$'); ax[1].set_ylabel('proved gain over CHSH [%]')
ax[1].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.10,1.0),ncol=2,fontsize=5.6,columnspacing=0.6,handlelength=1.8); ax[1].set_ylim(-0.5,17.5)
for a,t in zip(ax,'ab'): a.text(0.02,0.95,'(%s)'%t,transform=a.transAxes,va='top')
fig.tight_layout(pad=0.3); fig.savefig('figs/fig_eta_v14.pdf'); print('ok')
