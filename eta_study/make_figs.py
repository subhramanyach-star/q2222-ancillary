"""Figures for the eta study (reads out/*.jsonl; no computation)."""
import json, glob, os, numpy as np, matplotlib
matplotlib.use('Agg'); import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':9,'axes.labelsize':9,'legend.fontsize':7.5,'font.family':'serif','mathtext.fontset':'cm'})
from core import load_opt as load
def certs(V,guess=('A',1)):
    C={}
    if os.path.exists('out/verified_eta.json'):
        for d in json.load(open('out/verified_eta.json')):
            if abs(d['v']-V)<1e-9 and tuple(d['guess'])==tuple(guess): C[d['eta']]=d
    return C
STY={'chsh_maxent':('CHSH, $|\\Phi^+\\rangle$','0.55',':'),'chsh_opt':('CHSH, best state','k','--'),
     'tilted':('tilted CHSH $\\alpha\\langle A_0\\rangle+$CHSH','C0','-'),'cone':('best expression in $C_{\\theta}$','C3','-'),
     'full':('full statistics','C2','-.')}
def curve(R,fam):
    e=sorted(x for x,f in R if f==fam); return np.array(e),np.array([max(R[(x,fam)]['H'],0) for x in e])
def fig_hmin(V,fn,suf='',guess=('A',1),ylab='$H_{\\min}(A_1|X{=}1,E)$ [bits]'):
    R=load(V,suf); C=certs(V,guess)
    fig,ax=plt.subplots(1,2,figsize=(7.0,2.7))
    for fam in ('chsh_maxent','chsh_opt','tilted','cone','full'):
        e,h=curve(R,fam)
        if len(e): ax[0].plot(e,h,STY[fam][2],color=STY[fam][1],label=STY[fam][0],lw=1.2)
    if C:
        es=sorted(C); ax[0].plot(es,[C[x]['chsh_H_up'] for x in es],'kx',ms=4,label='proved upper bound, every CHSH protocol')
        for fam,mk in (('tilted','o'),('cone','s')):
            pts=[(x,C[x][fam][0]) for x in es if fam in C[x]]
            if pts: ax[0].plot(*zip(*pts),mk,mfc='none',color=STY[fam][1],ms=3.5,label='proved lower bound (%s)'%fam)
    ax[0].axvline(2/3,color='0.7',lw=0.6); ax[0].axvline(2*np.sqrt(2)-2 if False else 2/(1+np.sqrt(2)),color='0.7',lw=0.6)
    ax[0].set_xlabel('detector efficiency $\\eta$'); ax[0].set_ylabel(ylab); ax[0].set_yscale('log'); ax[0].set_ylim(1e-4,1.05)
    ax[0].legend(loc='lower right',frameon=False,fontsize=6.3,handlelength=2.2)
    eb,hb=curve(R,'chsh_opt')
    for fam in ('tilted','cone','full'):
        e,h=curve(R,fam); m=[(x,y/dict(zip(eb,hb))[x]) for x,y in zip(e,h) if x in dict(zip(eb,hb)) and dict(zip(eb,hb))[x]>1e-3]
        if m: ax[1].plot(*zip(*[(x,100*(r-1)) for x,r in m]),STY[fam][2],color=STY[fam][1],lw=1.2,label=STY[fam][0])
    if C:
        for fam,mk in (('tilted','o'),('cone','s')):
            pts=[(x,100*(C[x][fam][0]/C[x]['chsh_H_up']-1)) for x in sorted(C) if fam in C[x] and C[x]['chsh_H_up']>1e-6]
            if pts: ax[1].plot(*zip(*pts),mk,mfc='none',color=STY[fam][1],ms=3.5)
    ax[1].axhline(0,color='k',lw=0.6); ax[1].set_xlabel('detector efficiency $\\eta$'); ax[1].set_ylabel('gain over best CHSH protocol [%]')
    ax[1].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.0,0.9),fontsize=6.5); ax[1].set_xlim(0.66,1.0)
    for a,t in zip(ax,'ab'): a.text(0.02,0.95,'(%s)'%t,transform=a.transAxes,va='top')
    fig.tight_layout(); fig.savefig('figs/%s.pdf'%fn); fig.savefig('figs/%s.png'%fn,dpi=200); plt.close(fig)
def fig_vn(V,fn):
    B={}
    if not os.path.exists('out/vn_v%s.jsonl'%V): return
    for l in open('out/vn_v%s.jsonl'%V): d=json.loads(l); B[(d['eta'],d['fam'])]=d
    fig,ax=plt.subplots(1,2,figsize=(7.0,2.7))
    for fam,col,ls,lab in (('chsh','k','--','CHSH, best state'),('tilted','C0','-','tilted CHSH'),('full','C2','-.','full statistics')):
        e=sorted(x for x,f in B if f==fam)
        if e: ax[0].plot(e,[max(B[(x,fam)]['H_cert'],1e-6) for x in e],ls,color=col,marker='.',ms=3,lw=1.1,label=lab)
    ax[0].set_ylabel('$H(A_1|X{=}1,E)$ [bits] (certified)'); ax[0].set_yscale('log'); ax[0].set_ylim(1e-2,1.05)
    ax[0].legend(frameon=False,loc='lower right')
    e=sorted(x for x,f in B if f=='tilted' and (x,'chsh') in B)
    up=lambda x: B[(x,'chsh')].get('H_upper_every_chsh') or B[(x,'chsh')]['H_cert']
    ax[1].plot(e,[100*(B[(x,'tilted')]['H_cert']/up(x)-1) for x in e],'-o',color='C0',ms=3,mfc='none',lw=1.1,label='asymptotic (vs every CHSH protocol)')
    for k,ls,lab in ((1,'--','EAT, $n=10^{10}$ (vs best CHSH)'),(0,':','EAT, $n=10^{8}$ (vs best CHSH)')):
        ee=[x for x in e if B[(x,'chsh')]['eat'][k]>0 and B[(x,'tilted')]['eat'][k]>0]
        ax[1].plot(ee,[100*(B[(x,'tilted')]['eat'][k]/B[(x,'chsh')]['eat'][k]-1) for x in ee],ls,marker='.',ms=3,color='C0',lw=1.0,label=lab)
    ax[1].axhline(0,color='k',lw=0.6); ax[1].set_ylabel('gain of tilted CHSH [%]'); ax[1].legend(frameon=False,loc='lower center',fontsize=6.5)
    for a,t in zip(ax,'ab'):
        a.set_xlabel('detector efficiency $\\eta$'); a.text(0.02,0.95,'(%s)'%t,transform=a.transAxes,va='top')
    fig.tight_layout(); fig.savefig('figs/%s.pdf'%fn); fig.savefig('figs/%s.png'%fn,dpi=200); plt.close(fig)
def fig_prl(fn='fig_eta_prl'):
    """PRL figure: (a) H_min vs eta at v=1 (curves: numerical optimum per family; markers: proved bounds);
    (b) proved gain of the tilted and C_theta witnesses over every CHSH protocol, v=1 and v=0.99."""
    R=load(1.0); C1=certs(1.0); C2=certs(0.99)
    fig,ax=plt.subplots(1,2,figsize=(7.0,2.6))
    for fam in ('chsh_maxent','chsh_opt','tilted','cone','full'):
        e,h=curve(R,fam); k=h>0
        ax[0].plot(e[k],h[k],STY[fam][2],color=STY[fam][1],label=STY[fam][0],lw=1.1)
    es=sorted(C1); ax[0].plot(es,[C1[x]['chsh_H_up'] for x in es],'k_',ms=6,mew=1.0,label='upper bound, every CHSH protocol')
    for fam,mk in (('tilted','o'),('cone','s')):
        pts=[(x,C1[x][fam][0]) for x in es if fam in C1[x]]
        ax[0].plot(*zip(*pts),mk,mfc='none',color=STY[fam][1],ms=3,mew=0.8)
    for x0 in (2/3,2/(1+np.sqrt(2))): ax[0].axvline(x0,color='0.75',lw=0.6)
    ax[0].set_xlim(0.66,1.0); ax[0].set_yscale('log'); ax[0].set_ylim(1e-4,1.05)
    ax[0].set_xlabel('detector efficiency $\\eta$'); ax[0].set_ylabel('$H_{\\min}(A_1|X{=}1,E)$ [bits]')
    ax[0].legend(loc='lower right',frameon=False,fontsize=6.2,handlelength=2.0)
    for C,ls,lab in ((C1,'-','$v=1$'),(C2,'--','$v=0.99$')):
        for fam,mk in (('tilted','o'),('cone','s')):
            pts=[(x,100*(C[x][fam][0]/C[x]['chsh_H_up']-1)) for x in sorted(C) if fam in C[x] and C[x][fam][0]>C[x]['chsh_H_up']]
            if pts: ax[1].plot(*zip(*pts),ls,marker=mk,mfc='none',ms=3,mew=0.8,lw=0.9,color=STY[fam][1],label=('tilted CHSH' if fam=='tilted' else 'best in $C_{\\theta}$')+', '+lab)
    ax[1].axhline(0,color='k',lw=0.6); ax[1].set_xlim(0.66,1.0)
    ax[1].set_xlabel('detector efficiency $\\eta$'); ax[1].set_ylabel('proved gain over CHSH [%]')
    ax[1].legend(frameon=False,loc='upper left',bbox_to_anchor=(0.07,1.0),ncol=2,fontsize=6.3,columnspacing=0.8); ax[1].set_ylim(-0.5,14.5)
    for a,t in zip(ax,'ab'): a.text(0.02,0.95,'(%s)'%t,transform=a.transAxes,va='top')
    fig.tight_layout(); fig.savefig('figs/%s.pdf'%fn); fig.savefig('figs/%s.png'%fn,dpi=200); plt.close(fig)
if __name__=='__main__':
    fig_prl()
    fig_hmin(1.0,'fig_eta_hmin_v1'); fig_hmin(0.99,'fig_eta_hmin_v099')
    fig_vn(1.0,'fig_eta_vn_v1'); fig_vn(0.99,'fig_eta_vn_v099')
    if glob.glob('out/opt_AB11_v1.0*.jsonl'): fig_hmin(1.0,'fig_eta_hmin_AB_v1','_AB11',('AB',1,1),'$H_{\\min}(A_1B_1|11,E)$ [bits]')
    print('ok')
