import json, os, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
d=json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'lossycone','t4_summary.json')))
P=d['points']; eta=[p['eta'] for p in P]; Hc=[p['H_lo'] for p in P]; cap=[p['H_cap'] for p in P]
g=[100*p['gain'] for p in P]
# WAP (certified, own best device) and polygon C_theta' (certified), wapcone/REPORT.md
E6=[0.75,0.80,0.85,0.90,0.95,0.98]
chsh6=[0.009013,0.033444,0.081930,0.168312,0.330204,0.525646]
wap=[0.009041,0.033995,0.085560,0.182071,0.366937,0.580212]
pol=[0.009465,0.036046,0.089927,0.186183,0.361305,0.558502]
gw=[100*(w/c-1) for w,c in zip(wap,chsh6)]
gp=[100*(w/c-1) for w,c in zip(pol,chsh6)]
# interval guarantee
iv=[(0.75,0.80,5.3),(0.80,0.85,8.5),(0.85,0.90,10.7),(0.90,0.95,12.2),(0.95,0.98,8.0)]
plt.rcParams.update({'font.size':11,'axes.labelsize':12,'legend.fontsize':9.5,'xtick.labelsize':10.5,'ytick.labelsize':10.5})
fig,ax=plt.subplots(2,1,figsize=(3.4,4.6),sharex=True,gridspec_kw={'height_ratios':[1,1.15],'hspace':0.08})
a=ax[0]
a.semilogy(eta,Hc,'-',color='C3',lw=2,label='our witness')
a.semilogy(E6,wap,'s',color='C0',ms=5,label='asymmetric CHSH')
a.semilogy(eta,cap,'--',color='k',lw=1.4,label='every CHSH protocol')
a.set_ylabel(r'$H_{\min}$ (bits/round)')
a.legend(loc='upper left',frameon=False,handlelength=1.6)
a.text(0.9,0.06,'(a)',transform=a.transAxes,fontsize=12)
a.set_ylim(0.007,3)
b=ax[1]
b.plot(eta,g,'-',color='C3',lw=2,label='our witness')
b.plot(E6,gw,'s',color='C0',ms=5,label='asym. CHSH')
b.plot(E6,gp,'^',color='C2',ms=5,label=r'polygon $\mathcal{C}_{\theta^\prime}$')
for i,(l,h,v) in enumerate(iv):
    b.plot([l,h],[v,v],color='C3',lw=1.2,ls=':',label='every $\\eta$ (proved)' if i==0 else None)
b.axhline(0,color='k',lw=1.0,ls='--')
b.set_ylabel('gain over CHSH (%)'); b.set_xlabel(r'detection efficiency $\eta$')
b.set_ylim(-1,27); b.set_xlim(0.745,0.985)
b.legend(loc='upper left',frameon=False,handlelength=1.3,ncol=2,borderaxespad=0.2,columnspacing=0.6,fontsize=8.5)
b.text(0.9,0.06,'(b)',transform=b.transAxes,fontsize=12)
fig.savefig('../fig_gain_v16.pdf',bbox_inches='tight')
for e,a1,b1,c1 in zip(E6,gw,gp,chsh6): print(e,round(a1,2),round(b1,2))
