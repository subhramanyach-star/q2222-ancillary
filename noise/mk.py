"""Tables from out/opt_*.jsonl (float optima) and out/cert_*.jsonl (certificates) -> out/tables.md"""
import json, os
def load(fn,key):
    R={}
    if not os.path.exists(fn): return R
    for l in open(fn):
        d=json.loads(l); k=key(d)
        if fn.find('opt_')>=0:
            if k not in R or d['H']>R[k]['H']: R[k]=d
        else: R[k]=d
    return R
out=[]
for m,lab in (('white','white noise, visibility v (F=(1+3v)/4)'),('deph','dephasing, coherence p (F=1-(1-p)sin^2(2th)/2)')):
    O=load('out/opt_%s.jsonl'%m,lambda d:(d['eta'],d['q'],d['fam'])); C=load('out/cert_%s.jsonl'%m,lambda d:(d['eta'],d['q']))
    out.append('\n### %s\n'%lab)
    out.append('| eta | %s | best CHSH H_min | CHSH cert. <= | tilted >= | C_theta\' >= | gain tilted | gain C_theta\' | status |'%('v' if m=='white' else 'p'))
    out.append('|---|---|---|---|---|---|---|---|---|')
    for (eta,q) in sorted(set(k[:2] for k in O),key=lambda t:(t[0],-t[1])):
        c=O.get((eta,q,'chsh')); t=O.get((eta,q,'tilted')); o=O.get((eta,q,'cone'))
        if not (c and t and o): continue
        ce=C.get((eta,q)); Hc=c['H']
        if ce and 'tilted' in ce and 'cone' in ce:
            ht,hk=ce['tilted'][0],ce['cone'][0]; up=ce['H_chsh_up_cert']
            st=[]
            for nm,h in (('tilted',ht),('cone',hk)):
                if h>up: st.append(nm)
            status=('proved: '+', '.join(st)) if st else 'not proved'
            if m=='deph': status+=' (vs loose bound)'
            cs='%.6f'%up; ts='%.6f'%ht; ks='%.6f'%hk
        else:
            ht,hk=t['H'],o['H']; status='float only'; cs='–'; ts='%.6f*'%ht; ks='%.6f*'%hk
        g=lambda h: ('%+.1f%%'%(100*(h/Hc-1))) if Hc>1e-12 else ('CHSH=0' if h>1e-9 else '–')
        out.append('| %.2f | %.3f | %.6f | %s | %s | %s | %s | %s | %s |'%(eta,q,Hc,cs,ts,ks,g(ht),g(hk),status))
open('out/tables.md','w').write('\n'.join(out)+'\n'); print('\n'.join(out))
# ---- experiment-matched points (model wd: white v and dephasing p)
EXPS={(0.7515,(0.999,0.997)):'Shalm et al. PRL 2015 (NIST)',(0.755,(0.999,0.997)):'Bierhorst et al. Nature 2018 (NIST; vis. assumed = 2015 source)',
      (0.763,(0.999,0.997)):'Shalm et al. Nat. Phys. 2021 (NIST; vis. assumed)',(0.774,(0.99,1.0)):'Giustina et al. PRL 2015 (Vienna)',
      (0.7865,(0.995,0.99)):'Liu et al. Nature 2018 (USTC DIQRNG)',(0.841,(0.994,0.991)):'Liu et al. Nat. Phys. 2021 (USTC expansion)',
      (0.8749,(0.9936,1.0)):'Liu et al. PRL 2022 (USTC DIQKD), F as white noise',(0.8749,(1.0,0.9768)):'Liu et al. PRL 2022, F as dephasing'}
O=load('out/opt_wd.jsonl',lambda d:(d['eta'],tuple(d['q']),d['fam'])); C=load('out/cert_wd.jsonl',lambda d:(d['eta'],tuple(d['q'])))
out=['\n### experiment-matched points (white v and dephasing p)\n','| experiment | eta | v | p | best CHSH | CHSH cert. <= | tilted >= | C_theta\' >= | gain tilted | gain C_theta\' | status |','|---|---|---|---|---|---|---|---|---|---|---|']
for (eta,q),name in sorted(EXPS.items()):
    c=O.get((eta,q,'chsh')); t=O.get((eta,q,'tilted')); o=O.get((eta,q,'cone'))
    if not (c and t and o): continue
    ce=C.get((eta,q)); Hc=c['H']
    if ce and 'tilted' in ce and 'cone' in ce:
        ht,hk,up=ce['tilted'][0],ce['cone'][0],ce['H_chsh_up_cert']; st=[n for n,h in (('tilted',ht),('cone',hk)) if h>up]
        status=('proved: '+', '.join(st)) if st else 'not proved (CHSH side numerical)'
        cs,ts,ks='%.6f'%up,'%.6f'%ht,'%.6f'%hk
    else: ht,hk=t['H'],o['H']; status='float only'; cs='–'; ts,ks='%.6f*'%ht,'%.6f*'%hk
    out.append('| %s | %.4f | %.4f | %.4f | %.6f | %s | %s | %s | %+.1f%% | %+.1f%% | %s |'%(name,eta,q[0],q[1],Hc,cs,ts,ks,100*(ht/Hc-1),100*(hk/Hc-1),status))
open('out/tables.md','a').write('\n'.join(out)+'\n'); print('\n'.join(out))
