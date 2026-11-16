import numpy as np, json, warnings
warnings.filterwarnings('ignore')
from npa import pguess
from families import *
LEV='3'
VS=[0.999,0.995,0.99,0.985,0.98,0.97,0.96,0.95,0.94,0.92,0.90,0.88,0.86,0.84]
H=lambda pg:-np.log2(pg)
out={'v':VS}
th=np.pi/8
P,order,V=hexagon(th); gt=gtilt(th); gch=g_chsh/(g_chsh@P)
res={}
for nm,g in [('g_theta',gt),('CHSH',gch)]:
    res[nm]=[H(pguess(LEV,('A',1),g=g,omega=v*(g@P))[0]) for v in VS]; print(nm,res[nm],flush=True)
res['full']=[H(pguess(LEV,('A',1),pobs=v*P)[0]) for v in VS]; print('full',res['full'],flush=True)
out['th1']=res
d=np.pi/6
P,g,Y,W1,W2=Idelta(d); gch=g_chsh/(g_chsh@P)
res={}
for nm,gg in [('I_delta',g),('CHSH',gch)]:
    res[nm]=[H(pguess(LEV,('AB',0,0),g=gg,omega=v*(gg@P))[0]) for v in VS]; print(nm,res[nm],flush=True)
res['full']=[H(pguess(LEV,('AB',0,0),pobs=v*P)[0]) for v in VS]; print('full',res['full'],flush=True)
out['Id2']=res
json.dump(out,open('dense.json','w'),indent=1)
print('DONE')
