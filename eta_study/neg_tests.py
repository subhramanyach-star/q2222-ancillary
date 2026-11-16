"""Negative controls for verify_eta.py: tampered certificates must be rejected."""
import json, copy, verify_eta as V
from fractions import Fraction as F
d=json.load(open('certs/eta0.85_v1.0.json'))
def run(name,mod):
    x=copy.deepcopy(d); mod(x)
    try: r=V.verify(x); print('%-40s NOT CAUGHT  %s'%(name,{k:r[k] for k in ('chsh_H_up','tilted') if k in r}))
    except Exception as e: print('%-40s caught: %s'%(name,str(e)[:60] or type(e).__name__))
run('CHSH NPA certificate entry perturbed',lambda x: x['chsh']['cert']['Z'][1].__setitem__(1,str(F(x['chsh']['cert']['Z'][1][1])-F(1,10**6))))
run('S_up lowered',lambda x: x['chsh'].__setitem__('S_up',str(F(x['chsh']['S_up'])-F(1,1000))))
run('CHSH attack weights changed',lambda x: x['chsh']['attack']['weights'].__setitem__(0,str(F(x['chsh']['attack']['weights'][0])+F(1,100))))
run('tilted threshold raised (omega)',lambda x: x['tilted'].__setitem__('omega',str(F(x['tilted']['omega'])+F(1,100))))
run('tilted mu changed',lambda x: x['tilted']['upper'].__setitem__('mu',str(F(x['tilted']['upper']['mu'])*2)))
run('tilted alpha changed',lambda x: x['tilted'].__setitem__('alpha',str(F(x['tilted']['alpha'])+F(1,10))))
run('cone weights not convex',lambda x: x['cone']['lam'].__setitem__(0,str(F(x['cone']['lam'][0])+F(1,10))))
run('honest device angle changed',lambda x: x['tilted']['par_q'].__setitem__(1,str(F(x['tilted']['par_q'][1])+F(1,2))))
