import numpy as np, time
from bffcert import Rate
from families import ptheta, detpt
R1=Rate(m=8,xstar=1); R0=Rate(m=8,xstar=0)
s2=np.sqrt(2)
zT=np.array([0,0,0,0,1/s2,1/s2,1/s2,-1/s2])
for lab,R,z in [('Tsirelson v=0.999 x*=0',R0,0.999*zT),('Tsirelson v=0.999 x*=1',R1,0.999*zT),('deterministic ++++ (expect<=0)',R1,detpt(1,1,1,1)),('mixed det/2 (z=0) x*=1 (expect <=1, true=0? no: z=0 is uniform noise -> true inf 0)',R1,np.zeros(8))]:
    t0=time.time(); r=R.bound('full',z=z); print(lab,'cert %.6f num %.6f'%(r['H_cert'],r['H_num']),'%.0fs'%(time.time()-t0),flush=True)
for th in [np.pi/8, 27.76/180*np.pi]:
    p,c0,_=ptheta(th)
    for v in [0.9999,0.99]:
        t0=time.time(); r=R1.bound('full',z=v*p); print('theta/pi %.4f v=%.4f full-stat H(A1|E) cert %.6f num %.6f'%(th/np.pi,v,r['H_cert'],r['H_num']),'%.0fs'%(time.time()-t0),flush=True)
