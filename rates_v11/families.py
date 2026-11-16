import numpy as np, itertools
s2=np.sqrt(2)
DET=[np.array([A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1],float) for A0,A1,B0,B1 in itertools.product([1,-1],repeat=4)]
def detpt(A0,A1,B0,B1): return np.array([A0,A1,B0,B1,A0*B0,A0*B1,A1*B0,A1*B1],float)
def hL(g): return max(g@L for L in DET)
def ptheta(th):
    s=np.sin(2*th);C=np.cos(2*th);r=np.sqrt(1+s*s)
    p=np.array([C,0,C/r,C/r,1/r,1/r,s*s/r,-s*s/r]); c0=(1-s*s)/(1+s*s)
    return p,c0,(s,C,r)
def gammas(p,c0):
    G=np.array([[1,c0],[c0,1]]);Gi=np.linalg.inv(G)
    al=p[:2];be=p[2:4];E=p[4:].reshape(2,2)
    out={}
    for x in (0,1):
        for a in (1,-1):
            t=1+a*al[x]; w=be+a*E[x]; v=Gi@w/t
            g=np.zeros(8); g[x]=-a; g[2:4]=v; g[4+2*x:6+2*x]=a*v
            out[(a,x)]=g
    order=[(1,0),(1,1),(-1,0),(-1,1)]
    return [out[k] for k in order]
def hexagon(th):
    p,c0,(s,C,r)=ptheta(th); Gm=gammas(p,c0)
    m=np.array([1+C,-s*s,1-C,-s*s])
    facets={'V1':[(1,1,1,1),(1,-1,-1,1)],'V2':[(1,-1,-1,1),(1,1,1,-1)],'V3':[(1,-1,1,1),(-1,1,1,-1)],'V4':[(-1,1,1,-1),(-1,-1,-1,1)]}
    V={}
    for k,(L1,L2) in facets.items():
        M=np.array([m,[g@p for g in Gm],[g@detpt(*L1) for g in Gm],[g@detpt(*L2) for g in Gm]])
        mu=np.linalg.solve(M,[0,1,1,1]); V[k]=(sum(mu[i]*Gm[i] for i in range(4)),mu)
    sig=lambda g: np.array([g[0],-g[1],g[3],g[2],g[5],g[4],-g[7],-g[6]])  # A1->-A1, B0<->B1
    # check sigma: z=(a0,a1,b0,b1,E00,E01,E10,E11) -> (a0,-a1,b1,b0,E01,E00,-E11,-E10)
    order=[V['V1'][0],V['V2'][0],sig(V['V1'][0]),V['V3'][0],V['V4'][0],sig(V['V3'][0])]
    return p,order,V
def gtilt(th):
    p,c0,(s,C,r)=ptheta(th); al=2*C/r
    g=np.array([al,0,0,0,1,1,1,-1])/np.sqrt(8+2*al*al); return g
g_chsh=np.array([0,0,0,0,1,1,1,-1],float)
def Idelta(d):
    E=np.array([0,np.cos(d),np.cos(d),-np.sin(2*d)]); p=np.concatenate([np.zeros(4),E])
    IQ=2*np.cos(d)**3/(np.cos(2*d)*np.sin(d))
    g=np.array([0,0,0,0,1,1/np.sin(d),1/np.sin(d),-1/np.cos(2*d)])/IQ
    c0=-np.sin(d);G=np.array([[1,c0],[c0,1]]);Gi=np.linalg.inv(G)
    Y=[]
    for x in (0,1):
        y=np.zeros(8); y[2:4]=Gi@E[2*x:2*x+2]; y[x]=-1; Y.append(y)
    cd=np.cos(d);sd=np.sin(d);c2=np.cos(2*d)
    W1=np.array([-(1-cd)*c2/(2*cd**3), sd/(2*cd**3)]); W2=np.array([-sd*c2/(2*cd**3),(cd-c2)/(2*cd**3)])
    return p,g,Y,W1,W2
if __name__=='__main__':
    for th in [np.pi/8,3*np.pi/16]:
        p,order,V=hexagon(th)
        print('theta/pi',th/np.pi)
        for g in order: print('  g.p',round(g@p,12),'hL',round(hL(g),12),'tight',sum(abs(g@L-1)<1e-9 for L in DET))
        for k in V: print('  mu',k,V[k][1])
        s=np.sin(2*th);C=np.cos(2*th);r=np.sqrt(1+s*s);D=r+1
        k_=(C+r+2*s*s)/(4*D);n_=(r+2-C)/(4*D)
        print('  V2 closed',np.allclose(V['V2'][0],[(r+C)/(2*D),0,(C-1)/(2*D),(C-1)/(2*D),k_,k_,n_,-n_]))
        gt=gtilt(th);print('  gtilt.p',gt@p,'hL',hL(gt))
    p,g,Y,W1,W2=Idelta(np.pi/8)
    for d in [W1,W2,-W1,-W2,np.zeros(2)]:
        gg=g+d[0]*Y[0]+d[1]*Y[1]; print('Idelta vertex g.p',gg@p,'hL',hL(gg))
