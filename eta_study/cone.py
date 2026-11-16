import os
"""Normalised normal-cone section C_theta = {g in N_Q(p_theta): g.p_theta = 1} of Thm 2 (hexagon for theta<=theta*,
octagon above), computed as the polygon {g = sum_{a,x} mu_{a|x} gamma_{a|x}, mu.m = 0, g.p = 1, g.L <= 1 for all 16
deterministic L} (slice cone cut by the deterministic points, the description proved in Thm 2 / fullrange/REPORT.md).
Float only (used to generate witnesses; every witness actually used is re-checked: g.p_theta = 1 = h_Q(g) is not
needed for the security proof -- any g is a valid Bell witness -- the cone only labels the family)."""
import numpy as np, itertools, sys
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rigorous_rates'))
from families import ptheta, gammas, DET, hexagon
def polygon(th):
    p,c0,(s,C,r)=ptheta(th); Gm=np.array(gammas(p,c0)); m=np.array([1+C,-s*s,1-C,-s*s])
    # affine plane: mu with mu.m=0 and (mu.Gm).p=1  -> 4-2 = 2 dims
    A=np.vstack([m,Gm@p]); mu0=np.linalg.lstsq(A,[0,1],rcond=None)[0]
    _,_,Vt=np.linalg.svd(A); N=Vt[2:].T          # null space 4x2
    # inequalities (mu0+N u).Gm.L <=1
    Hm=np.array([(N.T@(Gm@L)) for L in DET]); hb=np.array([1-mu0@(Gm@L) for L in DET])
    verts=[]
    for i,j in itertools.combinations(range(16),2):
        M=Hm[[i,j]]
        if abs(np.linalg.det(M))<1e-12: continue
        u=np.linalg.solve(M,hb[[i,j]])
        if np.all(Hm@u<=hb+1e-9):
            if not any(np.linalg.norm(u-w)<1e-8 for w in verts): verts.append(u)
    verts=np.array(verts); c=verts.mean(0); ang=np.arctan2(*(verts-c).T[::-1]); verts=verts[np.argsort(ang)]
    G=np.array([(mu0+N@u)@Gm for u in verts])
    return p,G
if __name__=='__main__':
    for th in (np.pi/8,0.21*np.pi,np.pi/4):
        p,G=polygon(th); print(th/np.pi,len(G),[round(g@p,9) for g in G])
        if th<0.2*np.pi:
            _,order,_=hexagon(th)
            print(' match hexagon:',all(min(np.abs(G-o).max(1))<1e-8 for o in order))
        if th==np.pi/4: print(np.round(G*2*np.sqrt(2),4))
