"""Certify, at every grid efficiency of Theorem 'cover', that CHSH sign variants / deterministic no-click
assignments outside the canonical class are no better than the canonical CHSH with no-click -> +1.
Reduction (analytic): per-input binning a_x,b_y in {+-1} gives E'_xy = a_x b_y * E^{(+1)}_xy(relabelled device), so
every (variant, binning) pair equals some sign pattern s (prod s = -1) on +1-binned data of a relabelled device.
Input swaps of either party preserve Q and the +1 binning and permute the entries of s, so the 8 patterns form two
classes: one -1 (canonical, +++-) and three -1 (class 2, ---+). We certify class 2 with an exact-rational NPA L3
certificate at each grid point and compare with the canonical S_up stored in data/t4."""
import os, sys, json, glob, time
from fractions import Fraction as F
A=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..')
sys.path[:0]=[os.path.join(A,'rigorous_rates'),os.path.join(A,'eta_study'),os.path.join(A,'lossycone')]
import rates_core as R
from verify_all import check_cert
from cert_eta import loss_functional_F, G_CHSH_F
def functional(sv,eta):
    sa=[sv[0]+sv[1],sv[2]+sv[3]]; sb=[sv[0]+sv[2],sv[1]+sv[3]]
    c0=(1-eta)**2*sum(sv)
    c=[eta*(1-eta)*sa[0],eta*(1-eta)*sa[1],eta*(1-eta)*sb[0],eta*(1-eta)*sb[1]]+[eta**2*x for x in sv]
    return c0,c
out=[]
files=sorted(glob.glob(os.path.join(A,'lossycone','data','t4','*.json')))
if len(sys.argv)>1: files=files[:int(sys.argv[1])]
for f in files:
    es=os.path.basename(f)[:-5]; e=F(es); st=json.load(open(f))
    c0,ge=loss_functional_F(G_CHSH_F,e); b1,_=check_cert(st['chsh_cert'],ge); S1=c0+b1
    t=time.time()
    c02,c2=functional([-1,-1,-1,1],e)
    beta,L,cj,qv=R.certify_hq(c2)
    b2,_=check_cert(cj,c2)  # independent exact re-check
    S2=c02+b2
    ok=S2<=S1
    out.append(dict(eta=es,S1_canonical=float(S1),S2_class2=float(S2),qubit_class2=float(c02)+qv,margin=float(S1-S2),ok=ok,cert=cj))
    print('%s canonical S_up=%.8f class2 S_up=%.8f (qubit %.8f) margin %.3e %s  [%.1fs]'%(es,float(S1),float(S2),float(c02)+qv,float(S1-S2),'OK' if ok else 'FAIL',time.time()-t),flush=True)
json.dump(out,open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'class2_certs.json'),'w'))
print('ALL OK =',all(o['ok'] for o in out),len(out))
