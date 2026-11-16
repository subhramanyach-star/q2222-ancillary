import os
"""Does any other CHSH variant (sign pattern s_xy with prod s = -1) or any other deterministic binning of no-clicks
do better than the canonical CHSH (+++-) with no-click -> +1?  Binning a party's no-click to -1 equals relabelling
that party's binned outcome, which maps the canonical functional to another sign variant evaluated on +1-binned data,
so it suffices to compare the 8 variants (x 4 marginal sign choices absorbed).  For each variant s, the maximum of
S_s over all devices at efficiency eta is c0_s + h_Q(g_s,eta) (loss map exact).  We print the exact-rational NPA upper
bound for each variant and the canonical qubit value: every variant's certified bound must be <= canonical S_up."""
import sys, itertools, json, warnings; warnings.simplefilter('ignore')
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','rigorous_rates'))
from fractions import Fraction as F
import rates_core as R
from cert_eta import loss_functional_F
rows=[]
for eta in [F(7,10),F(3,4),F(4,5),F(17,20),F(9,10),F(19,20)]:
    out=[]
    for s in itertools.product((1,-1),repeat=3):
        sv=list(s)+[-s[0]*s[1]*s[2]]
        for ma in (1,-1):          # binned-outcome relabelling of Alice / Bob (no-click -> -1 on that side)
            for mb in (1,-1):
                m=ma*mb; sa=[sv[0]+sv[1],sv[2]+sv[3]]; sb=[sv[0]+sv[2],sv[1]+sv[3]]
                # binned A' = ma*(eta A + (1-eta)), B' = mb*(eta B + (1-eta)):  sum s E' = m[eta^2 sum s AB + eta(1-eta) sum s (A_x+B_y) + (1-eta)^2 sum s]
                c0=m*(1-eta)**2*sum(sv)
                c=[m*eta*(1-eta)*sa[0],m*eta*(1-eta)*sa[1],m*eta*(1-eta)*sb[0],m*eta*(1-eta)*sb[1]]+[m*eta**2*x for x in sv]
                beta,L,cj,qv=R.certify_hq(c)
                out.append((float(c0+beta),float(c0)+qv,tuple(sv),ma,mb))
    out.sort(key=lambda t:-t[0])
    can=[o for o in out if o[2]==(1,1,1,-1) and o[3]==1 and o[4]==1][0]
    print('eta=%s canonical S_up=%.8f; best over all variants/binnings: cert %.8f qubit %.8f (s=%s, binA=%+d, binB=%+d)'%(eta,can[0],out[0][0],out[0][1],out[0][2],out[0][3],out[0][4]),flush=True)
    rows.append(dict(eta=str(eta),canonical=can[0],best=out[0][0],best_qubit=out[0][1],ok=out[0][1]<=can[1]+1e-9 and out[0][0]<=can[0]+1e-9))
json.dump(rows,open('out/variants_check.json','w'),indent=1)
