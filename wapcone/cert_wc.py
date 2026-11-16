"""Rigorous H_min(A1|X=1,E) for the best single witness of each family (fam.jsonl best record), v=1:
witness g = dual combination of the family rows at the optimised device (g = lam.R_ge - nu.R_eq), normalised, rationalised
(denominators <= 1e12), then cert_eta.winner: honest score in arb >= omega (dyadic), exact-rational NPA duals (levels 1+AB/2/3)
+ arb => P_guess <= P_up for EVERY quantum device with g.z >= omega; explicit qubit attack => two-sided interval.
Any rational g is a valid Bell witness; family membership only labels where g came from.  Writes certs/<eta>_<fam>.json, cert.jsonl."""
import sys, json, os, warnings; warnings.simplefilter('ignore')
_H=os.path.dirname(os.path.abspath(__file__)); os.chdir(_H)
sys.path.insert(0,os.path.join(_H,'..','eta_study')); sys.path.insert(0,os.path.join(_H,'..','rigorous_rates'))
import numpy as np
from fractions import Fraction as F
import cert_eta as CE
from flint import arb
import fam as FM
def witness(eta,f,rec):
    x=np.array(rec['par']); z=FM.z_honest(x[:5],eta)
    if f=='wap': P=FM.pg('wap',0,2); P.H(z,Re=FM.R_wap); g=-P.duals()[1]@FM.R_wap
    elif f=='full': P=FM.pg('f',0,8); P.H(z,Re=FM.I8); g=-P.duals()[1]
    elif f=='tilted': P=FM.pg('t',1,1); P.H(z,Rg=FM.G_CHSH[None,:],Re=FM.I8[[0]]); dg,de=P.duals(); g=dg[0]*FM.G_CHSH-de[0]*FM.I8[0]
    else:
        n=FM.NROW[f]; R=FM.pad(FM.rows(f,x,eta),n); P=FM.pg(f,n,0); P.H(z,Rg=R); g=np.maximum(P.duals()[0],0)@R
    return g/np.abs(g).max()
def run(eta,f):
    rec=FM.best(eta,f); g=witness(eta,f,rec)
    gq=[F(float(t)).limit_denominator(10**12) for t in g]; ga=[CE.aq(t) for t in gq]
    w=CE.winner(rec['par'][:5],F(eta).limit_denominator(1000),F(1),ga,f)
    w['g']=[str(t) for t in gq]; w['g_float']=list(map(float,g)); w['H_float_level2']=rec['H']; w['fam_par']=rec['par']
    json.dump(w,open('certs/eta%s_%s.json'%(eta,f),'w'))
    s=dict(eta=eta,fam=f,H_lo=w['H_lo'],H_up=w['H_up'],H_float=rec['H'],g=list(map(float,g)))
    with open('cert.jsonl','a') as fo: fo.write(json.dumps(s)+'\n')
    print(json.dumps(s),flush=True)
if __name__=='__main__':
    for e in sys.argv[1].split(','):
        for f in sys.argv[2].split(','): run(float(e),f)
