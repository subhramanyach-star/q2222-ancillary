"""Exact combinatorial structure of the real NPA moment matrix for the CHSH (2222) scenario.
Binary observables A0,A1,B0,B1 (A_x^2=B_y^2=1, [A_x,B_y]=0).  Words are reduced (no repeated letter),
Alice letters before Bob letters.  Entry (u,v) of the moment matrix is <u^dagger v>; for a REAL relaxation
<w> and <w^dagger> are identified (Re of a complex moment matrix is feasible).  The canonical key identifies
(a,b) with (reverse(a), reverse(b)), reversing both party words simultaneously.  Self-contained (no solver),
used by both the certifier and the verifier.
Coordinates z=(a0,a1,b0,b1,E00,E01,E10,E11)."""
import itertools
OPS_A=('A0','A1'); OPS_B=('B0','B1')
def _red(s):
    out=[]
    for o in s:
        if out and out[-1]==o: out.pop()
        else: out.append(o)
    return tuple(out)
def reduce(word):
    return _red([o for o in word if o[0]=='A']), _red([o for o in word if o[0]=='B'])
def canon(word):
    a,b=reduce(word)
    return min((a,b),(tuple(reversed(a)),tuple(reversed(b))))
def words(level):
    if level=='1+AB':
        return [()]+[(o,) for o in OPS_A+OPS_B]+[(x,y) for x in OPS_A for y in OPS_B]
    L=int(level); W=set()
    for n in range(L+1):
        for w in itertools.product(OPS_A+OPS_B,repeat=n):
            a,b=reduce(w); W.add(a+b)
    return sorted(W,key=lambda w:(len(w),w))
class Struct:
    def __init__(s,level):
        s.level=level; s.W=words(level); n=len(s.W); s.n=n
        keys={}; s.idx=[[0]*n for _ in range(n)]
        for i,u in enumerate(s.W):
            for j,v in enumerate(s.W):
                k=canon(tuple(reversed(u))+v)
                if k not in keys: keys[k]=len(keys)
                s.idx[i][j]=keys[k]
        s.keys=keys; s.m=len(keys)
        s.k1=keys[((),())]
        s.kz=[keys[((x,),())] for x in OPS_A]+[keys[((),(y,))] for y in OPS_B]+[keys[((x,),(y,))] for x in OPS_A for y in OPS_B]
        s.cls=[[] for _ in range(s.m)]
        for i in range(n):
            for j in range(n): s.cls[s.idx[i][j]].append((i,j))
