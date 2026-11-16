import sys
sys.argv=['tgen3.py',sys.argv[1]]
import tgen3 as T
from flint import arb, fmpq
# boxes inside the float-negative zone near the edge: accept() must NOT accept the negative pair
for (ta,tb) in [(fmpq(0),fmpq(1,4096)),(fmpq(1,2000),fmpq(1,1900)),(fmpq(1,300),fmpq(1,290))]:
    for k in [6,8,10,12,14]:
        ca=fmpq(-1); cb=fmpq(-1)+fmpq(1,2**k)
        d=T.hats(arb(ta+(tb-ta)/2)); g=T.G(d,arb(cb-fmpq(1,2**(k+1))),T.pin(d))
        print(float(ta),k,'G_mid',[x.str(3) for x in g],'accept:',T.accept(ca,cb,ta,tb,{}))
