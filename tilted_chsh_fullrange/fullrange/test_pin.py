import sys
sys.argv=['vgen.py','c_Wa.json']
exec(open('vgen.py').read().split("def lambda_run")[0])
import random
FORCE_MANUAL = True
random.seed(1)
for trial in range(6):
    eps=random.choice([1,-1]); c0=fmpq(random.randint(-900,900),1000)
    a=[fmpq(random.randint(-500,500),1000) for _ in range(4)]
    c=arb_series([arb(c0),1],prec=3)
    rho=[a[0]+a[1]*c, a[2]+a[3]*c]
    direct=(1-eps*c)*gnorm(rho,c)
    # force manual: compute with helper internals by making a copy with cs non-finite impossible; compare via finite difference
    man=pinned_term(rho,c,eps)
    h=fmpq(1,10**12)
    f=lambda cc: float(((1-eps*cc)*gnorm([a[0]+a[1]*cc,a[2]+a[3]*cc],cc)).mid())
    fd=(f(arb(c0+h))-f(arb(c0-h)))/(2*float(h))
    print(eps, float(c0), coeffs(direct,3)[1].str(8), coeffs(man,3)[1].str(8), fd)
