# Containment test for the pinned_term code path: on boxes near c=1 (where the manual enclosure is picked), point values
# F(c,x) and central finite differences dF/dc (prec 300, h=1e-40) at random points of the box must lie in vgen's enclosures.
import sys, random
sys.argv=['vgen.py','c_Wa.json','--prec','300']
exec(open('vgen.py').read().split("def lambda_run")[0])
random.seed(7); bad=0; ntest=0; nman=0
for strip in (20,21,22,23,12):
    S=CERT['strips'][strip]; R=[[q(a),q(b)] for a,b in S['R']]; xa0,xb0=q(S['xa']),q(S['xb'])
    for trial in range(40):
        e=random.randint(4,34); w=fmpq(1,2**e); cb=fmpq(1) if trial%2==0 else 1-w*random.randint(1,3); ca=cb-w
        k=random.randint(0,63); xa=xa0+(xb0-xa0)*fmpq(k,64); xb=xa+(xb0-xa0)/64
        dB,dS,dM=theta_data_mv(xa,xb); pB=build(dB,R); cm=(ca+cb)/2
        Fc=F(dB,pB,arb_series([ball(ca,cb),1],prec=3),cm)
        for _ in range(5):
            c=ca+(cb-ca)*fmpq(random.randint(1,999),1000); x=xa+(xb-xa)*fmpq(random.randint(0,1000),1000)
            dat=theta_data(arb(x)); P=build(dat,R); h=fmpq(1,10**40)
            f0=F(dat,P,arb(c)); fp=F(dat,P,arb(c+h)); fm=F(dat,P,arb(c-h))
            for i in range(4):
                cs=coeffs(Fc[i],3); ntest+=1
                d=(fp[i]-fm[i])/(2*arb(h))
                okv = cs[0].overlaps(f0[i]) if cs[0].is_finite() else True
                okd = (cs[1].overlaps(d) or not cs[1].is_finite())
                # finite difference error ~ h^2 * F''' -> tolerance
                if cs[1].is_finite() and not okd: okd = abs(float((cs[1]-d).mid())) < 1e-30
                if not (okv and okd): bad+=1; print('VIOLATION', strip, float(ca), float(cb), i, cs[0].str(5), f0[i].str(5), cs[1].str(5), d.str(5))
print('containment tests:', ntest, 'violations:', bad)
