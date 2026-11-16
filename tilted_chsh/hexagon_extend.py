"""Goal 4 (hexagon part): rigorous interval proof of the hexagon structure on intervals beyond Theta."""
import time, arbeval, hexagon_m2 as H
from flint import arb, fmpq
H.theta_vars = lambda x: arbeval.theta_vars(x * arb.pi())   # boxes in x = theta/pi
for lo, hi in [(fmpq(1, 16), fmpq(1, 8)), (fmpq(1, 32), fmpq(1, 16)), (fmpq(1, 128), fmpq(1, 32)), (fmpq(1, 1024), fmpq(1, 128)),
               (fmpq(3, 16), fmpq(2022, 10000))]:
    t0 = time.time(); res = H.interval_proof(lo, hi, 200)
    if res[0]:
        w = sorted(res[2].items(), key=lambda kv: kv[1])[:2]
        print("[%s, %s] pi: hexagon PROVED (%d boxes, %.1fs); smallest lower bounds %s" % (lo, hi, res[1], time.time() - t0, [(k, '%.3g' % v) for k, v in w]), flush=True)
    else:
        print("[%s, %s] pi: FAILED at %s on %s" % (lo, hi, res[1], res[2]), flush=True)
