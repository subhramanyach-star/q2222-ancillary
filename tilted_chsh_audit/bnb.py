"""Independent branch-and-bound for Theorem 2, checks (d) coverage and (f) strictness.
Proves, for vertex V in {V1, V3} and every (c, s) in [-1,1] x S_Theta  (s = sin 2theta),
    F_{s0 s1}(c, s) = 1 - s0 a0 - s1 a1 - ||B_0(delta)+s0 e_0||_{G_c} - ||B_1(delta)+s1 e_1||_{G_c} >= 0,
with F > 0 unless c = c0(s), or c = +-1 and (s0,s1) is a tight sign pair at that endpoint.
Inputs: exact formulas parsed from prl.tex (model.py).  Proven exactly elsewhere (check_exact.py):
  F(c0(s),s) = dF/dc(c0(s),s) = 0 for all four pairs, and F(+-1, s) >= 0 (endpoint lemma).
Tests per (c-interval I, s-strip S, pair):
  T0  natural/centred enclosure of F on I x S has lower bound > 0
  T1  mean-value form F(cm) + F_c(I)(I-cm) > 0                       (needs A,B > 0 on I)
  T2  F_cc > 0 on hull(I u c0(S)) x S   (+ exact F=F_c=0 at c0)         (needs A,B > 0 on hull)
  T3  I touches c=+1 (resp. -1) and F_c < 0 (resp. > 0) on I x S  (+ exact F(+-1,s) >= 0)
Coverage: s-strips are exact rational intervals tiling S_Theta; c-leaves tile [-1,1]; both are
checked with exact rational arithmetic."""
import sys, json, time, os, argparse, random
from fractions import Fraction as Fr
from flint import arb, arb_poly, fmpq, ctx
import model as M
from kfun import KF, kderiv, ball
from kfield import K
ctx.prec = 128

ap = argparse.ArgumentParser()
ap.add_argument('vertex', choices=['V1', 'V3'])
ap.add_argument('--nstrips', type=int, default=64)
ap.add_argument('--only', type=str, default='')          # comma list of strip indices
ap.add_argument('--sample', type=int, default=0)         # random sample of strips
ap.add_argument('--seed', type=int, default=1)
ap.add_argument('--out', type=str, default='')
ap.add_argument('--minw', type=float, default=2.0**-34)
args = ap.parse_args()
name = args.vertex
Vd = M.vertex_from_tex('V_1' if name == 'V1' else 'V_3')
g = Vd['g']
# work in the G_c-orthogonal basis f_+ = (1,1), f_- = (1,-1):  ||u||_{G_c}^2 = (1+c)/2 (u.f_+)^2 + (1-c)/2 (u.f_-)^2
QTY = {'a0': g[0], 'a1': g[1], 'bp': g[2] + g[3], 'bm': g[2] - g[3],
       'e0p': g[4] + g[5], 'e0m': g[4] - g[5], 'e1p': g[6] + g[7], 'e1m': g[6] - g[7],
       'd0p': Vd['d0'][0] + Vd['d0'][1], 'd0m': Vd['d0'][0] - Vd['d0'][1],
       'd1p': Vd['d1'][0] + Vd['d1'][1], 'd1m': Vd['d1'][0] - Vd['d1'][1], 'c0': M.c0}
FUN = {k: (KF(v), KF(kderiv(v))) for k, v in QTY.items()}
R = M.TABLE[0] if name == 'V1' else M.TABLE[1]
Rq = {jl: (fmpq(v[0].p, v[0].q), fmpq(v[1].p, v[1].q)) for jl, v in R.items()}
kdeg = max(j for j, l in R)

# tight pairs at the endpoints (from the facet points; confirmed exactly in check_exact.py)
TIGHT = {1: set(), -1: set()}
for L in M.FACETS[name]:
    TIGHT[L[2]*L[3]].add((L[0], L[1]))
PAIRS = [(1, 1), (1, -1), (-1, 1), (-1, -1)]

import mpmath as mp
mp.mp.dps = 50
S_LO = Fr(int(mp.floor(mp.sin(mp.pi/4) * 10**30)), 10**30)
S_HI = Fr(int(mp.ceil(mp.sin(3*mp.pi/8) * 10**30)), 10**30)

def q2f(x): return fmpq(x.numerator, x.denominator)

def enclose(key, S, sm):
    f, fp = FUN[key]
    cen = f(sm) + fp(S) * (S - sm)
    nat = f(S)
    lo = max(cen.lower(), nat.lower()); hi = min(cen.upper(), nat.upper())
    return ball(lo, hi) if lo <= hi else cen


HALF = arb(fmpq(1, 2))
class Quad:
    """A(c) = (1+c)/2 p(c)^2 + (1-c)/2 m(c)^2 = ||u(c)||_{G_c}^2, u.f_+ = p, u.f_- = m."""
    def __init__(self, p, m):
        self.p = (p, p.derivative(), p.derivative().derivative())
        self.m = (m, m.derivative(), m.derivative().derivative())
    def _pm(self, x):
        return [f(x) for f in self.p], [f(x) for f in self.m]
    def val(self, x):
        (p, _, _), (m, _, _) = self._pm(x)
        return HALF*(1 + x)*p*p + HALF*(1 - x)*m*m
    def d1(self, x):
        (p, p1, _), (m, m1, _) = self._pm(x)
        return HALF*(p*p - m*m) + (1 + x)*p*p1 + (1 - x)*m*m1
    def d2(self, x):
        (p, p1, p2), (m, m1, m2) = self._pm(x)
        return 2*(p*p1 - m*m1) + (1 + x)*(p1*p1 + p*p2) + (1 - x)*(m1*m1 + m*m2)
    def __getitem__(self, i):
        return (self.val, self.d1, self.d2)[i]

def strip_polys(slo, shi):
    S = ball(q2f(slo), q2f(shi)); sm = arb(q2f((slo + shi) / 2))
    E = {k: enclose(k, S, sm) for k in QTY}
    sig = (S - arb(fmpq(83, 100))) * 10
    c = arb_poly([0, 1])
    h = c - E['c0']
    Pc = [arb_poly([0]), arb_poly([0])]
    for (j, l), (R0, R1) in Rq.items():
        mon = arb_poly([0]*j + [1]) * (sig ** l)
        Pc[0] = Pc[0] + mon * arb(R0); Pc[1] = Pc[1] + mon * arb(R1)
    dp = E['d0p'] + E['d1p']*h + h*h*(Pc[0] + Pc[1])        # delta.f_+
    dm = E['d0m'] + E['d1m']*h + h*h*(Pc[0] - Pc[1])        # delta.f_-
    comp = {}
    for a in (1, -1):
        up = E['bp']/2 + dp + a*E['e0p']; um = E['bm']/2 + dm + a*E['e0m']
        vp = E['bp']/2 - dp + a*E['e1p']; vm = E['bm']/2 - dm + a*E['e1m']
        comp[('A', a)] = (up, um)
        comp[('B', a)] = (vp, vm)
    T = {(s0, s1): 1 - s0*E['a0'] - s1*E['a1'] for s0, s1 in PAIRS}
    polys = {}
    for key, (p, m) in comp.items():
        polys[key] = Quad(p, m)
    return polys, T, E['c0']

def sqrtp(x):
    """sqrt of an enclosure of a quantity known to be >= 0"""
    if x.lower() < 0:
        return ball(arb(0), x.upper().sqrt().upper()) if x.upper() >= 0 else arb(0)
    return x.sqrt()

def pos_lower(x): return x.lower() > 0

def evalq(polys, key, I, cm):
    Q = polys[key]
    nat = Q.val(I); cen = Q.val(cm) + Q.d1(I) * (I - cm)
    lo = max(nat.lower(), cen.lower()); hi = min(nat.upper(), cen.upper())
    return ball(lo, hi) if lo <= hi else nat

def test_box(polys, T, c0enc, cl, ch, pair):
    I = ball(q2f(cl), q2f(ch)); cm = arb(q2f((cl + ch) / 2))
    s0, s1 = pair
    A = evalq(polys, ('A', s0), I, cm); B = evalq(polys, ('B', s1), I, cm)
    Tv = T[pair]
    # T0
    F = Tv - sqrtp(A) - sqrtp(B)
    if F.lower() > 0: return 'T0'
    if not (A.lower() > 0 and B.lower() > 0):
        return None
    dA = polys[('A', s0)][1](I); dB = polys[('B', s1)][1](I)
    sA, sB = A.sqrt(), B.sqrt()
    Fp = -dA / (2*sA) - dB / (2*sB)
    # T3 endpoint monotonicity
    if ch == 1 and Fp.upper() < 0: return 'T3+'
    if cl == -1 and Fp.lower() > 0: return 'T3-'
    # T3b second-order endpoint test: F(c) >= F(e) + F'(e)(c-e) + (1/2) min(0, F''(J)) (c-e)^2,
    #     J = hull(e, I), e = +-1, using F(e, s) >= 0 for all s (exact, endpoint lemma)
    for e in (-1, 1):
        dist = (Fr(ch) + 1) if e == -1 else (1 - Fr(cl))
        if dist > Fr(1, 8): continue
        J = ball(q2f(Fr(-1)), q2f(ch)) if e == -1 else ball(q2f(cl), q2f(Fr(1)))
        AJ = polys[('A', s0)][0](J); BJ = polys[('B', s1)][0](J)
        if not (AJ.lower() > 0 and BJ.lower() > 0): continue
        ce = arb(e)
        Ae = polys[('A', s0)][0](ce); Be = polys[('B', s1)][0](ce)
        Fpe = -polys[('A', s0)][1](ce) / (2*Ae.sqrt()) - polys[('B', s1)][1](ce) / (2*Be.sqrt())
        slope = Fpe if e == -1 else -Fpe          # outward-to-inward slope, must be > 0
        dAJ = polys[('A', s0)][1](J); dBJ = polys[('B', s1)][1](J)
        d2AJ = polys[('A', s0)][2](J); d2BJ = polys[('B', s1)][2](J)
        sAJ, sBJ = AJ.sqrt(), BJ.sqrt()
        FppJ = -(d2AJ/(2*sAJ) - dAJ*dAJ/(4*AJ*sAJ)) - (d2BJ/(2*sBJ) - dBJ*dBJ/(4*BJ*sBJ))
        curv = FppJ.lower()
        if curv > 0: curv = arb(0)
        # for t = |c - e| in (0, dist]: F >= t (slope + curv t / 2) >= t * min(slope + curv*dist/2, slope)
        lb = slope.lower() + curv * arb(q2f(dist)) / 2
        if lb > 0: return 'T3b%+d' % e
    # T1
    Am = polys[('A', s0)][0](cm); Bm = polys[('B', s1)][0](cm)
    if Am.lower() > 0 and Bm.lower() > 0:
        Fm = Tv - Am.sqrt() - Bm.sqrt()
        F1 = Fm + Fp * (I - cm)
        if F1.lower() > 0: return 'T1'
    # T2 near c0
    lo = min(Fr(cl), Fr(0)) if False else None
    c0l, c0u = c0enc.lower(), c0enc.upper()
    if arb(q2f(cl)) - c0u < 0.3 and c0l - arb(q2f(ch)) < 0.3:
        Hlo = arb(q2f(cl)) if arb(q2f(cl)) < c0l else c0l
        Hhi = arb(q2f(ch)) if arb(q2f(ch)) > c0u else c0u
        H = ball(Hlo.lower(), Hhi.upper())
        AH = polys[('A', s0)][0](H); BH = polys[('B', s1)][0](H)
        if AH.lower() > 0 and BH.lower() > 0:
            dAH = polys[('A', s0)][1](H); dBH = polys[('B', s1)][1](H)
            d2AH = polys[('A', s0)][2](H); d2BH = polys[('B', s1)][2](H)
            sAH, sBH = AH.sqrt(), BH.sqrt()
            Fpp = -(d2AH/(2*sAH) - dAH*dAH/(4*AH*sAH)) - (d2BH/(2*sBH) - dBH*dBH/(4*BH*sBH))
            if Fpp.lower() > 0: return 'T2'
    return None

def run_strip(slo, shi, minw):
    polys, T, c0enc = strip_polys(slo, shi)
    stats = {'T0': 0, 'T1': 0, 'T2': 0, 'T3+': 0, 'T3-': 0, 'T3b-1': 0, 'T3b+1': 0}
    leaves = 0
    stack = [(Fr(-1), Fr(1), tuple(PAIRS))]
    covered = []          # c-leaves in order of acceptance (checked for tiling)
    while stack:
        cl, ch, pend = stack.pop()
        rest = []
        for pr in pend:
            t = test_box(polys, T, c0enc, cl, ch, pr)
            if t is None: rest.append(pr)
            else: stats[t] += 1
        if not rest:
            leaves += 1; covered.append((cl, ch)); continue
        if ch - cl < minw:
            return False, {'fail_c': [float(cl), float(ch)], 'pairs': rest, 'stats': stats}
        cm = (cl + ch) / 2
        stack.append((cm, ch, tuple(rest))); stack.append((cl, cm, tuple(rest)))
    # coverage of [-1,1] by leaves (exact): leaves are produced in increasing order by the DFS
    covered.sort()
    ok = covered[0][0] == -1 and covered[-1][1] == 1 and all(covered[i][1] == covered[i+1][0] for i in range(len(covered)-1))
    return ok, {'leaves': leaves, 'stats': stats, 'cover_ok': ok}

def strip_bounds(n):
    w = (S_HI - S_LO) / n
    return [(S_LO + w*i, S_LO + w*(i+1)) for i in range(n)]

def prove_strip(slo, shi, depth=0, log=None):
    """recursively split the s-strip until every c-box is accepted; returns list of accepted sub-strips"""
    t0 = time.time()
    ok, info = run_strip(slo, shi, args.minw)
    if ok:
        return [(slo, shi, info['leaves'], info['stats'], time.time() - t0)], 0
    if depth > 18:
        raise RuntimeError('FAILED strip %s %s: %s' % (float(slo), float(shi), info))
    sm = (slo + shi) / 2
    a, na = prove_strip(slo, sm, depth+1); b, nb = prove_strip(sm, shi, depth+1)
    return a + b, 1 + na + nb

if __name__ == '__main__':
    strips = strip_bounds(args.nstrips)
    idx = list(range(args.nstrips))
    if args.only: idx = [int(x) for x in args.only.split(',')]
    if args.sample: random.seed(args.seed); idx = sorted(random.sample(idx, args.sample))
    out = args.out or ('bnb_%s.jsonl' % name)
    done = set()
    if os.path.exists(out):
        for line in open(out):
            try: done.add(json.loads(line)['strip'])
            except Exception: pass
    print('vertex', name, 'tight pairs', {k: sorted(v) for k, v in TIGHT.items()}, 'strips', len(idx), 'done', len(done), flush=True)
    T0 = time.time()
    for i in idx:
        if i in done: continue
        slo, shi = strips[i]
        t0 = time.time()
        try:
            subs, nsplit = prove_strip(slo, shi)
            ok = True; err = ''
        except RuntimeError as ex:
            subs, nsplit, ok, err = [], -1, False, str(ex)
        # exact tiling of the strip by accepted sub-strips
        tile_ok = ok and subs[0][0] == slo and subs[-1][1] == shi and all(subs[k][1] == subs[k+1][0] for k in range(len(subs)-1))
        agg = {}
        for sb in subs:
            for kk, vv in sb[3].items(): agg[kk] = agg.get(kk, 0) + vv
        rec = dict(strip=i, s_lo=str(slo), s_hi=str(shi), theta_over_pi=[float(mp.asin(mp.mpf(slo.numerator)/slo.denominator)/(2*mp.pi)), float(mp.asin(mp.mpf(shi.numerator)/shi.denominator)/(2*mp.pi))],
                   ok=bool(ok and tile_ok), substrips=len(subs), splits=nsplit, cboxes=sum(sb[2] for sb in subs), tests=agg,
                   seconds=round(time.time() - t0, 2), err=err,
                   finished=time.strftime('%Y-%m-%d %H:%M:%S %Z'))
        with open(out, 'a') as fh: fh.write(json.dumps(rec) + '\n')
        print(json.dumps(rec), flush=True)
    print('total %.1f s' % (time.time() - T0))
