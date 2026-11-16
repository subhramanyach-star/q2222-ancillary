#!/usr/bin/env python3
"""Stand-alone re-check (no SDP solver) of the certified LOWER bounds on H_min(A1|X=1,E) stored in
  ../wapcone/certs/*.json   (asymmetric-CHSH and polygon columns of Table S2 of prl_v18)
  ../noise/certs/*.json     (tilted and polygon lower bounds, and the CHSH cap, of Table S4 of prl_v18).
For each stored witness block {par_q, g, omega, upper}:
 (1) honest side: the stored rational qubit device (angles par_q) under the loss (+ noise) model reaches g.z >= omega,
     evaluated in 200-bit ball arithmetic;
 (2) adversary side: P_guess <= max_e[c0_e + beta_e + mu*|g - gtilde|_1] - mu*omega, with every beta_e from an exact
     rational NPA dual certificate re-checked by exact LDL^T (verify_all.check_cert); H_lo = -log2 P_up must reproduce the
     stored value (to 1e-9, rounding direction checked).
For noise/ certificates also (3) the CHSH cap: S_up = c0 + v*(p*beta + (1-p)*L) (white noise v, dephasing p), beta
re-checked exactly, L = 4 eta(1-eta) + 2 eta^2, compared with the stored S_up exactly, and H_cap = f(S_up) upper end.
Prints one line per block and 'ALL OK = True' at the end."""
import os, sys, glob, json
from fractions import Fraction as Fr
_H = os.path.dirname(os.path.abspath(__file__)); sys.path[:0] = [os.path.join(_H, '..', 'rigorous_rates'), os.path.join(_H, '..', 'eta_study')]
from flint import arb, ctx
ctx.prec = 200
from verify_all import check_cert, guess_funcs
from cert_eta import loss_functional_F, G_CHSH_F
def aq(q): q = Fr(q); return arb(q.numerator)/q.denominator
def up(b): return float(b.mid()) + float(b.rad())
def lo(b): return float(b.mid()) - float(b.rad())
LOG2 = arb(2).log()
def honest(par, eta, v=1, pc=1):
    th = aq(par[0]); a = [aq(x) for x in par[1:5]]; C = (2*th).cos(); S = (2*th).sin(); eta = aq(eta); v = aq(v); pc = aq(pc)
    A = [v*a[0].cos()*C, v*a[1].cos()*C]; B = [v*a[2].cos()*C, v*a[3].cos()*C]
    E = [v*(a[x].cos()*a[2+y].cos() + pc*a[x].sin()*a[2+y].sin()*S) for x in (0, 1) for y in (0, 1)]
    A2 = [eta*t + (1-eta) for t in A]; B2 = [eta*t + (1-eta) for t in B]
    E2 = [eta*eta*E[2*x+y] + eta*(1-eta)*(A[x]+B[y]) + (1-eta)**2 for x in (0, 1) for y in (0, 1)]
    return A2 + B2 + E2
def witness(blk):
    if 'g' in blk: return [aq(Fr(x)) for x in blk['g']]
    if 'alpha' in blk: return [aq(Fr(blk['alpha']))] + [arb(0)]*3 + [arb(1), arb(1), arb(1), arb(-1)]
    # polygon witness: rational convex weights over the vertices of C_theta' (vertices solved in ball arithmetic)
    from cert_eta import polygon_arb
    p, verts = polygon_arb(Fr(blk['thetap'])); lam = [Fr(x) for x in blk['lam']]
    assert all(l >= 0 for l in lam) and sum(lam) == 1
    return [sum((aq(l)*V[j] for l, V in zip(lam, verts)), arb(0)) for j in range(8)]
def lower(blk, eta, v=1, pc=1):
    g = witness(blk); om = Fr(blk['omega']); par = [Fr(x) for x in blk['par_q']]
    z = honest(par, eta, v, pc); sc = sum((gi*zi for gi, zi in zip(g, z)), arb(0)) - aq(om)
    ok_h = lo(sc) >= 0
    u = blk['upper']; mu = Fr(u['mu']); gt = [Fr(x) for x in u['gtilde']]
    err = sum((abs(x - aq(y)) for x, y in zip(g, gt)), arb(0)); terms = []
    for (c0, vv, lab), ce in zip(guess_funcs(('A', 1)), u['per_e']):
        c = [vv[j] + mu*gt[j] for j in range(8)]
        beta, L = check_cert(ce['cert'], c)
        terms.append(aq(c0) + aq(beta) + aq(mu)*err - aq(mu)*aq(om))
    Pup = max(up(t) for t in terms)
    Hlo = lo(-(aq(Fr(Pup))).log()/LOG2) if Pup < 1 else 0.0
    return ok_h, Hlo
def P(x): return json.loads(x) if isinstance(x, str) else x
allok = True; n = 0
for f in sorted(glob.glob(os.path.join(_H, '..', 'wapcone', 'certs', '*.json'))):
    d = json.load(open(f)); es = os.path.basename(f).split('_')[0][3:]
    if 'upper' not in d: continue
    ok_h, Hlo = lower(d, Fr(es)); ok = ok_h and abs(Hlo - float(d['H_lo'])) < 1e-9 and Hlo <= float(d['H_lo']) + 1e-12
    allok &= ok; n += 1
    print('wapcone %-22s honest>=omega %s  H_lo >= %.6f (stored %.6f)  %s' % (os.path.basename(f), ok_h, Hlo, float(d['H_lo']), 'OK' if ok else 'FAIL'), flush=True)
for f in sorted(glob.glob(os.path.join(_H, '..', 'noise', 'certs', '*.json'))):
    d = json.load(open(f)); e = Fr(d['eta']); q = Fr(d['q']); model = d['model']
    if model == 'wd':   # white noise v and dephasing p together (the eight rows of Table S4); (v, p) from the file name
        v, pc = [Fr(x) for x in os.path.basename(f)[:-5].split('_q')[1].split('-')]
    else:
        v, pc = (q, 1) if model == 'white' else (1, q)
    ch = d['chsh']; c0, ge = loss_functional_F(G_CHSH_F, e); beta, _ = check_cert(ch['cert'], ge)
    Lloc = 4*e*(1-e) + 2*e*e
    S = c0 + v*(pc*beta + (1-pc)*Lloc)
    ok_s = S == Fr(ch['S_up']); Sa = aq(S)
    Hcap = up(-((arb(1)/2 + (arb(2) - Sa*Sa/4).sqrt()/2).log()/LOG2)) if S > 2 else 0.0
    line = 'noise %-28s S_up exact %s  every CHSH <= %.6f' % (os.path.basename(f), ok_s, Hcap); ok = ok_s
    for k in ('tilted', 'cone'):
        if k not in d or not d[k]: continue
        b = d[k]; ok_h, Hlo = lower(b, e, v, pc)
        okk = ok_h and (abs(Hlo - float(b['H_lo'])) < 1e-9 if Hlo > 0 else float(b['H_lo']) <= 1e-9)  # no claim if P_up >= 1
        ok &= okk; line += '  %s >= %.6f%s' % (k, Hlo, '' if okk else ' FAIL')
    allok &= ok; n += 1
    print(line + ('  OK' if ok else '  FAIL'), flush=True)
print('ALL OK =', allok, '(%d files)' % n)
