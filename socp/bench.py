"""Benchmark: SOCP (NDWU dual) vs closed forms from the literature vs NPA levels 1+AB, 2, 3 (Clarabel), plus
an explicit two-qubit lower bound. Output: bench.json and bench.md. Numerical (floating point) benchmark."""
import json, time, warnings, numpy as np
warnings.filterwarnings('ignore')
from socp_core import hQ_socp, hQ_qubit, hL
from npa import npa_value

def gigena_root(al, be):
    t0 = al**6*(27*be**2-8)+al**4*(-54*be**4+48*be**2+32)+al**2*(27*be**6+48*be**4-400*be**2+128)-8*be**6+32*be**4+128*be**2-512
    t1 = -168*al**3*be**3+60*al**5*be+160*al**3*be+60*al*be**5+160*al*be**3-576*al*be
    t2 = -6*al**4*be**2-6*al**2*be**4-64*al**2*be**2+20*al**4+96*al**2+20*be**4+96*be**2+320
    t3 = 8*al**3*be**3-24*al**3*be-24*al*be**3+96*al*be
    t4 = 11*al**2*be**2-16*al**2-16*be**2-64
    t5 = 4*al*be; t6 = 4
    r = np.roots([t6, t5, t4, t3, t2, t1, t0]); r = r[np.abs(r.imag) < 1e-9].real
    return float(r.max())
def gigena_quartic(al):
    r = np.roots([1, 4-al**2, 11/4*al**4-12*al**2-4, 2*al**6-al**4-20*al**2-32, 5*al**6-21*al**4+16*al**2-32])
    r = r[np.abs(r.imag) < 1e-9].real; return float(r.max())

cases = []
def add(fam, label, g, cf, src): cases.append(dict(family=fam, label=label, g=list(map(float, g)), closed=cf, source=src))
add('CHSH', 'CHSH', [0,0,0,0,1,1,1,-1], 2*np.sqrt(2), 'Tsirelson')
for a in (0.5, 1.0, 1.5):
    add('tilted CHSH [AMP]', f'alpha={a}', [a,0,0,0,1,1,1,-1], np.sqrt(8+2*a*a), 'sqrt(8+2a^2)')
for (al, be) in [(0.5,0.5),(0.9,0.9),(0.99,0.99),(0.999,0.999),(1.0,0.5),(0.3,1.2),(1.5,0.3),(1.9,0.05)]:
    cf = gigena_quartic(al) if al == be else None
    add('doubly tilted CHSH [Gigena]', f'alpha={al},beta={be}', [al,0,be,0,1,1,1,-1], cf, 'quartic Eq.(18) (alpha=beta); no usable closed form otherwise')
    cases[-1]['printed_sextic_largest_root'] = gigena_root(al, be)
for (a, ph) in [(0.8, np.pi/4), (1.2, 1.0), (1.6, 0.5)]:
    c2 = 1+np.cos(2*ph)
    y2 = a*a*np.sin(ph)*(96-16*a*a-a**4*c2)/((4-a*a)*(32-16*a*a+a**4*c2))
    cf = np.sqrt(2*(4-a*a)*(32-a**4*c2)/(32-16*a*a+a**4*c2)) if abs(y2) <= 1 else None
    add('MK double-tilted', f'alpha={a},phi={ph:.4f}', [a*np.cos(ph/2), a*np.sin(ph/2),0,0,1,1,1,-1], cf, 'Mikos-Nuszkiewicz-Kaniewski Eq.(39)' + ('' if cf else ' (|y2|>1: formula not applicable)'))
for (a0, a1) in [(0.2, 0.5), (-0.3, 0.8), (0.5, 1.0)]:
    xp = (a0*a1+np.sqrt((a1*a1+4-4*a0*a0)/(2-a0*a0)))/(2*(1-a0*a0))
    lp = (a0*a1+np.sqrt((4+a1*a1-4*a0*a0)*(2-a0*a0)))/(1-a0*a0)
    ok = (a1/2 < xp < 1)
    add('generalised Wolfe-Yelin [MK]', f'alpha0={a0},alpha1={a1}', [a0,a0,a1,0,1,1,1,-1], lp if ok else None, 'MK Eq.(54)-(57)' + ('' if ok else ' (conditions (56) fail)'))
for ph in (np.pi/6, 0.3):
    add('WBC I_phi', f'phi={ph:.4f}', [0,0,0,0,1,1/np.sin(ph),1/np.sin(ph),-1/np.cos(2*ph)], 2*np.cos(ph)**3/(np.cos(2*ph)*np.sin(ph)), 'WBC22')
for gm in (0.0, 0.1):
    mu = np.cos(gm+np.pi/6); k = 4*mu*mu-1
    add('WBC J_gamma', f'gamma={gm}', [0,0,0,0,1, k, k, -k], 8*mu**3, 'WBC22')
for (th, ph) in [(np.pi/8, np.pi/8), (0.3, 0.2), (np.pi/8, -0.3)]:
    t2 = 1/(np.sin(2*th)**2/np.tan(ph)**2 - np.cos(2*th)**2)
    g = [0,0,t2*np.cos(2*th)/np.cos(ph),t2*np.cos(2*th)/np.cos(ph),1/np.cos(ph),1/np.cos(ph),t2*np.sin(2*th)/np.sin(ph),-t2*np.sin(2*th)/np.sin(ph)]
    add('extended tilted CHSH [BSB, MPW]', f'theta={th:.4f},phi={ph:.4f}', g, 2*(1+t2), 'eta_Q = 2(1+tau^2)')
for a in (0.5, 2.0):
    add('asymmetric CHSH [WAP]', f'a={a}', [0,0,0,0,1,-1,a,a], 2*np.sqrt(1+a*a), '2 sqrt(1+a^2)')
rho = 1-1/np.sqrt(2)
add('BB octagon vertex g0', 'g0', [rho/np.sqrt(2),rho/np.sqrt(2),-rho,0]+[1/(2*np.sqrt(2))]*3+[-1/(2*np.sqrt(2))], 1.0, 'Barizien-Bancal')

out = []
for cs in cases:
    g = cs['g']; r = hQ_socp(g)
    cs['socp'] = r['value']; cs['socp_c'] = float(r['c']); cs['socp_time_s'] = r['time']
    t = time.time(); cs['qubit_lower'] = hQ_qubit(g, starts=16); cs['qubit_time_s'] = time.time()-t
    cs['local'] = hL(g)
    for lev in ('1+AB', 2, 3):
        t = time.time()
        try: v = npa_value(g, lev)
        except Exception as e: v = None
        cs[f'npa_{lev}'] = v; cs[f'npa_{lev}_time_s'] = time.time()-t
    out.append(cs); print(json.dumps({k: cs[k] for k in ('family','label','closed','socp','qubit_lower','npa_1+AB','npa_2','npa_3','socp_time_s')}), flush=True)
json.dump(out, open('bench.json','w'), indent=1, default=float)
