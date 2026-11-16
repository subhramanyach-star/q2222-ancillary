"""Ideal targets for global randomness, 2 x 3 scenario. State cos(th)|00>+sin(th)|11>, A0=sz, A1=sx,
B_{0,1} = cos(mu) sz +- sin(mu) sx with cos(mu)=1/r (tilted-CHSH optimal; th=pi/4 gives CHSH Tsirelson),
key setting B2 = cos(phi) sz + sin(phi) sx.  z = (A0,A1,B0,B1,B2, A0B0,A0B1,A0B2, A1B0,A1B1,A1B2)."""
import numpy as np
from families import ptheta, gtilt
def target(th, phi):
    p, c0, (s, C, r) = ptheta(th)
    B2 = C * np.cos(phi); A0B2 = np.cos(phi); A1B2 = s * np.sin(phi)
    return np.array([p[0], p[1], p[2], p[3], B2, p[4], p[5], A0B2, p[6], p[7], A1B2])
def lift(g8):
    """2222 functional (a0,a1,b0,b1,E00,E01,E10,E11) -> 2x3 coordinates."""
    g = np.zeros(11); g[[0, 1, 2, 3, 5, 6, 8, 9]] = g8; return g
def pab(z, xs, ys, nB=3):
    a = z[xs]; b = z[2 + ys]; e = z[2 + nB + nB * xs + ys]
    return np.array([(1 + sa * a + sb * b + sa * sb * e) / 4 for sa in (1, -1) for sb in (1, -1)])
def H(p): p = p[p > 1e-15]; return float(-(p * np.log2(p)).sum())
def best_key(th, v=1.0):
    """(xs, phi) maximising the Shannon entropy of p(ab|xs,2) at visibility v (grid search)."""
    best = None
    for xs in (0, 1):
        for phi in np.linspace(0, np.pi, 721):
            h = H(pab(v * target(th, phi), xs, 2))
            if best is None or h > best[0] + 1e-12: best = (h, xs, phi)
    return best
if __name__ == '__main__':
    for lab, th in [('0.1pi', .1*np.pi), ('pi/8', np.pi/8), ('27.76', 27.76/180*np.pi), ('3pi/16', 3*np.pi/16), ('0.22pi', .22*np.pi), ('pi/4', np.pi/4)]:
        h, xs, phi = best_key(th); print(lab, 'ideal H(AB)=%.4f xs=%d phi/pi=%.4f' % (h, xs, phi / np.pi))
