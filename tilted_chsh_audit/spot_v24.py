"""Floating-point spot check of the hand certificates for V2, V4 (Appendix C), theta in (0, pi/4)."""
import numpy as np
def Phi(g, d, c):
    G = np.array([[1, c], [c, 1]]); b = g[2:4]; e = [g[4:6], g[6:8]]; B = [b/2 + d, b/2 - d]
    return sum(max(a*g[x] + np.sqrt(max(0.0, (B[x] + a*e[x]) @ G @ (B[x] + a*e[x]))) for a in (1, -1)) for x in (0, 1))
worst = -1; eqpts = []
for th in np.linspace(0.01, np.pi/4 - 0.01, 60):
    s = np.sin(2*th); C = np.cos(2*th); r = np.sqrt(1 + s*s); D = r + 1; c0 = (1 - s*s)/(1 + s*s)
    for which in ('V2', 'V4'):
        if which == 'V2':
            a0 = (r + C)/(2*D); be = (C - 1)/(2*D); ka = (C + r + 2*s*s)/(4*D); nu = (r + 2 - C)/(4*D)
        else:
            a0 = (C - r)/(2*D); be = (C + 1)/(2*D); ka = (-C + r + 2*s*s)/(4*D); nu = (r + 2 + C)/(4*D)
        g = np.array([a0, 0, be, be, ka, ka, nu, -nu])
        P = np.array([C, 0, C/r, C/r, 1/r, 1/r, s*s/r, -s*s/r]); assert abs(g @ P - 1) < 1e-12
        u1 = abs(a0)/ka; c1 = u1*u1/2 - 1
        for c in np.concatenate([np.linspace(-1, 1, 2001), [c0]]):
            u = np.sqrt(2 + 2*c)
            tau = -be/2 - a0/u if c >= c1 else -be/2 - np.sign(a0)*ka
            ph = Phi(g, np.array([tau, tau]), c)
            worst = max(worst, ph)
            if ph > 1 - 1e-9: eqpts.append(round(c - c0, 6) if abs(c - c0) < 1e-3 else c)
print('max Phi over grid (theta in (0,pi/4), both vertices):', worst)
print('near-equality points (c, or c-c0 near c0):', sorted(set(eqpts))[:10], '...')
