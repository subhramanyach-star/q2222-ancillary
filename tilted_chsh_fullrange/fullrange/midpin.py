"""Kink candidates of the endpoint function D -> Phi_eps at c = eps (eps = +-1); works for numpy floats and arb."""
def cands(g, sg):
    a0, a1 = g[0], g[1]; bp = (g[2] + sg * g[3]) / 2
    gam0 = g[4] + sg * g[5]; gam1 = g[6] + sg * g[7]
    p0, q0, p1, q1 = bp + gam0, bp - gam0, bp + gam1, bp - gam1
    K0 = [-p0, -q0, -(p0 + q0) / 2 - a0, -(p0 + q0) / 2 + a0]
    K1 = [p1, q1, (p1 + q1) / 2 + a1, (p1 + q1) / 2 - a1]
    return K0 + K1
