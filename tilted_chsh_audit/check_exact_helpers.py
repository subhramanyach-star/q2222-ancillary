import model as M
from kfield import K
def vertices_and_mu():
    V = {}
    V['V1'] = M.vertex_from_tex('V_1'); V['V3'] = M.vertex_from_tex('V_3')
    V2, V4 = M.V2V4(); V['V2'] = V2; V['V4'] = V4
    V['V1s'] = dict(g=M.sigma_g(V['V1']['g'])); V['V3s'] = dict(g=M.sigma_g(V['V3']['g']))
    MU = {}
    for name, W in V.items():
        g = W['g']; mu = {}
        for x in (0, 1):
            vp, vm = M.VV[(1, x)], M.VV[(-1, x)]
            for comp in (0, 1):
                det = vp[comp] - vm[comp]
                if not det.iszero(): break
            r1, r2 = -g[x], g[4 + 2*x + comp]
            mu[(1, x)] = (r1*(-vm[comp]) + r2) / det; mu[(-1, x)] = (r2 - vp[comp]*r1) / det
        MU[name] = mu
    return V, MU
