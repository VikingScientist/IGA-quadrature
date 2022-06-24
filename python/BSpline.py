from splipy import BSplineBasis
from splipy.state import state

def BSpline(knot, p, t):
    # augment knot vector so we always have open knot vector
    ext_knot = [knot[0]]*p + list(knot) + [knot[-1]]*p
    with state(knot_tolerance=1e-18):
        b  = BSplineBasis(p+1, ext_knot)
        N  = b(t, sparse=True).T
        dN = b(t,1, sparse=True).T
    # return N, dN
    return N[p:-p,:], dN[p:-p,:]
