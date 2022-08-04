import numpy as np
from BSpline import BSpline
from scipy import sparse
import warnings

def getOptimalQuadPoints(knot, p, **kwargs):
# intended use: w, x = getOptimalQuadPoints(knot,p)
    n  = len(knot)-p-1  # dimension of our spline space
    
    if n % 2==1: # need to have a space of even dimension
        i = np.where(np.diff(knot) == np.max(np.diff(knot)))  
        i = i[0][len(i[0])//2]                               # insert new knot in middle of
        knot = sorted(knot + [np.mean(knot[i:i+2])]) # the largest, centermost knot span
        n += 1 

    # compute all greville points and integrals (used for initial guess)
    greville       = np.zeros(n)
    exact_integral = np.zeros(n)
    for i in range(n):
        greville[i]       = np.sum(knot[i+1:i+p+1]) / p 
        exact_integral[i] = (knot[i+p+1]-knot[i])/(p+1) 

    if 'x0' in kwargs: # if initial guess is provided, use these
        x = kwargs['x0']
        w = kwargs['w0']
    else:           # else compute them based on greville points and integrals
        w   = (exact_integral[0::2] + exact_integral[1::2])   
        x   = (      greville[0::2] +       greville[1::2])/2 
    # counter variable to count the number of recursive calls
    rec = kwargs['rec'] if 'rec' in kwargs else 1

    newton_tol       = 1e-10           # convergence tolerance
    newton_max_it    = 15              # max iterations before divergence
    while True:                        # recursive loop from algorithm 2
        for it in range(newton_max_it):  # newton iteration loop
            N, dN = BSpline(knot, p, x)
            F  = N @ w - exact_integral 

            ### This enumeration of the degrees of freedom crashes with umfpack solvers
            # dF = sparse.hstack((N, dN @ sparse.diags(w)))

            ### We need interleaved dofs. This is a readable version of what we want
            # dF = sparse.csr_matrix((n, n))
            # dF[:, ::2]  = N
            # dF[:, 1::2] = dN @ sparse.diags(w)

            ### Stupid sparse matrices need unreadable code. The below is just a more efficient
            ### implementaiton of the three lines above
            dN = dN @ sparse.diags(w)
            Ni = sorted(list(N.indptr) + list(N.indptr[1:]))
            N.resize(n,n)
            N.indptr = np.array(Ni)
            dNi = sorted(list(dN.indptr) + list(dN.indptr[:-1]))
            dN.resize(n,n)
            dN.indptr = np.array(dNi)
            dF = N + dN


            warnings.filterwarnings("error")
            try:
                dx = sparse.linalg.spsolve(dF,-F, use_umfpack=True)
            except Warning: # singular matrix
                warnings.filterwarnings("ignore")
                break
            except RuntimeError: # singular matrix
                warnings.filterwarnings("ignore")
                break
            warnings.filterwarnings("ignore")
            ### umfpack solvers have problems with sequential enumeration (matlab code uses this)
            # w = w + dx[0:n//2]
            # x = x + dx[n//2: ]
            ### interleaved dofs work much better
            w = w + dx[0::2]
            x = x + dx[1::2]

            # test for diverging (coarse heuristic, see section 3.3)
            if( np.min(x)<knot[ 0]): break
            if( np.max(x)>knot[-1]): break

            # test for converging 
            if(np.linalg.norm(dx)<newton_tol):  return (w, x, rec, it)

        # at this point, newton iteration has diverged. solve recursively on easier knot
        if 'knot0' in kwargs:
            w, x, rec, it = getOptimalQuadPoints((kwargs['knot0'] + knot)/2, p, **kwargs)
            kwargs['rec'] = rec + 1
            kwargs['knot0'] = (kwargs['knot0'] + knot)/2 
        else:
            uniformKnot = np.linspace(knot[0],knot[-1], n+p+1) 
            w, x, rec, it = getOptimalQuadPoints(uniformKnot, p, rec=1) 
            kwargs['rec'] = rec
            kwargs['knot0'] = uniformKnot 
        kwargs['rec'] += 1 
        kwargs['x0'] = x 
        kwargs['w0'] = w 
    # end loop up and start newton iteration with better initial guess
