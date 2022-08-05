from getOptimalQuadPoints import getOptimalQuadPoints
import numpy as np
import sys

def get_knot(p,k,nel):
    ans = [0] * (k+1)
    for i in range(nel+1):
        ans += [i] * (p-k)
    ans += [nel] * (k+1)
    return ans

nel = 128
for p in range(8,17):
    for k in range(p):
        knot = get_knot(p,k,nel)
        w,x,rec,it = getOptimalQuadPoints(knot,p)
        print(f'{rec:3d}  ', end='')
    print('')

