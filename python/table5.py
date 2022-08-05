from getOptimalQuadPoints import getOptimalQuadPoints
import numpy as np
import sys

def get_knot(p,k,nel):
    alpha = 9/10
    ans = [1] * (k+1)
    for i in range(nel+1):
        ans += [alpha**i] * (p-k)
    ans += [ans[-1]] * (k+1)
    return ans[::-1]

nel = 64
for p in range(8,15):
    for k in range(p):
        knot = get_knot(p,k,nel)
        w,x,rec,it = getOptimalQuadPoints(knot,p)
        print(f'{rec:3d}  ', end='')
    print('')

