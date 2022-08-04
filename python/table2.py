from getOptimalQuadPoints import getOptimalQuadPoints
import numpy as np
import sys

def get_knot(p,k,nel):
    ans = [0] * (k+1)
    for i in range(nel+1):
        ans += [i] * (p-k)
    ans += [nel] * (k+1)
    return ans
sys.setrecursionlimit(25) # only interested in the cases that converge on first try

nel = 128
for p in range(1,16):
    for k in range(p):
        knot = get_knot(p,k,nel)
        try:
            w,x,rec,it = getOptimalQuadPoints(knot,p)
            if rec == 1: print(f'{it+1:3d}  ', end='')
            else:        print('  -  ',        end='')
        except RecursionError:
            print('  -  ', end='')
    print('')

