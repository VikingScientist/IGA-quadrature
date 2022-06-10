from getOptimalQuadPoints import getOptimalQuadPoints
import numpy as np

def get_knot(p,k,nel):
    ans = []
    for i in range(nel):
        ans += [i] * (p-k)
    return ans

nel = 128
for p in range(1, 14):
    for k in range(p-1):
        knot = get_knot(p,k,nel)
        w,x,rec,it = getOptimalQuadPoints(knot,p)
        if rec == 1:
            print(f'{it+1:3d}', end='')
            if((len(knot)-p-1) % 2 == 1): print('* ', end='')
            else: print('  ', end='')
        else:
            print('  -  ', end='')
    print('')

