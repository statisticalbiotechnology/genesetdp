#/usr/bin python
import numpy as np
from scipy.special import comb

def find_maxscore(k,q):
    max_s, ix, rem = 0, len(k)-1, q
    while rem > k[ix]:
        max_s += k[ix]*ix
        rem -= k[ix]
        ix -= 1
    max_s += rem*ix
    return max_s

def genesetdp(k,q):
    max_s = find_maxscore(k,q)
    #print(max_s,q)
    N = np.zeros((max_s+1,q+1))
    N[0,0] = 1

    for a in range(len(k)):
        for s in range(max_s,-1,-1):
            for c in range(q,-1,-1):
                for b in range(k[a],0,-1):
                    if c-b>=0 and s-a*b>=0:
                        N[s,c] += comb(k[a], b, exact=False) * N[s-a*b,c-b]
    return N[:,q]

if __name__ == "__main__":
    # Query size q, number of links k[a]
    q = 2
    k = np.array([4,2,1,1],dtype='u4')
    #k = np.array([3,3,3],dtype='u4')
    N=genesetdp(k,q)
    p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)
    x,sk,ix=[],float(sum(k)),0
    for ki in k:
        for _ in range(ki):
            x.insert(0,(ix+0.5)/(sk))
            ix += 1
    print(N)
    print(p_vals)
    print(x)
