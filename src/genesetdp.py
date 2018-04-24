#/usr/bin python
import numpy as np
import operator as op
import functools

def ncr(n, r):
    r = min(r, n-r)
    numer = functools.reduce(op.mul, range(n, n-r, -1), 1)
    denom = functools.reduce(op.mul, range(1, r+1), 1)
    return numer//denom

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
    N = np.zeros((max_s+1,q+1),dtype='u4')
    N[0,0] = 1

    for a in range(len(k)):
        for s in range(max_s,-1,-1):
            for c in range(q,-1,-1):
                for b in range(k[a],0,-1): # Skipping zero (leftover)
                    if c-b>=0 and s-a*b>=0:
                        #print("pre\n",N,"a=",a,"s=",s,"c=",c,"b=",b,"k[a]=",k[a],"N[s-a*b,c-b]=",N[s-a*b,c-b],s-a*b,c-b)
                        N[s,c] += ncr(k[a], b) * N[s-a*b,c-b]
                        #print("post\n",N,"a=",a,"s=",s,"c=",c,"b=",b,"k[a]=",k[a])
                    #else:
                        #print("skip", N,"a=",a,"s=",s,"c=",c,"b=",b,"k[a]=",k[a],s-a*b,c-b)
    return N[:,q]

if __name__ == "__main__":
    # Query size q, number of links k[a]
    q = 8
    k = np.array([8,5,3,3,4],dtype='u4')
    #k = np.array([3,3,3],dtype='u4')
    print(genesetdp(k,q))
