import pandas as pd
import numpy as np
from genesetdp import find_maxscore

def genesetmc(k,q, n_runs):

    maxscore = find_maxscore(k,q)
    N_mc = np.zeros(maxscore+1)
    gene_links = np.array([ix for ix in range(len(k)) for _ in range(k[ix])])

    for i in range(n_runs):
        sample = sum(np.random.choice(gene_links, q, replace=False))
        N_mc[sample] += 1.0


    p_vals = (np.flipud(np.cumsum(np.flipud(N_mc))) - N_mc/2)/sum(N_mc)

    return p_vals
