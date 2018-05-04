import pandas as pd
import numpy as np
from genesetdp import find_maxscore

def genesetmc(k,q, n_runs):

    maxscore = find_maxscore(k,q)
    N_mc = np.zeros(maxscore+10)

    for i in range(n_runs):
        sample = sum(np.random.choice(np.arange(len(k)), q, p = np.array(k)/sum(k)))
        N_mc[sample] += 1


    p_vals = (np.flipud(np.cumsum(np.flipud(N_mc))) - N_mc/2)/sum(N_mc)

    return p_vals
