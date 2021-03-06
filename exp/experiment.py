import pandas as pd
import numpy as np
from scipy.special import comb

import sys
sys.path.insert(0, '../src')
import generate_k as genk
import genesetdp

import genesetmc as mc

if len(sys.argv)>1:
    Q=int(sys.argv[1])
else:
    Q=25

## Especify pathway here
pathway_name = 'GLYCOLYSIS_/_GLUCONEOGENESIS_-_HOMO_SAPIENS_(HUMAN)'

## example path
path = "~/git/binox/example/"

## read original patway information file
pathways_file = 'Pathways.tsv'
pathways = pd.read_csv(path + pathways_file, sep='\t', header=None)
pathways.columns = ['gene','pathway']

pathway_genes = pathways.where(pathways.pathway == pathway_name).dropna().gene.values

k = genk.generate_k(pathway_genes)

N = genesetdp.genesetdp(k,Q)
#assert(float(sum(N))==float(comb(sum(k),Q,exact=False)))

p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)

n_runs = 100000 ## number of monte carlo samples
p_vals_mc = mc.genesetmc(k, Q, n_runs)

groups = pd.read_csv('example-random/Random_Signatures.tsv', sep='\t')

network = genk.network

linked_genes = network.index[network.gene2.isin(pathway_genes)].tolist()

results = pd.DataFrame(columns = ['n','p','p_mc'])
## start for
for i in range(1,10001):
    groupname = 'group' + str(i)
    genes = groups.loc[groups.GroupName==groupname,'#Node']
    n = sum(linked_genes.count(x) for x in genes)
    p = p_vals[n]
    p_mc = p_vals_mc[n]
    results.loc[groupname] = [n,p, p_mc]

results.to_csv('example-random/genesetDP_results_Q' + str(Q) + '.tsv', sep = '\t')
