import pandas as pd
import numpy as np
from scipy.special import comb

import sys
sys.path.insert(0, '../src')
import generate_k as genk
import genesetdp

if len(sys.argv)>1:
    Q=int(sys.argv[1])
else:
    Q=25

## example path
path = "~/git/binox/example/"

## read original patway information file
pathways_file = 'Pathways.tsv'
pathways = pd.read_csv(path + pathways_file, sep='\t', header=None)
pathways.columns = ['gene','pathway']

pathway_name = 'GLYCOLYSIS_/_GLUCONEOGENESIS_-_HOMO_SAPIENS_(HUMAN)'

pathway_genes = pathways.where(pathways.pathway == pathway_name).dropna().gene.values

k = genk.generate_k(pathway_genes)

N = genesetdp.genesetdp(k,Q)
#assert(float(sum(N))==float(genesetdp.ncr(sum(k),Q)))
assert(float(sum(N))==float(comb(sum(k),Q,exact=False)))

p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)

groups = pd.read_csv('example-random/Random_Signatures.tsv', sep='\t')

network = genk.network

# linked_genes = network.loc[pathway_genes].dropna().gene2.values.tolist()
linked_genes = network.index[network.gene2.isin(pathway_genes)].tolist()

results = pd.DataFrame(columns = ['n','p'])
## start for
for i in range(1,10001):
    groupname = 'group' + str(i)
    genes = groups.loc[groups.GroupName==groupname,'#Node']
    n = sum(linked_genes.count(x) for x in genes)
    p = p_vals[n]
    results.loc[groupname] = [n,p]

results.to_csv('example-random/genesetDP_results.tsv', sep = '\t')
