import pandas as pd
import numpy as np

import sys
sys.path.insert(0, '../src')
import generate_k as genk
import genesetdp

## example path
path = "~/git/binox/example/"

## read original patway information file
pathways_file = 'Pathways.tsv'
pathways = pd.read_csv(path + pathways_file, sep='\t', header=None)
pathways.columns = ['gene','pathway']

pathway_name = 'GLYCOLYSIS_/_GLUCONEOGENESIS_-_HOMO_SAPIENS_(HUMAN)'

pathway_genes = pathways.where(pathways.pathway == pathway_name).dropna().gene.values

k = genk.generate_k(pathway_genes)

N = genesetdp.genesetdp(k,25)

p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)

####

import matplotlib.pyplot as plt
import seaborn as sns

sns.set_style("ticks")
sns.despine()

plt.step(np.arange(len(N)), N, where='mid')
plt.ylabel('$N(s)$')
plt.xlabel('Score, $s$')

plt.savefig('score_distribuition.png')


####

groups = pd.read_csv('../exp/example-random/Random_Signatures.tsv', sep='\t')

network = genk.network

linked_genes = network.loc[pathway_genes].dropna().gene2.values.tolist()

results = pd.DataFrame(columns = ['n','p'])
## start for
for i in range(1,10001):
    groupname = 'group' + str(i)
    genes = groups.loc[groups.GroupName==groupname,'#Node']
    n = sum(linked_genes.count(x) for x in genes)
    p = p_vals[n]
    results.loc[groupname] = [n,p]

results.to_csv('genesetDP_results.tsv', sep = '\t')
