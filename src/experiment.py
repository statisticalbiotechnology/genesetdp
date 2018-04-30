import pandas as pd
import numpy as np

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

N = genesetdp.genesetdp(k,10)

p = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)

import matplotlib.pyplot as plt
plt.step(np.arange(len(N))+1, N, where='mid')
plt.xlabel('score')
plt.ylabel('N(s)')

plt.savefig('F1_score_distribuition.png')

