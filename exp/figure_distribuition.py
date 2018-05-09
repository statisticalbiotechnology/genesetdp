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

N10 = genesetdp.genesetdp(k,10)
N20 = genesetdp.genesetdp(k,20)
N30 = genesetdp.genesetdp(k,30)

####

import matplotlib.pyplot as plt
import seaborn as sns

fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)

sns.set_style("ticks")
sns.despine()

ax1.set_ylabel('$N(s)$', fontsize=16)
ax2.set_ylabel('$N(s)$', fontsize=16)
ax3.set_ylabel('$N(s)$', fontsize=16)
ax3.set_xlabel('Score, $s$', fontsize=16)

ax1.set_title('                                 $Q=10$', fontsize=14, horizontalalignment='left')
ax2.set_title('                                 $Q=20$', fontsize=14, horizontalalignment='left')
ax3.set_title('                                 $Q=30$', fontsize=14, horizontalalignment='left')

for ax in [ax1,ax2,ax3]:
    ax.xaxis.set_tick_params(labelsize=12)
    ax.yaxis.set_tick_params(labelsize=12)


ax1.step(np.arange(len(N10)), N10, where='mid')
ax1.plot(len(N10),1, '.', color='r', markersize=12)
ax2.step(np.arange(len(N20)), N20, where='mid')
ax2.plot(len(N20),1, '.', color='r', markersize=12)
ax3.step(np.arange(len(N30)), N30, where='mid')
ax3.plot(len(N30),1, '.', color='r', markersize=12)

plt.tight_layout()
plt.savefig('score_distribuition_multiple.png')
# plt.show()



