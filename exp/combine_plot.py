import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


binox = pd.read_csv('example-random/BinoX_results.tsv', sep='\t')

###### only enriched
binox.loc[binox['#6:relationType'] == '-', '#4:p.value'] = 1 - binox.loc[binox['#6:relationType'] == '-', '#4:p.value']
######



binox = binox.set_index('#3:NameGroupB')

binox.index = [x.lower() for x in binox.index]

binox = binox['#4:p.value']

genesetdp = pd.read_csv('example-random/genesetDP_results.tsv', sep='\t', index_col=0)

binox = binox.loc[genesetdp.index]
binox.loc[np.isnan(binox)] = 1

binox_sorted = binox.sort_values()
binox_sorted.index = (np.arange(len(binox_sorted.index))+1)/len(binox_sorted.index)

genesetdp_sorted = genesetdp['p'].sort_values()
genesetdp_sorted.index = (np.arange(len(genesetdp_sorted.index))+1)/len(genesetdp_sorted)


x_y = pd.DataFrame(np.arange(0,1,0.001), index = np.arange(0,1,0.001))
x_2y = pd.DataFrame(np.arange(0,1,0.001)*2, index = np.arange(0,1,0.001))
x_halfy = pd.DataFrame(np.arange(0,0.5,0.001), index = np.arange(0,0.5,0.001)*2)

sns.set_style("ticks")
sns.despine()


plt.loglog(binox_sorted, '.')
plt.loglog(genesetdp_sorted, '.')
plt.loglog(x_y, '--', color = 'r', alpha = 0.5)
plt.loglog(x_2y, ':',color = 'b', alpha = 0.5)
plt.loglog(x_halfy, ':',color = 'b', alpha = 0.5)

plt.legend(['BinoX', 'GeneSetDP', 'y=x', 'y=0.5x,y=2x'])
plt.ylabel('$p$ value')
plt.xlabel('Normalized rank')

plt.savefig('calibration.png')
