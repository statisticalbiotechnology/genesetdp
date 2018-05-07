import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

Q = 30

binox = pd.read_csv('example-random/BinoX_results_Q' + str(Q) + '.tsv', sep='\t')

###### only enriched
binox.loc[:, 'enrichment'] = binox.loc[:, '#4:p.value']
binox.loc[binox['#6:relationType'] == '-', 'enrichment'] = 1 - binox.loc[binox['#6:relationType'] == '-', '#4:p.value']
######

binox = binox.set_index('#3:NameGroupB')

binox.index = [x.lower() for x in binox.index]

binox_both = binox['#4:p.value']
binox_enrich = binox['enrichment']

genesetdp = pd.read_csv('example-random/genesetDP_results_Q' + str(Q) + '.tsv', sep='\t', index_col=0)
genesetdp_sorted = genesetdp['p'].sort_values()
genesetdp_sorted_mc = genesetdp['p_mc'].sort_values()
genesetdp_ix = (np.arange(len(genesetdp_sorted.index))+0.5)/float(len(genesetdp_sorted))

binox_both = binox_both.loc[genesetdp.index]
binox_enrich = binox_enrich.loc[genesetdp.index]
binox_both.loc[np.isnan(binox_both)] = 1
binox_enrich.loc[np.isnan(binox_enrich)] = 1

binox_both_sorted = binox_both.sort_values()
binox_enrich_sorted = binox_enrich.sort_values()

binox_ix = (np.arange(len(binox_both_sorted.index))+0.5)/float(len(binox_both_sorted.index))

x_y = pd.DataFrame(np.arange(0,1,0.001), index = np.arange(0,1,0.001))
x_2y = pd.DataFrame(np.arange(0,1,0.001)*2, index = np.arange(0,1,0.001))
x_halfy = pd.DataFrame(np.arange(0,0.5,0.001), index = np.arange(0,0.5,0.001)*2)

sns.set_style("ticks")
sns.despine()

plt.rc('font', size=12)
plt.loglog(binox_ix, binox_both_sorted, '.')
plt.loglog(binox_ix, binox_enrich_sorted, '.')
plt.loglog(genesetdp_ix, genesetdp_sorted,'.')
plt.loglog(genesetdp_ix, genesetdp_sorted_mc,'.')
plt.loglog(x_y, '--', color = 'r', alpha = 0.5)
plt.loglog(x_2y, ':',color = 'b', alpha = 0.5)
plt.loglog(x_halfy, ':',color = 'b', alpha = 0.5)

plt.ylabel('$p$ value')
plt.xlabel('Normalized rank')

plt.savefig('calibration_Q' + str(Q) + '.png')
plt.show()
