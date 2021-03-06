import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr

def one_plot(ax, Q):

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

    combined = pd.merge(binox, genesetdp, right_index = True, left_index = True)
    combined = combined[['p','enrichment']]
    
    print(spearmanr(combined['enrichment'], combined['p']))


    x_y = pd.DataFrame(np.arange(0,1,0.001), index = np.arange(0,1,0.001))
    x_2y = pd.DataFrame(np.arange(0,1,0.001)*2, index = np.arange(0,1,0.001))
    x_halfy = pd.DataFrame(np.arange(0,0.5,0.001), index = np.arange(0,0.5,0.001)*2)

    sns.set_style("ticks")
    sns.despine()

    ax.loglog(combined['enrichment'], combined['p'],'.')
    # ax.loglog(binox_ix, binox_both_sorted, '.')

    ax.loglog(x_y, '--', color = 'r', alpha = 0.5)
    ax.loglog(x_2y, ':',color = 'b', alpha = 0.5)
    ax.loglog(x_halfy, ':',color = 'b', alpha = 0.5)


fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True)

one_plot(ax1, 10)
one_plot(ax2, 20)
one_plot(ax3, 30)

ax1.set_title('Q = 10', fontsize=18)
ax2.set_title('Q = 20', fontsize=18)
ax3.set_title('Q = 30', fontsize=18)

ax1.set_ylabel('GeneSetDP', fontsize = 16)
ax1.set_xlabel('BinoX', fontsize = 16)
ax2.set_xlabel('BinoX', fontsize = 16)
ax3.set_xlabel('BinoX', fontsize = 16)

for ax in [ax1,ax2,ax3]:
    ax.xaxis.set_tick_params(labelsize=12)
    ax.yaxis.set_tick_params(labelsize=12)

# lgd = fig.legend(['GeneSetDP', 'GeneSetMC', 'BinoX'], loc='center', ncol =3, bbox_to_anchor=(0.5,0), fontsize=16)
# fig.legend(['GeneSetDP', 'GeneSetMC', 'BinoX'], fontsize=15, ncol=3)


# fig.set_size_inches(15,4)
fig.savefig('rank_multiple.png', bbox_inches='tight')
# plt.show()

