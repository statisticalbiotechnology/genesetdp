import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns


def one_plot(ax, Q):

    genesetdp = pd.read_csv('example-random/genesetDP_results_Q' + str(Q) + '.tsv', sep='\t', index_col=0)
    sns.set_style("ticks")
    sns.despine()

    dp = genesetdp['p']
    mc = genesetdp['p_mc']


    ax.loglog(mc,dp,'.')

    x_y = pd.DataFrame(np.arange(0,1,0.001), index = np.arange(0,1,0.001))
    x_2y = pd.DataFrame(np.arange(0,1,0.001)*2, index = np.arange(0,1,0.001))
    x_halfy = pd.DataFrame(np.arange(0,0.5,0.001), index = np.arange(0,0.5,0.001)*2)

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

ax1.set_ylabel('GeneSetDP $p$-value', fontsize = 16)
ax1.set_xlabel('GeneSetMC $p$-value', fontsize = 16)
ax2.set_xlabel('GeneSetMC $p$-value', fontsize = 16)
ax3.set_xlabel('GeneSetMC $p$-value', fontsize = 16)

for ax in [ax1,ax2,ax3]:
    ax.xaxis.set_tick_params(labelsize=12)
    ax.yaxis.set_tick_params(labelsize=12)

fig.set_size_inches(15,4)
fig.savefig('scatter_dp_mc.png', bbox_inches='tight')

