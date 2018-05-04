import sys
import pandas as pd
import numpy as np

if len(sys.argv)>1:
    Q=int(sys.argv[1])
else:
    Q=25

original = pd.read_csv('~/git/binox/example/Network.tsv', sep='\t')
## manipulate network
network_treshold = 0.7
network_cut = original['#0:PFC'] >= network_treshold
network_one = original.loc[network_cut]
network = network_one.iloc[:,[0,1]].reset_index(drop=True)
network.columns = ['gene1','gene2']
network = network.set_index('gene1')

genes = np.unique(network.index)
print("Num genes=",len(genes))

random_signatures = pd.DataFrame()

n_samples = 10000
n_genes = Q

for i in np.arange(1,n_samples+1,1):
    sample_genes = np.random.choice(genes, n_genes, replace=False)
    new_group = pd.DataFrame(np.transpose([sample_genes, np.repeat('group' + str(i), n_genes)]))
    random_signatures = random_signatures.append(new_group)

random_signatures.columns = ['#Node', 'GroupName']


random_signatures.set_index('#Node').to_csv('example-random/Random_Signatures.tsv', sep = '\t')


#####

pathways = pd.read_csv('~/git/binox/example/Pathways.tsv', sep='\t', header=None)

pathway_name = 'GLYCOLYSIS_/_GLUCONEOGENESIS_-_HOMO_SAPIENS_(HUMAN)'

one_pathway = pathways.where(pathways.iloc[:,1] == pathway_name).dropna()

one_pathway.to_csv('example-random/Pathway_One.tsv', sep='\t', header=None, index=False)
