import pandas as pd
import numpy as np

original = pd.read_csv('Signatures.tsv', sep='\t')

genes = np.unique(original.index)

random_signatures = pd.DataFrame()

n_samples = 10000
n_genes = 25

for i in np.arange(1,n_samples+1,1):
    sample = np.random.choice(len(genes),n_genes, replace=False)
    sample_genes = genes[sample]
    new_group = pd.DataFrame(np.transpose([sample_genes, np.repeat('group' + str(i), n_genes)]))
    random_signatures = random_signatures.append(new_group)

random_signatures.columns = ['#Node', 'GroupName']


random_signatures.set_index('#Node').to_csv('Random_Signatures_25.tsv', sep = '\t')


#####

pathways = pd.read_csv('Pathways.tsv', sep='\t', header=None)

pathway_name = 'GLYCOLYSIS_/_GLUCONEOGENESIS_-_HOMO_SAPIENS_(HUMAN)'

one_pathway = pathways.where(pathways.iloc[:,1] == pathway_name).dropna()

one_pathway.to_csv('Pathway_One.tsv', sep='\t', header=None, index=False)

