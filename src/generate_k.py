import pandas as pd
import numpy as np

## example path
path = "~/git/binox/example/"


## read original network file provided by BinoX example
network_file = 'Network.tsv'
network_proto = pd.read_csv(path + network_file, sep='\t')

## manipulate network
network_treshold = 0.99
network_cut = network_proto['#0:PFC'] >= network_treshold
network_one = network_proto.loc[network_cut]
network = network_one.iloc[:,[0,1]].reset_index(drop=True)
network.columns = ['gene1','gene2']
network = network.set_index('gene1')

def generate_k(pathway_genes):
    links = network.loc[pathway_genes].dropna()

    k_proto = links.gene2.value_counts().value_counts()

    k = k_proto.loc[list(range(max(k_proto.index)+1))]
    k[k.isna()] = 0
    k = np.array(k[1:], dtype='u4')
    
    return k