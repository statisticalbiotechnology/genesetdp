import pandas as pd
import numpy as np
from genesetdp import find_maxscore

def genesetmc(k,q, n_runs):

    maxscore = find_maxscore(k,q)
    N_mc = np.zeros(maxscore+1)
    gene_links = np.array([ix for ix in range(len(k)) for _ in range(k[ix])])

    for i in range(n_runs):
        sample = sum(np.random.choice(gene_links, q, replace=False))
        N_mc[sample] += 1.0


    p_vals = (np.flipud(np.cumsum(np.flipud(N_mc))) - N_mc/2)/sum(N_mc)

    return p_vals


if __name__ == "__main__":
    import argparse
    import sys
    import pandas as pd
    import generate_k as genk

    parser = argparse.ArgumentParser(description='Calculates statistics for network enrichment of pathways.')
    parser.add_argument('networkfile', type=argparse.FileType('r'),
        help='Filename of Network file')
    parser.add_argument('pathwayfile', type=argparse.FileType('r'),
        help='Filename of Pathway file')
    parser.add_argument('querygenesfile', nargs='?', type=argparse.FileType('r'), default=sys.stdin,help='Filename of File with query genes. If omitted the gene names will be read from STDIN')
    parser.add_argument('-s','--singlepathway', metavar='PATHWAY',default='',help='The examined pathway. If omitted, all pathways in the pathway file will be tested.')
    parser.add_argument('-n','--networktreshold', default=0.7, metavar='float', type=float)

    parser.add_argument('-r','--montecarloruns', default=100000, metavar='N', type=int) ##

    args = parser.parse_args()

    query = pd.read_csv(args.querygenesfile, sep='\t', header=None)
    query.rename(columns = {0: 'gene'}, inplace =True)
    q = query.shape[0]

    pathways = pd.read_csv(args.pathwayfile, sep='\t', header=None)
    pathways.columns = ['gene','pathway']

    network_proto = pd.read_csv(args.networkfile, sep='\t')
    n_genes = len(np.unique(network_proto['#2:Gene1']))
    network_cut = network_proto['#0:PFC'] >= args.networktreshold
    network_one = network_proto.loc[network_cut]
    network = network_one.iloc[:,[0,1]].reset_index(drop=True)
    network.columns = ['gene1','gene2']
    network = network.set_index('gene1')

    ps = []

    if args.singlepathway:
        examinedpathways = [args.singlepathway]
    else:
        examinedpathways = np.unique(pathways['pathway'])

    for pathway_name in examinedpathways:
        pathway_genes = pathways.where(pathways.pathway == pathway_name).dropna().gene.values
        k = genk.generate_k(pathway_genes,network,n_genes)
        linked_genes = network.index[network.gene2.isin(pathway_genes)].tolist()
        s = sum(linked_genes.count(x) for x in query['gene'])

        # N = genesetdp(k,q)
        # p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)
        p_vals = genesetmc(k,q, args.montecarloruns)
        p = p_vals[s]
        ps.append((pathway_name,p))

    print("{}\t{}".format("Pathway","p-value"))
    for pathway_name,p in ps:
        print("{}\t{}".format(pathway_name,p))
