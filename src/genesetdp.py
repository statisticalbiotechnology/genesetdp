#/usr/bin python
import numpy as np
from scipy.special import comb

def find_maxscore(k,q):
    max_s, ix, rem = 0, len(k)-1, q
    while rem > k[ix]:
        max_s += k[ix]*ix
        rem -= k[ix]
        ix -= 1
    max_s += rem*ix
    return max_s

def genesetdp(k,q):
    max_s = find_maxscore(k,q)
    #print(max_s,q)
    N = np.zeros((max_s+1,q+1))
    N[0,0] = 1

    for a in range(len(k)):
        for s in range(max_s,-1,-1):
            for c in range(q,-1,-1):
                for b in range(k[a],0,-1):
                    if c-b>=0 and s-a*b>=0:
                        N[s,c] += comb(k[a], b, exact=False) * N[s-a*b,c-b]
    return N[:,q]

def calculate_p(N):
    p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)
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

    args = parser.parse_args()

    query = pd.read_csv(args.querygenesfile, sep='\t', header=None)
    query.rename(columns ={0: ‘gene’}, inplace =True)
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

        N = genesetdp(k,q)
        p_vals = (np.flipud(np.cumsum(np.flipud(N))) - N/2)/sum(N)
        p = p_vals[s]
        ps.append((pathway_name,p))

    print("{}\t{}".format("Pathway","p-value"))
    for pathway_name,p in ps:
        print("{}\t{}".format(pathway_name,p))
