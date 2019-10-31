import sys, getopt
import os
import pandas as p
import numpy as np
import scipy.stats as ss
import scipy as sp
import scipy.misc as spm
import math
import argparse


from collections import defaultdict

def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("anvio_file", help="ANVIO format frequencies")

    parser.add_argument("hmm_file", help="ANVIO hmm calls")

    args = parser.parse_args()
    core_hmm = set()
    with open(args.hmm_file) as g:
        next(g)
        for line in g:
    
            line = line.rstrip()
            tokens = line.split('\t')
            core_hmm.add(tokens[1])

    anvio_file = args.anvio_file
#    import ipdb; ipdb.set_trace()
    sample_pos = defaultdict(lambda: defaultdict(dict))
    samples = set()
    genes = set()
    posns = defaultdict(set)
    bFirst = True
    mapPosns = {}
    
    vals = ('split_name','corresponding_gene_call','pos','A','C','G','T','sample_id')
    with open(args.anvio_file) as f:
        for line in f:
            line = line.rstrip()
            tokens = line.split('\t')

            if bFirst:
                header = tokens
                bFirst = False
            
                for val in vals:
                    mapPosns[val] = tokens.index(val)
            else:
                sample_id = tokens[mapPosns['sample_id']]
                gene_id = tokens[mapPosns['split_name']]
                pos = int(tokens[mapPosns['pos']])
                
                gene_call = tokens[mapPosns['corresponding_gene_call']]
                if gene_call in core_hmm:
                    
                    if sample_id not in samples:
                        samples.add(sample_id)
            
                    posns[gene_id].add(pos)
                    genes.add(gene_id)
                    freqs = (int(tokens[mapPosns['A']]),int(tokens[mapPosns['C']]), int(tokens[mapPosns['G']]),int(tokens[mapPosns['T']]))
                    sample_pos[gene_id][pos][sample_id] = freqs

    samples = list(samples)
    samples.sort()
    genes = list(genes)
    genes.sort()
    
    #sString = ",".join(samples)
    
    expanded = []
    for name in samples:
        expanded.append(name + "-A")
        expanded.append(name + "-C")
        expanded.append(name + "-G")
        expanded.append(name + "-T")

    eString = ",".join(expanded)
    print("Gene,Posn," + eString)
    for gene in genes:
        posns_gene = list(posns[gene])    
        posns_gene.sort()

        for posn in posns_gene:
            outString = gene + "," + repr(posn) 
            for sample in samples: 
                if sample in sample_pos[gene][posn]:
                    freqs = sample_pos[gene][posn][sample]
                else:
                    freqs = (0,0,0,0)
                freqs2 = [repr(f) for f in freqs]
                fstring = ",".join(freqs2)
                outString = outString + "," + fstring
            print(outString)
if __name__ == "__main__":
    main(sys.argv[1:])
