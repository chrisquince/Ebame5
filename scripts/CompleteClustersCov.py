import sys, getopt
import os
import numpy as np
import argparse
import math


def rchop(thestring, ending):
  if thestring.endswith(ending):
    return thestring[:-len(ending)]
  return thestring


def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("scg_file", help="scg frequencies tsv")

    parser.add_argument("cov_file", help="coverages csv")

    args = parser.parse_args()
    comp = {}
    for line in open(args.scg_file):
        line = line.rstrip()
    
        tokens = line.split("\t")

        count=0

        for tok in tokens[3:]:
            if tok == "1":
                count=count+1
        if count/36. > 0.75:
            comp[tokens[0]]=1
        else:
            comp[tokens[0]]=0
    #import ipdb; ipdb.set_trace()

    bFirst = True
    for line in open(args.cov_file):
        line = line.rstrip()

        if bFirst == False:
            tokens = line.split(",")
            id = tokens.pop(0)
            covArray = np.asarray([float(x) for x in tokens])
            count=0
            covSum = np.sum(covArray)
    #        print("Cluster" + id + "," + str(covSum))
            if covSum > 50.0 and comp[id] == 1:
                 print("Cluster" + id)
        else:
            bFirst = False

    #import ipdb; ipdb.set_trace()
if __name__ == "__main__":
    main(sys.argv[1:])
