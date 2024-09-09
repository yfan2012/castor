from sklearn.metrics.pairwise import paired_distances
from scipy.spatial import KDTree
from scipy.cluster.hierarchy import DisjointSet
from sklearn.metrics.pairwise import pairwise_distances
from functions_fov import *
from functions_dedup import *
import multiprocessing as mp
import itertools
import argparse
import time
start_time = time.time()

def parseArgs() :
    parser = argparse.ArgumentParser(
            description='extract tx info for a single fov or fov/gene pair')
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-i','--idxfile', type=str, required=True,
            help = "index of tabular transcript csv file")
    parser.add_argument('-f','--fov', type=str, required=True,
            help = "name of fov")
    parser.add_argument('-g','--gene', type=str, required=False,
            help = "gene to filter for")
    parser.add_argument('-o','--outfile', type=str, required=True,
            help = "deduplicated counts output file of cell zeros (unassigned transcripts)")
    args = parser.parse_args()
    return args


def main():
    '''
    take txfile and idx file
    write a tx file for a specific fov/gene combination
    '''
    args = parseArgs()
    
    ## headers the output file
    with open(args.txfile, 'r') as f, open(args.outfile, 'w') as g:
        header = f.readline()
        g.write(header)

    ## read the fov of interest
    fovidx = read_idx(args.idxfile)
    fovidxinfo = [x for x in fovidx if x[0] in [args.fov, args.fov+'_0']]

    fovcontent = []
    for i in fovidxinfo:
        fovcontent += filter(None, read_fov(args.txfile, i))
        
    if args.gene:
        towrite = [x for x in fovcontent if x.split(',')[8]==args.gene]
    else:
        towrite = fovcontent
    
    with open(args.outfile, 'a') as f:
        f.write('\n'.join(towrite)+'\n')
    
    
if __name__ == "__main__":
   main()
