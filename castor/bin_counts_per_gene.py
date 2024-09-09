from scipy.spatial import KDTree
from scipy.cluster.hierarchy import DisjointSet
from functions_fov import *
from functions_dedup import *
import multiprocessing as mp
import itertools
import argparse
import time


def parseArgs() :
    parser = argparse.ArgumentParser(
            description='divide fovs into equally sized square bins and count the transcripts per bin')
    parser.add_argument('-p','--threads', type=int, required=False, default=12, 
            help = "number of parallel processes (default:12)")
    parser.add_argument('-b','--binsize', type=int, required=False, default=21, 
            help = "number of pixels on one side of the bin square")
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-i','--idxfile', type=str, required=True,
            help = "index of tabular transcript csv file")
    parser.add_argument('-o','--outfile', type=str, required=False,
            help = "output file, a csv")
    args = parser.parse_args()
    return args


def process_fov_bins(txfile, fovidxpair, binsize, q):
    '''
    takes fov
    returns bin counts
    '''
    fovcontent = read_fov(txfile, fovidxpair[0]) + read_fov(txfile, fovidxpair[1])
    fovcontent = list(filter(None, fovcontent))
    bincounts = fov_bin_counts(fovcontent, binsize)

    towrite = '\n'.join([','.join(map(str,x)) for x in bincounts]) + '\n'
    q.put(towrite)
    
                     
def main():
    '''
    take txfile and idx file
    write duduped version of txfile, and deduped version of counts file
    '''
    args = parseArgs()
    
    ## headers on output files
    with open(args.outfile, 'w') as f:
        header = 'fov,gene,xmin,xmax,ymin,ymax,count' + '\n'
        f.write(header)

    ## set up multiprocessing
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    watcher_q = pool.apply_async(listener, (args.outfile, q))
    jobs = []

    ## split up jobs
    fovidxlist = []
    fovidx = read_idx(args.idxfile)
    for i in fovidx:
        if i[0].endswith('_0'):
            fovset = [x for x in fovidx if x[0] == i[0].replace('_0', '')]
            fovset.append(i)
            fovidxlist.append(fovset)
    
    for i in fovidxlist:
        job = pool.apply_async(process_fov_bins, (args.txfile, i, args.binsize, q))
        jobs.append(job)

    ## run
    for job in jobs:
        job.get()

    ## tidy
    q.put('kill')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
   main()
