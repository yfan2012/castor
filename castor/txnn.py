from scipy.spatial import KDTree
from scipy.cluster.hierarchy import DisjointSet
from castor.functions_dedup import *
from castor.functions_txnn import *
import multiprocessing as mp
import itertools
import argparse
import time
start_time = time.time()

def parseArgs() :
    parser = argparse.ArgumentParser(
            description='report distance to the nearest like-neighbor for each transcript')
    parser.add_argument('-i','--idxfile', type=str, required=True,
            help = "index of tabular transcript csv file")
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-o','--outfile', type=str, required=True,
            help = "csv file of transcripts and the distance of their nearest neighbor")
    parser.add_argument('-p','--threads', type=int, required=False, default=12, 
            help = "number of parallel processes (default:12)")

    args = parser.parse_args()
    return args


def process_fov_txnn(txfile, fovidxpair, q):
    '''
    takes fov index
    returns each transcript with distance of nearest neighbor
    '''
    fovcontent = read_fov(txfile, fovidxpair[0]) + read_fov(txfile, fovidxpair[1]) 
    txnn_dists = txnn(fovcontent)
    write_txnn = '\n'.join([','.join(x) for x in txnn_dists])+'\n'
    q.put(write_txnn)


def main(args):
    '''
    take txfile and idx file
    return transcripts and distance to nearest like neighbor
    '''
    ## headers on output files
    with open(args.outfile, 'w') as f:
        header = ','.join(['fov', 'x', 'y', 'target', 'distance']) + '\n'
        f.write(header)


    ## set up multiprocessing
    manager = mp.Manager()
    q = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    watcher_q = pool.apply_async(listener, (args.outfile, q))
    jobs = []

    ## split up jobs
    fovidxpairs = []
    fovidx = read_idx(args.idxfile)
    for i in fovidx:
        if i[0].endswith('_0'):
            fovset = [x for x in fovidx if x[0] == i[0].replace('_0', '')]
            fovset.append(i)
            fovidxpairs.append(fovset)
    
    for fovidxpair in fovidxpairs:
        job = pool.apply_async(process_fov_txnn, (args.txfile, fovidxpair, q))
        jobs.append(job)

    ## run
    for job in jobs:
        job.get()

    ## tidy
    q.put('kill')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
    args = parseArgs()
    main(args)
