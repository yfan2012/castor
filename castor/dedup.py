from sklearn.metrics.pairwise import paired_distances
from scipy.spatial import KDTree
from scipy.cluster.hierarchy import DisjointSet
from sklearn.metrics.pairwise import pairwise_distances
from castor.functions_fov import *
from castor.functions_dedup import *
import multiprocessing as mp
import itertools
import argparse


def parseArgs() :
    parser = argparse.ArgumentParser(
            description='deduplicate Nanostring tabular transcript files and corresponding counts files')
    parser.add_argument('-i','--idxfile', type=str, required=True,
            help = "index of tabular transcript csv file")
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-e','--exprMatfile', type=str, required=True,
            help = "tabular counts csv file")
    parser.add_argument('-d','--deduptxfile', type=str, required=True,
            help = "deduplicated tabular transcript csv output file")
    parser.add_argument('-r','--removedtxfile', type=str, required=True,
            help = "tabular transcript csv output file containing removed transcript groups")
    parser.add_argument('-c','--countsfile', type=str, required=True,
            help = "deduplicated counts output  file")
    parser.add_argument('-z','--zerocountsfile', type=str, required=False,
            help = "deduplicated counts output file of cell zeros (unassigned transcripts)")
    parser.add_argument('-m','--maxdist', type=float, required=False, default=1.9, 
            help = "distance to consider transcripts linked. use 1.5 to consider adjacent/diagonally adjacent pixels (default: 1.9)")
    parser.add_argument('-p','--threads', type=int, required=False, default=12, 
            help = "number of parallel processes (default:12)")
    args = parser.parse_args()
    return args


def process_fov_voxel(txfile, fovidxinfo, maxdist, exprMatfile, zerocountsfile, q, r, s, t):
    '''
    takes txfile and index of one fov
    calculates deduped tx and puts into q
    calculates counts mat and puts into r
    '''
    fovcontent = read_fov(txfile, fovidxinfo)
    keeptx, rmtx = dedup_voxel(fovcontent, maxdist)
    write_keeptx = '\n'.join(keeptx)+'\n'
    q.put(write_keeptx)

    ##if zeroth cell is asked for, write it
    if fovidxinfo[0].endswith('_0'):
        if len(zerocountsfile)>0:
            geneset = get_geneset(exprMatfile)
            countsmat = get_fov_cellxgene(keeptx, geneset)
            write_countsmat = '\n'.join([','.join(x) for x in countsmat])+'\n'
            t.put(write_countsmat)
    else:
        geneset = get_geneset(exprMatfile)
        countsmat = get_fov_cellxgene(keeptx, geneset)
        write_countsmat = '\n'.join([','.join(x) for x in countsmat])+'\n'
        r.put(write_countsmat)

    ##write info on removed groups
    write_rmtx = '\n'.join(rmtx)+'\n'
    s.put(write_rmtx)
    

def main(args):
    '''
    take txfile and idx file
    write duduped version of txfile, and deduped version of counts file
    '''

    ## headers on output files
    with open(args.txfile, 'r') as f, open(args.deduptxfile, 'w') as g:
        header = f.readline()
        g.write(header)

    with open(args.txfile, 'r') as f, open(args.removedtxfile, 'w') as g:
        header = f.readline().rstrip()+',group'+'\n'
        g.write(header)
        
    with open(args.exprMatfile, 'r') as f, open(args.countsfile, 'w') as g:
        header = f.readline()
        g.write(header)

    ## set up multiprocessing
    manager = mp.Manager()
    q = manager.Queue()
    r = manager.Queue()
    s = manager.Queue()
    pool = mp.Pool(processes=args.threads)
    watcher_q = pool.apply_async(listener, (args.deduptxfile, q))
    watcher_r = pool.apply_async(listener, (args.countsfile, r))
    watcher_s = pool.apply_async(listener, (args.removedtxfile, s))
    jobs = []

    if args.zerocountsfile:
        zerocountsfile=args.zerocountsfile
        with open(args.exprMatfile, 'r') as f, open(args.zerocountsfile, 'w') as g:
            header = f.readline()
            g.write(header)
        t = manager.Queue()
        watcher_t = pool.apply_async(listener, (zerocountsfile, t))
    else:
        zerocountsfile=''
        t = ''

    ## split up jobs
    fovidx = read_idx(args.idxfile)
    for i in fovidx:
        job = pool.apply_async(process_fov_voxel, (args.txfile, i, args.maxdist, args.exprMatfile, zerocountsfile, q, r, s, t))
        jobs.append(job)

    ## run
    for job in jobs:
        job.get()

    ## tidy
    q.put('kill')
    r.put('kill')
    s.put('kill')
    t.put('kill')
    pool.close()
    pool.join()

    
if __name__ == "__main__":
    args = parseArgs()
    main(args)
