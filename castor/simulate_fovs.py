from functions_fov import *
from functions_simulations import *
import multiprocessing as mp
import itertools
import argparse


def parseArgs() :
    parser = argparse.ArgumentParser(
            description='simulate bin counts')
    parser.add_argument('-p','--threads', type=int, required=False, default=12, 
            help = "number of parallel processes (default:12)")
    parser.add_argument('-b','--binsize', type=int, required=False, default=21, 
            help = "number of pixels on one side of the bin square")
    parser.add_argument('-c','--countbins', type=str, required=True,
            help = "counts per bin csv file from bin_counts_per_gene")
    parser.add_argument('-o','--outfile', type=str, required=False,
            help = "output file, a csv")
    args = parser.parse_args()
    return args


def get_sim_inputs(countbins):
    '''
    read countbins csv file
    return a list of counts
    '''
    with open(countbins, 'r') as f:
        content = f.read().split('\n')
    counts = [int(x.split(',')[6]) for x in content[1:] if int(x.split(',')[6]) > 1 and len(x)>0]

    return counts
    
                     
def main():
    '''
    take bin_counts csv
    return csv of simulated bin groupsizes
    '''
    args = parseArgs()
    
    ## headers on output files
    with open(args.outfile, 'w') as f:
        header = 'simID,simcount,groupsize,count' + '\n'
        f.write(header)

    ## read the real bins file
    sim_inputs = get_sim_inputs(args.countbins)
        
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
