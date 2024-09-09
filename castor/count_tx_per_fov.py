from functions_fov import *
import argparse
import multiprocessing as mp

def parseArgs() :
    parser = argparse.ArgumentParser(
        description='generate some fov transcript summary info')
    parser.add_argument('-t','--txfile', type=str, required=True,
                        help = "tabular transcript csv file")
    parser.add_argument('-i','--idxfile', type=str, required=True,
                        help = "index of tabular transcript csv file")
    parser.add_argument('-o', '--outfile', type=str, required=True,
                        help = "output csv with lines: fov,gene,count")
    parser.add_argument('-p', '--processes', type=int, required=True,
                        help = "number of processes to use")
    args = parser.parse_args()
    return args


def count_genes_in_fov(txfile, fovidxinfo, q):
    '''
    read in fovcontent
    calculate gene count per fov
    write results to q
    '''
    fovcontent = read_fov(txfile, fovidxinfo)
    countslist = get_fov_gene_counts(fovcontent)
    write_counts = '\n'.join([','.join([fovidxinfo[0], x[0], str(x[1])]) for x in countslist])+'\n'
    q.put(write_counts)


def listener(outfile, q):
    '''
    listens to q
    writes contents to outfile
    '''
    with open(outfile, 'a') as f:
        while True:
            m = q.get()
            if m == 'kill':
                break
            f.write(m)
            f.flush()


def main():
    '''
    take the tx file and the idx file
    write a summary of gene counts per fov
    '''
    args = parseArgs()

    ##header on output
    with open(args.outfile, 'w') as f:
        f.write('fov,gene,count'+'\n')
    
    ## set up multiprocessing
    manager = mp.Manager()
    pool = mp.Pool(processes=args.processes)
    q = manager.Queue()
    watcher_q = pool.apply_async(listener, (args.outfile, q))
    jobs = []

    fovidx = read_idx(args.idxfile)
    for i in fovidx:
        job = pool.apply_async(count_genes_in_fov, (args.txfile, i, q))
        jobs.append(job)

    for job in jobs:
        job.get()

    q.put('kill')

    pool.close()
    pool.join()


if __name__ == "__main__":
   main()
    
