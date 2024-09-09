import argparse
from castor.index import main as index_main
from castor.txnn import main as txnn_main
from castor.dedup import main as dedup_main

def main():
    parser = argparse.ArgumentParser(description="castor CLI Tool")
    subparsers = parser.add_subparsers(dest='command')

    index_parser = subparsers.add_parser('index',
                                         help='index a transcripts file by fov')
    index_parser.add_argument('-t', '--txfile', type=str, required=True,
                              help='input decompressed tx file')
    index_parser.add_argument('-o', '--outidx', type=str, required=True,
                              help='index file to write')

    txnn_parser = subparsers.add_parser('txnn',
                                        help='report distance to the nearest like-neighbor for each transcript')
    txnn_parser.add_argument('-i', '--idxfile', type=str, required=True,
                             help = "index of tabular transcript csv file") 
    txnn_parser.add_argument('-t', '--txfile', type=str, required=True,
                             help = "tabular transcript csv file")
    txnn_parser.add_argument('-o','--outfile', type=str, required=True,
                             help = "csv file of transcripts and the distance of their nearest neighbor")
    txnn_parser.add_argument('-p', '--threads', type=int, required=False, default=12,
                             help = "number of parallel processes (default:12)")

    dedup_parser = subparsers.add_parser('dedup',
                                         help='deduplicate Nanostring tabular transcript files and corresponding counts files')
    dedup_parser.add_argument('-i', '--idxfile', type=str, required=True,
                              help = "index of tabular transcript csv file") 
    dedup_parser.add_argument('-t', '--txfile', type=str, required=True,
                              help = "tabular transcript csv file")
    dedup_parser.add_argument('-e', '--exprMatfile', type=str, required=True,
                              help = "tabular counts csv file")
    dedup_parser.add_argument('-c','--countsfile', type=str, required=True,
                              help = "deduplicated counts output  file")
    dedup_parser.add_argument('-d','--deduptxfile', type=str, required=True,
                              help = "deduplicated tabular transcript csv output file")
    dedup_parser.add_argument('-r','--removedtxfile', type=str, required=True,
                              help = "tabular transcript csv output file containing removed transcript groups")
    dedup_parser.add_argument('-z','--zerocountsfile', type=str, required=False,
                              help = "deduplicated counts output file of cell zeros (unassigned transcripts)")
    dedup_parser.add_argument('-m','--maxdist', type=float, required=False, default=1.9,
                              help = "distance to consider transcripts linked. (default:1.9)") 
    dedup_parser.add_argument('-p', '--threads', type=int, required=False, default=12,
                              help = "number of parallel processes (default:12)")

    args = parser.parse_args()

    if args.command == 'index':
        index_main(args)
    elif args.command == 'txnn':
        txnn_main(args)
    elif args.command == 'dedup':
        dedup_main(args)

if __name__ == "__main__":
    main()
