import argparse

def parseArgs() :
    parser = argparse.ArgumentParser(
            description='')
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-o','--outfile', type=str, required=False,
            help = "output file containing only unqiue transcripts")
    args = parser.parse_args()
    return args


def main():
    '''
    take txfile and idx file
    write duduped version of txfile, and deduped version of counts file
    '''
    args = parseArgs()
    header = 'fov,cell_ID,cell,x_local_px,y_local_px,x_global_px,y_global_px,z,target,CellComp'
    
    with open(args.txfile, 'r') as f:
        content = f.read().split('\n')
    content.remove(header)
    unique = list(set(content))
    unique.sort()
    
    ## headers on output file
    with open(args.outfile, 'w') as g:
        g.write(header)
        for i in unique:
            g.write(i+'\n')


if __name__ == "__main__":
   main()
