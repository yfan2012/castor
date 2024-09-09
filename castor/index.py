import argparse

def parseArgs():
    parser=argparse.ArgumentParser(description='index a transcripts file by fov')
    parser.add_argument('-t', '--txfile', type=str, required=True, help='input decompressed tx file')
    parser.add_argument('-o', '--outidx', type=str, required=True, help='index file to write')
    args=parser.parse_args()
    return args


def make_idx(intxfile, outidx):
    '''
    index the tx file by fov
    returns a csv with fov, byte_offset, byte_len, numcells
    cell 0 denoted specially with fov_0, byte_offset, byte_len, 1
    inspiried by my megalodon mod basecalls script which snatched from nanocompore
    '''
    with open(outidx, 'w') as f:
        f.write(','.join(['fov', 'byte_offset', 'byte_len', 'numcells']) + '\n')
        txinfo = open(intxfile, 'r')
        fov = ''
        cell = ''
        byteoff = 0
        bytelen = 0
        numcells = 1
        for line in txinfo:
            ## not the first line
            if byteoff!=0:
                ## first fov
                if fov=='':
                    fov = line.split(',')[0]
                    cell = line.split(',')[1]
                    if cell!='0':
                        numcells += 1
                    bytelen += len(line)
                ## still on the same fov and cell
                elif fov==line.split(',')[0] and cell==line.split(',')[1]:
                    bytelen += len(line)
                ## same fov but different cell
                elif fov==line.split(',')[0] and cell!=line.split(',')[1]:
                    ## if the cell was a 0 cell, write out info and reset
                    if cell=='0':
                        f.write(','.join([fov+'_0', str(byteoff), str(bytelen), str(numcells)])+'\n')
                        byteoff += bytelen
                        fov = line.split(',')[0]
                        cell = line.split(',')[1]
                        bytelen = len(line)
                        numcells = 1
                    ## if the cell was not a 0 cell, reset the cell and continue adding to the byte length for this fov
                    else:
                        bytelen += len(line)
                        cell = line.split(',')[1]
                        numcells += 1
                ## new fov
                else:
                    f.write(','.join([fov, str(byteoff), str(bytelen), str(numcells)])+'\n')
                    byteoff += bytelen
                    fov = line.split(',')[0]
                    cell = line.split(',')[1]
                    bytelen = len(line)
                    numcells = 1
            ## first line
            else:
                ##add header
                byteoff+=len(line)

        ##write last fov
        f.write(','.join([fov, str(byteoff), str(bytelen), str(numcells)])+'\n')
        f.close()

def main(args):
    make_idx(args.txfile, args.outidx)

if __name__ == "__main__":
    args = parseArgs()
    main()



    
    
                        
