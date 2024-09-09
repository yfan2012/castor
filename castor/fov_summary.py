from dedup import read_idx, read_fov, get_geneset
import argparse

def parseArgs() :
    parser = argparse.ArgumentParser(
            description='generate some fov transcript summary info')
    parser.add_argument('-t','--txfile', type=str, required=True,
            help = "tabular transcript csv file")
    parser.add_argument('-i','--idxfile', type=str, required=True,
            help = "index of tabular transcript csv file")
    parser.add_argument('-e','--exprMatfile', type=str, required=True,
            help = "tabular counts csv file")
    parser.add_argument('-c','--countsfile', type=str, required=True,
            help = "deduplicated tabular transcript csv output file")
    args = parser.parse_args()
    return args


def get_fovxgene(txfile, idxfile, exprMatfile):
    '''
    generate an fovxgene file
    '''
    fovidx = read_idx(idxfile)
    geneset = get_geneset(exprMatfile)
    countsdict = {}
    for idx, i in enumerate(fovidx):
        print(idx)
        if len(i)>0:
            fovcount = []
            fovcontent = read_fov(txfile, i)
            fov = i[0]
            countsdict[fov] = {x:0 for x in geneset}
            for j in fovcontent:
                if len(j)>0:
                    gene = j.split(',')[8]
                    countsdict[fov][gene]+=1
    return countsdict
    
