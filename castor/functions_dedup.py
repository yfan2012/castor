from sklearn.metrics.pairwise import paired_distances
from scipy.spatial import KDTree
from scipy.cluster.hierarchy import DisjointSet
from sklearn.metrics.pairwise import pairwise_distances
from castor.functions_fov import *
import multiprocessing as mp
import itertools
import time

def merge_overlapping_sublists(list_of_lists):
    '''
    from chatgpt
    '''
    merged = []
    while list_of_lists:
        # Start with the first sublist
        first, *rest = list_of_lists
        first = set(first)
        # Check if it overlaps with any other sublist
        length_before_merge = -1
        while length_before_merge != len(first):
            length_before_merge = len(first)
            rest2 = []
            for sublist in rest:
                sublist = set(sublist)
                if first & sublist:
                    first |= sublist  # Merge the sublists
                else:
                    rest2.append(sublist)
            rest = rest2
        merged.append(list(first))
        list_of_lists = rest
    return merged


def merge_overlapping(nearpairs):
    '''    
    take lists of pairs of points close together
    return lists of grouped points (transitive closeness)
    use disjointset/union find for speed
    '''
    dset = DisjointSet(nearpairs[0])
    dset.merge(nearpairs[0][0], nearpairs[0][1])
    for i in nearpairs[1:]:
        dset.add(i[0])
        dset.add(i[1])
        dset.merge(i[0], i[1])
        
    neargroups = []
    for i in dset.subsets():
        neargroups.append(list(i))
        
    return neargroups


def find_pairs_within_distance(genecoords, maxdist):
    '''
    Take list of lists of coordinates and a maximum tolerated distance
    return pairs of points within the maxdist of each other
    each point is identified by its index in the original coordinate list
    '''
    tree = KDTree(genecoords)
    pairs = tree.query_pairs(maxdist)
    nearpairs = [list(pair) for pair in pairs]
    
    return nearpairs


def dedup(fovcontent, maxdist):
    '''
    for a given fov list of transcript info (csv strings)
    return a deduped list in the same format
    consider transcripts within maxdist of each other to be the same
    allow chaining (if a is prox to b, and b is prox to c, all are collapsed)
    '''
    coords = []
    genes = []
    zplanes = []
    for i in fovcontent:
        if len(i)>0:
            coords.append([int(i.split(',')[3]), int(i.split(',')[4])])
            genes.append(i.split(',')[8])
            zplanes.append(int(i.split(',')[7]))

    geneset = list(set(genes))

    keeptx = []
    rmtx = []
    for gene in geneset:
        geneidx = [x for x in range(len(genes)) if genes[x] == gene]
        genecoords = [coords[x] for x in geneidx]
        nearpairs = find_pairs_within_distance(genecoords, maxdist)
        if not nearpairs:
            for i in geneidx:
                keeptx.append(fovcontent[i])
            continue
        neargroups = merge_overlapping(nearpairs)

        ## check for zplane groups
        ## could probs implement this with union join as well. maybe later
        for idx, group in enumerate(neargroups):
            zs = [[zplanes[geneidx[x]]] for x in group]
            zdists = pairwise_distances(zs, metric='euclidean')
            znear = []
            for i in zdists:
                znearidx = i <= 1
                znewset = [k for k, x in enumerate(znearidx) if x]
                if znewset not in znear:
                    znear.append(znewset)
            zmulti = []
            for i in znear:
                if len(i)==1:
                    keeptx.append(fovcontent[geneidx[int(group[int(i[0])])]])
                else:
                    txsubset = [set(i) <= set(x) for x in znear]
                    ## not subset of another clump
                    if sum(txsubset)==1:
                        zmulti.append(i)
            zgroups = merge_overlapping_sublists(zmulti)
            for zidx, i in enumerate(zgroups):
                keeptx.append(fovcontent[geneidx[int(group[int(i[0])])]]) ## keep one representative
                rmtx+=[fovcontent[geneidx[int(group[int(x)])]]+',g'+str(idx)+'_z'+str(zidx) for x in i] ## record group
                    
        ## deal with singletons
        notsingle = list(itertools.chain.from_iterable(neargroups))
        single = list(set(range(len(geneidx))) - set(notsingle))
        for i in single:
            keeptx.append(fovcontent[geneidx[i]])

    return keeptx, rmtx


def get_geneset(exprMatfile):
    '''
    extracts the header from an exprMatfile to get the full expected geneset
    returns list of gene names as strings
    '''
    with open(exprMatfile) as f:
        content = f.read().split('\n')
    geneset = content[0].split(',')[2:]
    
    return geneset


def process_fov(txfile, fovidxinfo, maxdist, exprMatfile, zerocountsfile, q, r, s, t):
    '''
    takes txfile and index of one fov
    calculates deduped tx and puts into q
    calculates counts mat and puts into r
    '''
    fovcontent = read_fov(txfile, fovidxinfo)
    keeptx, rmtx = dedup(fovcontent, maxdist)
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

    
## https://stackoverflow.com/questions/13446445/python-multiprocessing-safely-writing-to-a-file
def listener(outfile, q):
    '''
    listens to q
    writes contents to deduptxfile
    '''
    with open(outfile, 'a') as f:
        while True:
            m = q.get()
            if m == 'kill':
                break
            f.write(m)
            f.flush()

            
def dedup_voxel(fovcontent, maxdist):
    '''
    for a given fov list of transcript info (csv strings)
    return a deduped list in the same format
    consider transcripts within maxdist of each other to be the same
    allow chaining (if a is prox to b, and b is prox to c, all are collapsed)
    do this with voxels
    '''
    coords = []
    genes = []
    for i in fovcontent:
        if len(i)>0:
            coords.append([int(i.split(',')[3]), int(i.split(',')[4]), int(i.split(',')[7])])
            genes.append(i.split(',')[8])

    geneset = list(set(genes))

    keeptx = []
    rmtx = []
    for gene in geneset:
        geneidx = [x for x in range(len(genes)) if genes[x] == gene]
        genecoords = [coords[x] for x in geneidx]
        nearpairs = find_pairs_within_distance(genecoords, maxdist)
        if not nearpairs:
            for i in geneidx:
                keeptx.append(fovcontent[i])
            continue
        neargroups = merge_overlapping(nearpairs)

        for idx, i in enumerate(neargroups):
            keeptx.append(fovcontent[geneidx[i[0]]])
            rmtx+=[fovcontent[geneidx[x]]+',g'+str(idx) for x in i]

        ## deal with singletons
        notsingle = list(itertools.chain.from_iterable(neargroups))
        single = list(set(range(len(geneidx))) - set(notsingle))
        for i in single:
            keeptx.append(fovcontent[geneidx[i]])

    return keeptx, rmtx


