from scipy.spatial import KDTree
from decimal import *

def find_nearest_like_neighbors(genecoords):
    '''
    Take list of lists representing coordinates [x,y]
    return list of distances from nearest neighbor
    '''
    tree = KDTree(genecoords)
    d, i = tree.query(x=genecoords, k=2)
    dists = [x[1] for x in d]

    return dists


def txnn(fovcontent):
    '''
    takes list of tx info
    returns list of tx info with [fov, x, y, target, txnn]
    '''
    fov = fovcontent[0].split(',')[0]
    
    coords = []
    genes = []
    for i in fovcontent:
        if len(i)>0:
            coords.append([int(i.split(',')[3]), int(i.split(',')[4]), int(i.split(',')[7])])
            genes.append(i.split(',')[8])

    geneset = list(set(genes))

    nninfo = []
    for gene in geneset:
        geneidx = [x for x in range(len(genes)) if genes[x] == gene]
        genecoords = [coords[x] for x in geneidx]
        dists = find_nearest_like_neighbors(genecoords)

        for i in range(len(genecoords)):
            nninfo.append([fov, str(genecoords[i][0]), str(genecoords[i][1]), gene, str(round(dists[i], 4))])

    return nninfo
    

