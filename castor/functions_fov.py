import numpy as np
import time

def read_idx(idxfile):
    ''' 
    take path of index file
    return list of lists, where each sublist has fov index info
    '''
    with open(idxfile, 'r') as f:
        content = f.read().split('\n')

    fovidx = []
    for i in content[1:]:
        if len(i)>0:
            fovinfo = i.split(',')
            fovidx.append([fovinfo[0], int(fovinfo[1]), int(fovinfo[2]), int(fovinfo[3])])

    return fovidx


def read_fov(txfile, fovidxinfo):
    '''
    take a txfile and fov index line
    return the fov info from the tx file
    '''
    with open(txfile, 'r') as f:
        f.seek(fovidxinfo[1], 0)
        fovcontent = f.read(fovidxinfo[2]).split('\n')
        f.close()

    return fovcontent


def get_fov_cellxgene(keeptx, geneset):
    '''
    take transcript info for an fov
    return a cell/gene matrix as list of lists
    '''
    countsdict = {}
    fov = keeptx[0].split(',')[0]
    for i in keeptx:
        cell = i.split(',')[1]
        gene = i.split(',')[8]
        if cell not in countsdict:
            countsdict[cell] = {x:0 for x in geneset}
        countsdict[cell][gene]+=1

    countsmat = []
    for cell in countsdict:
        cellmat = [fov, cell]
        for gene in geneset:
            cellmat.append(str(countsdict[cell][gene]))
        countsmat.append(cellmat)

    return(countsmat)


def get_fov_gene_counts(fovcontent):
    '''
    take raw fov tx info read straight from the tx file
    return list of [gene, count]
    '''
    countinfo = {}
    for i in fovcontent:
        if len(i)>0:
            if i.split(',')[8] not in countinfo:
                countinfo[i.split(',')[8]] = 1
            else:
                countinfo[i.split(',')[8]] += 1

    countlist = [[x, countinfo[x]] for x in countinfo]

    return countlist


def count_tx_in_bin(points, xmin, xmax, ymin, ymax):
    '''
    chatgptism
    take in coords as list of tuples
    return number that fall in the bin defined by the mins/maxes
    '''
    count = 0
    for x, y in points:
        if xmin <= x <= xmax and ymin <= y <= ymax:
            count += 1
    return count


def fov_bin_counts(fovcontent, binsize):
    '''
    divide fovs into square bins of binsize x binsize
    MAKE SURE THAT THE FOVCONTENT INCLUDES BOTH CLASSIFIED AND UNCLASSIFIED TX
    for each gene, report bin coords and counts
    takes a couple minutes for binsize=21
    '''
    fovcontent = list(filter(None, fovcontent))
    fov = fovcontent[0].split(',')[0]
    genes = [x.split(',')[8] for x in fovcontent]
    coords = [(int(x.split(',')[3]), int(x.split(',')[4])) for x in fovcontent]

    max_y = max(coords, key=lambda coord: coord[1])[1]
    max_x = max(coords, key=lambda coord: coord[0])[0]

    binrange_x = max_x//binsize ## floor division
    binrange_y = max_y//binsize ## floor division

    genebins = []
    for gene in set(genes):
        geneidx = [x for x in range(len(genes)) if genes[x] == gene]
        genecoords = [coords[x] for x in geneidx]
        points = np.array(genecoords)
        sorted_points = points[np.lexsort((points[:, 1], points[:, 0]))]
        ##loop thru bins each bin
        for xbin in range(binrange_x):
            xmax = (xbin+1)*binsize
            xmin = xbin*binsize + 1
            xmaxidx = np.searchsorted(sorted_points[:, 0], xmax, side='right')
            xminidx = np.searchsorted(sorted_points[:, 0], xmin, side='right')
            xbincoords = sorted_points[xminidx:xmaxidx]
            xbincoords_sorted = xbincoords[np.lexsort((xbincoords[:, 0], xbincoords[:, 1]))]
            for ybin in range(binrange_y):
                ymax = (ybin+1)*binsize
                ymin = ybin*binsize + 1
                ymaxidx = np.searchsorted(xbincoords_sorted[:, 1], ymax, side='right')
                yminidx = np.searchsorted(xbincoords_sorted[:, 1], ymin, side='right')
                count = ymaxidx - yminidx
                if count>0:
                    genebins.append([fov, gene, xmin, xmax, ymin, ymax, count])
    return genebins


def point_in_polygon(x, y, prepared_polygon):
    return prepared_polygon.contains(Point(x, y))


def count_tx_in_polygon(points, poly):
    '''
    see how many points are in a polygon
    '''
    preppoly = prep(poly)
    vectorized_function = np.vectorize(point_in_polygon)
    result = vectorized_function(points[:, 0], points[:, 1], preppoly)
    return np.sum(result)
