import random
from scipy.cluster.hierarchy import DisjointSet
from collections import Counter

def simulate_bin(binsize, count):
    '''
    given a binsize (px) and tx count
    simulate randomly distributed tx in a binsize x binsize x 8 px square
    return list of simulated points
    '''
    x = [random.randrange(binsize)+1 for i in range(count)]
    y = [random.randrange(binsize)+1 for i in range(count)]
    z = [random.randrange(8)+1 for i in range(count)]    

    points = list(zip(x, y, z))

    enough = False
    while not enough:
        if len(set(points)) == count:
            enough = True
        else:
            points.append((random.randrange(binsize)+1, random.randrange(binsize)+1, random.randrange(8)+1))

    return list(set(points))


def count_adjacent_tx(points):
    '''
    given a list of tx coordinates
    find which ones form adjacent groups
    return counts of groupsizes in the form [groupsize, count]
    '''
    tree = KDTree(points)
    pairs = tree.query_pairs(1.9)
    nearpairs = [list(pair) for pair in pairs]

    dset = DisjointSet(nearpairs[0])
    dset.merge(nearpairs[0][0], nearpairs[0][1])
    for i in nearpairs[1:]:
        dset.add(i[0])
        dset.add(i[1])
        dset.merge(i[0], i[1])
        
    neargroups = []
    for i in dset.subsets():
        neargroups.append(list(i))

    groupsizes = [len(x) for x in neargroups]
    sizefreqdict = Counter(groupsizes)
    sizefreq = [[x, sizefreqdict[x]] for x in sizefreqdict]
    
    return neargroups


