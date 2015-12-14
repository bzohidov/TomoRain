#-------------------------------------------------------------------------------
# Name:     quality_map_tools.py
# Purpose:  Provide quality map and link representativity in the area (map)
#-------------------------------------------------------------------------------

import numpy as np
from itertools import chain

def area_quality(lij, pixs, map_height, map_width):
    '''
    Generates quality map.
    Here, quality map is generated as follows:
        - First, Intersected part of the each link over each pixel is divided
                 by each link length, i.e.: pixel weights
        - Then, corresponding pixel weights in each pixel are summed up
    Inputs:
        lij        - array of arrays. Intersected parts of link lengths
        pixs       - array of arrays. Intersected pixel indices of links
        map_height - Height of quality map
        map_weight - Weight of quality map
    Outputs:
        qualmap - 2D array, (M,N). Quality map
    '''
    qualmap = np.zeros((map_height, map_width)).flatten()
    # representativity of pixel
    for l, p in zip(lij, pixs):
        qualmap[p] += np.array(l)/sum(l)
    return qualmap.reshape(map_height, map_width)


def numberlink_nesting(dx, dy, border, nesting_lij, nesting_pij):
    '''
    Returns pixel quality map in which each pixel value represents the number of
    link crossed this area (pixel)
    '''
    xsize = border[0][1] - border[0][0]
    ysize = border[1][1] - border[1][0]
    h = float(ysize) / dy[0]
    w = float(xsize) / dx[0]

    out1 = np.ones((h, w))*0.
    out2 = np.ones((h, w))*0.
    for i in range(len(dx)):
        hi = ysize / dy[i]
        wi = xsize / dx[i]
        mat = np.zeros((hi, wi))
        # 1. Q - matrix which shows number of links crossed each pixel
        # get only crossed pixels
        dd = np.array(list(chain.from_iterable(nesting_pij[i])))
        unknowns = np.unique( dd)
        b = np.unique(dd, return_index=True)
        Q = np.array([ len(dd[dd==j]) for j in b[0] ])   # number of crossed links in each pixel

        # 2. Q_weight - matrix which shows average weight of links crossed each pixel
        Q_weight = area_quality(nesting_lij[i], nesting_pij[i], hi, wi)
        Q_weight = Q_weight[Q_weight<>0.0].flatten()

        if i>0:
            ratio = dx[i-1] / dx[i]
            out1 = np.kron(out1, np.ones((ratio, ratio)))
            out2 = np.kron(out2, np.ones((ratio, ratio)))

        # out1 for Q
        out1 = out1.flatten()
        out1[[unknowns]] = Q
        out1 = out1.reshape(mat.shape)

        # out2 for Q_weight
        out2 = out2.flatten()
        out2[[unknowns]] = Q_weight
        out2 = out2.reshape(mat.shape)
    return out1, out2

