#-------------------------------------------------------------------------------
# Name:        representativity.py
# Purpose:     To provide pixel representativity map
#-------------------------------------------------------------------------------

import numpy as np
from itertools import chain

# my lib
from forward_process import discretize



def weight(lij):
    '''
    Returns weights for each crossed part of the link.

    Inputs:
        pij    - array of array; Crossed pixel index, e.g.: [[2,3], [4,5,2]...]
    Outputs:
        weight - array of arrays; Weights of crossed part of each link
                                  on each considered pixel.
    '''
    return np.array([np.array(each)/float(sum(each)) for each in lij])


def limit_weight_based(border, dx, dy, pij, lij, min_weight=None):
    '''
    Returns link index if its weight in considered pixel is equal or
    greater than *min_weight*.

    Inputs:
        border  -  (2,2) tuple. Area border: ((x0, x1),(y0, y1)), [km].
        dx      -  list.        Resolutions for retrieval in x axis, [km].
        dy      -  list.        Resolutions for retrieval in y axis, [km].
        Trx     -  1D array.    Transmitter coordinate in x axis, [km].
        Try     -  1D array.    Transmitter coordinate in y axis, [km].
        Recx    -  1D array.    Receiver coordinate in x axis, [km].
        Recy    -  1D array.    Receiver coordinate in y axis, [km].
    Optional:
        min_weight - scalar; Minimum value for links weights, ranges from 0 to 1.
                             By default, 'None' which means all weights of
                             the link crossed by the pixel are considered.
    Outputs:
        index  - array of arrays; link IDs where its weight is greater than
                                  given min_weight.
    '''
    wij = weight(lij)
    index = np.arange(len(pij))
    link_ind = []
    nomalum = np.unique( list(chain.from_iterable( pij )))
    for each in nomalum:
    	s = np.array([ [ i1, i3[i2 == each][0]]
    			         for i1, i2, i3 in zip(index, pij, wij) if each in i2 ])
    	# Get anten 'ID' if only its corresponding
    	# weight is greater than *min_weight*
    	k = np.where(s[:, 1] >= min_weight)[0]
    	if len(k) == 0:
    	    k = np.where(s[:, 1])[0]
    	link_ind.append( list(s[:, 0][k]) )

    return  link_ind


def maxweight_based( border, dx, dy, pij, lij ):
    '''
    Returns index of maximum weighted link for each considered pixel.

    Inputs:
        border  -  (2,2) tuple. Area border: ((x0, x1),(y0, y1)), [km].
        dx      -  list.        Resolutions for retrieval in x axis, [km].
        dy      -  list.        Resolutions for retrieval in y axis, [km].
        Trx     -  1D array.    Transmitter coordinate in x axis, [km].
        Try     -  1D array.    Transmitter coordinate in y axis, [km].
        Recx    -  1D array.    Receiver coordinate in x axis, [km].
        Recy    -  1D array.    Receiver coordinate in y axis, [km].
    Outputs:
        max_ind - 1D list; link indcies where its weight is maximum
    '''
    wij = weight(lij)
    index = np.arange(len(pij))
    max_ind = []
    nomalum = np.unique( list(chain.from_iterable( pij )))
    for each in nomalum:
        s = np.array([[ i4, i3[i1 == each][0], i2[i1 == each][0] ]
                        for i1, i2, i3, i4 in zip(pij, lij, wij, index)
                                                               if each in i1])
        # Find maximum weighted link index on each pixel
        max_wij = np.where(s[:, 1]<>1, max(s[:, 1]), max(s[:, 2]))[0]
        max_ind.append( s[:, 0][s[:, 1] == max_wij][0] )
    return max_ind

