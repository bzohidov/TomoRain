#-------------------------------------------------------------------------------
# Name:        forward_process.py
# Purpose:     To provide set of functions which serve as discretization of an
#               area as well as used in retrieval process.
#-------------------------------------------------------------------------------


import numpy as np
from itertools import chain


def flat_(x):
    '''
    Flattens [[seq1],..[seq2],..[seqn]] into [seq1, seq2,.., seq3].

    Inputs:
        x - sequence (array or list); e.g.: seq = [[...],[...],[...],[...]]
    Outputs:
        1D array;
    '''
    return np.array( list(chain.from_iterable(x)) )


def pixel_distmat(rain, dx, dy, pij):
    '''
    Returns distance between pixel centres in a matrix form.

    Inputs:
        rain     - 2D array;      Rainfall map with no values.
        dx       - scalar;        Pixel resolution in x axis.
        dy       - scalar;        Pixel resolution in y axis.
        pij      - list of lists, (M,). Pixel indicies.
    Output:
        dist_mat - (N, N) array.  Distance matrix representing distance between
                                  pixel centers.
    '''
    height, width = np.shape(rain)
    pix = np.unique( flat_(pij) )
    # define all pixel centres as a unit of length in rain
    x = np.arange(0.5*dx, width*dx, dx)
    y = np.arange(0.5*dy, height*dy, dy)
    # convert 1D pixel coord. to 2D pixel coord.
    cols = np.mod(pix, width)
    rows = np.divide(pix, width)
    pix_coor = np.vstack((x[cols], y[rows])).T
    # compute distance matrix
    dist_mat =np.array( [ [np.sqrt((row[0]-col[0])**2 + (row[1]-col[1])**2)
                                              for col in pix_coor]
                                                         for row in pix_coor])
    return dist_mat


def get_intrsctd_pnts(Tr, Rec, dx, dy):
    '''
    Returns all intersected points of Transmitter and
    Receiver pair line at [dx,dy] discretization.
    Inputs:
        Tr     - (1,2) array;  Transmitter coordinate (x,y).
        Rec    - (1,2) array;  Receiver coordinate (x,y).
        dx     - scalar;       Discretization step in x axis.
        dy     - scalar;       Discretization step in y axis.
    Output:
        output - (N,2) array;  Intersected points of each link.
    '''
    # length between two points ( straight line) in x and y axis
    len_x = Rec[0] - Tr[0]
    len_y = Rec[1] - Tr[1]
    # Beginning and ending (min,max) points of each line
    lx_min = np.ceil(  min(Tr[0], Rec[0]) / dx ) * dx
    lx_max = np.floor( max(Tr[0], Rec[0]) / dx ) * dx
    ly_min = np.ceil(  min(Tr[1], Rec[1]) / dy ) * dy
    ly_max = np.floor( max(Tr[1], Rec[1]) / dy ) * dy
    AB=[]
    if len_y<>0. and len_x<>0.:
        #slope
        m = len_y/len_x
        # intercept
        b = Tr[1] - m*Tr[0]
        # calculate all intersected x coordinates
        x = np.arange(lx_min, lx_max + dx, dx)
        y = m * x + b
        AB = zip(x, y)
        # calculate all intersected y coordinates
        y = np.arange(ly_min, ly_max + dy, dy)
        x = (y - b) / m
        AB.extend(zip(x, y))
    elif len_x==0. and len_y<>0.:
        # calculate all intersected y coordinates
        y = np.arange(ly_min, ly_max + dy, dy)
        x = Tr[0]*np.ones(len(y))
        AB = zip(x,y)
    elif len_x<>0. and len_y==0.:
        # calculate all intersected x coordinates
        x = np.arange(lx_min, lx_max + dx, dx)
        y = Tr[1]*np.ones(len(x))
        AB = zip(x,y)
    else:
        raise IOError('Such a line does not exist with coordinates of '+ str((Tr,Rec)))
    AB.append((Tr[0],  Tr[1]))
    AB.append((Rec[0], Rec[1]))
    AB.sort()
    AB = np.asarray(AB)
    # If two points are in the same place (distance less than 1e-20),
    # one of them will be removed.
    idx = np.where(np.sum(np.diff(AB, axis=0)**2, -1) < 1e-20)[0]+1
    output = np.delete(AB, idx, axis=0)

    return output


def get_pixel(AB, w, dx, dy):
    '''
    Returns pixel indicies if AB (sets of lines segments) is crossed
    by discretization interval [dx, dy].
    Inputs:
        AB  - (N,2);   Array points AB(x,y) pairs.
        w   - scalar;  Area width.
        dx  - scalar;  Discretization step in x axis.
        dy  - scalar;  Discretization step in y axis.
    Outputs:
        pij - list of lists (M, ). Pixel indicies
    '''
    # AB - sets of lines segments; A- beginning, B-ending point
    pij = []

    pij = [ int( np.floor( min(AB[j - 1][0], AB[j][0]) /dx ) +
                 np.floor( min(AB[j - 1][1], AB[j][1]) /dy ) * w)
                                                  for j in range(1, len(AB)) ]
    return np.array(pij)


def get_sparse(pij, lij):
    '''
    Returns sparse matrix (M, N),
    where M - data number, N - unknown pixel number.

    Inputs:
        pij      - list of lists (M, ). Crossed pixel index,
                                        e.g.: [[2,3],[4], ...]
        lij      - list of lists (M, ). Crossed links part,
                                        e.g.: [[1.22, 3.433],[1.22], ...]
     Outputs:
        a_sparse - (M, N) array, Sparse matrix filled with link length.
    NOTICE:
        Each list len(pij) = len(lij), strictly.
    '''
    x = np.unique( flat_(pij) )
    a_sparse = np.zeros((len(pij), max(x)+1))
    for row, p, l in zip(a_sparse, pij, lij):
        row[p] = l
    return a_sparse[:, x]


def discretize(border, dx, dy, Trx, Try, Recx, Recy):
    '''
    Discretizes antenna pairs (Tr, Rec) in the area(border) with
    resolution dx, dy, then, calculates lengths of all crossed part of links.
    After, pij crossed by link part is found.

    Inputs:
        border - (2,2) tuple;  Area border: ((x_min, x_max),(y_min, y_max))
        dx     - scalar;       Discretization step in x axis.
        dy     - scalar;       Discretization step in y axis.
        Trx    - (M,1) array;  M number of Transmitter coordinates in x axis.
        Try    - (M,1) array;  M number of Transmitter coordinates in y axis.
        Recx   - (M,1) array;  M number of Receiver coordinates in x axis.
        Recy   - (M,1) array;  M number of Receiver coordinates in y axis.

    Outputs:
        pij    - list of lists; Crossed pixel index, e.g.: [[2,3],[4,5,2]...]
        lij    - list of lists; Crossed links part, e.g.: [[1.22, 3.433],..]
    '''
    ((x0, x1), (y0, y1)) = border
    Tr = np.vstack([Trx, Try]).T
    Rec = np.vstack([Recx, Recy]).T
    h, w = float(y1-y0)/dy, float(x1-x0)/dx
    pij = []
    lij = []
    for trans , recei in zip(Tr, Rec):
            # in case, study domain size is changed.
            trans = trans - np.array([ x0, y0 ])
            recei = recei - np.array([ x0, y0 ])
            AB = get_intrsctd_pnts(trans, recei, dx, dy)
            # Distance (lij) of each part of AB.
            each_len = np.sqrt(np.diff(AB[:,0])**2 + np.diff(AB[:,1])**2)
            each_pxl = get_pixel(AB, w, dx, dy)
            lij.append(each_len)
            pij.append(each_pxl)
    return np.array(lij), np.array(pij)



def _test():
    import doctest, forward_process
    doctest.testmod(forward_process)

if __name__ == '__main__':
    _test()