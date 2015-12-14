#-------------------------------------------------------------------------------
# Name:        circle_detector1.py
# Purpose:     To provide detector that detects number of line in a given area
#              for each considered pixel.
#-------------------------------------------------------------------------------
import numpy as np

# my lib
from forward_process import discretize, flat_



def detect( A, B, C, r ):
    '''
    Returns True boolean if line segment "(A, B)" lies on a given
    circle "C" with radius "r", otherwise False.

    Inputs:
        A   -  1D array; Starting coordinate [x,y].
        B   -  1D array; Ending coordinate [x,y].
        C   -  list;     Circle coordinate, [x,y].
        r   -  scalar;   Circle radius.
    Output:
        boolean - True or False. True if given (A,B) lies on a circle.
    '''
    OA = complex( *A )
    OB = complex( *B )
    OC = complex( *C )
    # Now let's translate into a coordinate system where A is the origin
    AB = OB - OA
    AC = OC - OA
    # Before we go further let's cover one special case:  if either A or B is actually in
    # the circle,  then mark it as a detection
    BC = OC - OB
    if abs( BC ) < r or abs( AC ) < r: return True
    # Project C onto the line to find P, the point on the line that is closest to the circle centre
    AB_normalized = AB / abs( AB )
    AP_distance = AC.real * AB_normalized.real  +  AC.imag * AB_normalized.imag    # dot product (scalar result)
    AP = AP_distance * AB_normalized   # actual position of P relative to A (vector result)
    # If AB intersects the circle, and neither A nor B itself is in the circle,
    # then P, the point on the extended line that is closest to the circle centre, must be...
    # (1) ...within the segment AB:
    AP_proportion = AP_distance / abs( AB )   # scalar value: how far along AB is P?
    in_segment =  0 <= AP_proportion <= 1
    # ...and (2) within the circle:
    CP = AP - AC
    in_circle = abs( CP ) < r
    detected = in_circle and in_segment
    return  detected


#TODO documentation needed
def circle_coordinates(border, dx, dy, pix):
    '''
    Returns circle coordinates of an area (border) for only crossed
    pixels at a size of (dx, dy).
    Inputs:
        border   -  (2,2) tuple;   Area border: ((x_min, x_max),(y_min, y_max))
        dx       -  scalar;        Discretization in x axis in [km];
        dy       -  scalar;        Discretization in y axis in [km];
        pix      -  lists;         1 dimensional pixel indicies.
    Outputs:
        circle   -  (N,2) array;   Number of circle coordinates(x,y)
                                   for crossed pixels.

    **NOTICE**:  M - number of links; N - number of crossed pixels
    '''

    # Intersected lengths and pixels
    ((x0, x1), (y0, y1)) = border
    xticks = np.arange(x0, x1, dx) + dx/2
    yticks = np.arange(y0, y1, dy) + dy/2
    w = float(x1 - x0)/dx
    circle = [ [ [x,y] for i, x in enumerate(xticks)
                           if int(i + j*w) in pix] for j, y in enumerate(yticks) ]
    circle = sum(circle, [])
    return np.array(circle)


def constant_radius_based(border, dx, dy, Trx, Try, Recx, Recy, pix, radius):
    '''
    Returns link index for each crossed pixel which lies
    on a circle at a given radius.

    Inputs:
        border   -  (2,2) tuple;    Area border: ((x_min, x_max),(y_min, y_max))
        dx       -  scalar;         Discretization in x axis in [km];
        dy       -  scalar;         Discretization in y axis in [km];
        Trx      -  1D array(M,);   Transmitter coordinate in x axis.
        Try      -  1D array(M,);   Transmitter coordinate in y axis.
        Recx     -  1D array(M,);   Receiver coordinate in x axis.
        Recy     -  1D array(M,);   Receiver coordinate in y axis.
        pix      -  lists;          1 dimensional pixel indicies.
        radius   -  scalar;         Constant circle radius for all crossed pixels.
    Outputs:
        selected - list of lists.   Selected link indices in a given circle radius.

    **NOTICE**:  M - number of links; N - number of crossed pixels.
    '''
    circle = circle_coordinates(border, dx, dy, pix)
    Trn = np.vstack([Trx, Try]).T
    Recn = np.vstack([Recx, Recy]).T
    if len(Trx) == len(Try) == len(Recx) == len(Recy):
        ID = np.arange(len(Trx))
    selected = []
    for C in circle:
        count = [i for i, A, B in zip( ID, Trn, Recn) if detect(A, B, C, radius)]
        # if no link lies on a circle
        if len(count)==0:
            count.extend(ID)
        selected.append(count)
    return selected


#TODO documentation needed
def linknumber_based(border, dx, dy, Trx, Try, Recx, Recy, pix, link_num):
    '''
    Returns number of closest links indices for all crossed pixels
    Inputs:
        border   - (2,2) tuple; Area border: ((x_min, x_max),(y_min, y_max))
        dx       - scalar;      Discretization in x axis in [km];
        dy       - scalar;      Discretization in y axis in [km];
        Trx      - 1D array;    Transmitter coordinate in x axis.
        Try      - 1D array;    Transmitter coordinate in y axis.
        Recx     - 1D array;    Receiver coordinate in x axis.
        Recy     - 1D array;    Receiver coordinate in y axis.
        pix      -  lists;      1 dimensional pixel indicies.
        link_num - scalar;      Number of links that must be considered
                                for each crossed pixel area.
    Outputs:
        selected - list of lists. Selected link indices in a given circle radius.

    **NOTICE**:  M - number of links; N - number of crossed pixels.
    '''
    circle = circle_coordinates(border, dx, dy, pix)
    ((x0, x1), (y0, y1)) = border
    Trn = np.vstack([Trx, Try]).T
    Recn = np.vstack([Recx, Recy]).T
    if len(Trx) == len(Try) == len(Recx) == len(Recy):
        IDn = np.arange(len(Trx))
    selected = []

    # Min and Max size of radius
    rmin = 0
    rmax = np.ceil((x1 - x0)**2 + (y1 - y0)**2)
    for C in circle:
        new_ID = []
        ans_ID = []
        ID = copy(IDn)
        Tr = copy(Trn)
        Rec = copy(Recn)
        radius = copy(rmin)
        while len(new_ID) <> link_num and radius <= rmax:
            new_ID = [ i for i, A, B in zip( ID, Tr, Rec) if detect(A, B, C, radius)]
            if len(new_ID) <> 0:
                x = copy(IDn)
                x[new_ID] = -1
                ID = IDn[x<>-1]
                Tr = Trn[x<>-1]
                Rec = Recn[x<>-1]
                ans_ID += new_ID
                indexes = np.unique(ans_ID, return_index=True)[1]
                ans_ID = [ans_ID[index] for index in sorted(indexes)]
                if len(ans_ID) >= link_num:
                    ans_ID = ans_ID[:link_num]
                    break
            radius += 0.01
        selected.append(IDn[ans_ID].tolist())
    return selected
