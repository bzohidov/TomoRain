#-------------------------------------------------------------------------------
# Name:        process_antenna.py
# Purpose:     Provides processing functions for link data.
#              "link" means Transmitter and Receiver pair.
#-------------------------------------------------------------------------------

import numpy as np
import pandas as pn
from copy import deepcopy



def get_link_at_freq(link, frequency):
    '''
    Returns links if their frequencies are at a given frequency
    Inputs:
        link     - DataFrame. Antenna information, at least, containing
                              'Frequency' column.
        frequency - list. Frequency of antenna, [Hz], E.g: [18*1e9] --> 18 GHz
    '''
    a = np.concatenate([ link[ link['Frequency'] == i] for i in frequency])
    new_link = pn.DataFrame(a, columns = link.columns)
    return new_link


def get_link_at_interval(link, interval):
    '''
    Returns links if their column 'Length' are within a given link *interval*.
    Inputs:
        link      - DataFrame. Antenna information.
        interval  - 1D list;   Min and Max value of link lengths, [km].
                               (i.e.: [min, max] ), e.g.: [1, 5].
    '''
    return link[np.logical_and( link['Length']>= interval[0],
                                link['Length']<= interval[1])]


def get_link_length(link):
    '''
    Returns lengths of links.
    Inputs:
        link  - DataFrame. Information containing Trx, Try, Recx, Recy
                           coordinates.
    '''
    return np.sqrt( (link['Trx'] - link['Recx'])**2 +
                    (link['Try'] - link['Recy'])**2)


def crop(link, box):
    '''
    Returns links inside the box.
    Inputs:
        links - Data Frame. Antenna data.
        box   - (2,2) list/tuple. Interest area. ((x0, x1), (y0, y1)), [km]
    Output:
        new_link - Data Frame. Antenna data
    '''
    ((x0, x1), (y0, y1)) = box
    inBox = lambda p: x0 <= p[0] <= x1 and y0 <= p[1] <= y1
    new_link = [ each for indx, each in link.iterrows()
                                    if inBox([ each['Trx'],  each['Try'] ]) &
                                       inBox([ each['Recx'], each['Recy'] ])]
    return pn.DataFrame(new_link, columns = link.columns)

#TODO Way of cutting the area for general case is needed.
def fix_position(link, border):
    '''
    Returns links with new coordinates in a new position.
    However, Length of links will not be changed.
    Inputs:
        link    - DataFrame. Antenna information.
        border  - tuple. Area border. ((x0, x1),(y0, y1)), [km].
    Outputs:
        links   - DataFrame. Links fixed in border with new positions.
    '''
    # Just to prevent messing up in system
    clink = deepcopy(link)
    x_coor = np.concatenate(( clink['Trx'], clink['Recx'] ))
    y_coor = np.concatenate(( clink['Try'], clink['Recy'] ))
    x_min, x_max = np.floor( min(x_coor)), np.ceil( max(x_coor))
    y_min, y_max = np.floor( min(y_coor)), np.ceil( max(y_coor))

    bx_diff = border[0][1] - border[0][0]
    by_diff = border[1][1] - border[1][0]

    # Cut the area with box min/max size
    diff_x = 0.5*abs(bx_diff - (x_max - x_min))
    diff_y = 0.5*abs(by_diff - (y_max - y_min))
    x_min = deepcopy(x_min) - diff_x
    y_min = deepcopy(y_min) - diff_y

    # Move link coordinate to new coordinates which begins in (0,0)
    clink[['Trx', 'Recx']] = clink[['Trx', 'Recx']] - x_min
    clink[['Try', 'Recy']] = clink[['Try', 'Recy']] - y_min

    return clink


def convert_coorunit(link, coor_unit = None):
    '''
    Converts link coordinates unit.
    Inputs:
        link       - DataFrame.
    Optional:
        coor_unit  - String.    Link coordinates unit in output.
                                Options: "km2meter" or "meter2km".
                                By default, None.
    Output:
        link       - DataFrame. Link with new coordinate unit
    '''
    clink = deepcopy(link)
    if coor_unit == 'meter2km':
            clink[['Recx', 'Recy',
                   'Trx',  'Try']] = 1e-3 * clink[[ 'Recx', 'Recy',
                                                    'Trx',  'Try']]
    elif coor_unit == 'km2meter':
            clink[[ 'Recx', 'Recy',
                    'Trx',  'Try']] = 1e3 * clink[[ 'Recx', 'Recy',
                                                    'Trx', 'Try']]
    return clink



def changefixcrop(link, border, box, coor_unit=None):
    #---------------------------------------------------------------------------
    #                          PROCESS MICROWAVE LINKS
    #---------------------------------------------------------------------------
    l2 = convert_coorunit(link, coor_unit)
    l3 = fix_position(l2, border)
    link = crop(l3, box)
    return link


