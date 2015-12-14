#-------------------------------------------------------------------------------
# Name: raindata.py
# Purpose: This module provides two functions (i) "read_rain" (ii) "average_map"
#
#-------------------------------------------------------------------------------




from numpy import loadtxt, array, shape, zeros
import copy
import os
import fnmatch




def read_rain(path_name, map_num=None):
    ''' Returns rainmap names and its values (2D array). All files in the
        path (path_name) are read, but empty files are not.
    Inputs:
        path_name - file path; e.g., '../my_folder/..'
        map_num  - scalar. Number of maps in the shown folder.
               By default, it is None which means all maps in a folder.
    Output:
        raster - (map_num, h, w) array. Here, h and w are height and width
                 of rain maps, respectively.
        nme    - list of strings. Each rain map file name (*.txt).
    Requirement:
    Files must be in the format of `*.txt`. '''
    rmap = []
    rname = []
    for filename in os.listdir( path_name ):
        # reads all *.txt files in the folder
        if fnmatch.fnmatch(filename, '*.txt'):
                    ID_file = path_name + '/' + filename
                    # Check if the file is empty or not
                    if os.path.getsize( ID_file ) == 0:
                          print ID_file + '  is empty file'
                    else:
                          rmap.append( loadtxt( ID_file ) )
                          rname.append(filename)
    return array(rmap)[:map_num], rname[:map_num]


def average_map(rain, dx_in, dx_out):

    ''' Returns average rain map at a resolution `dx_out` by averaging
        given `rain` values for dx_out/dx_in times.
    Inputs:
        rain   - 2D array, rain map
        dx_in  - scalar, rain map resolution, [km].
        dx_out - scalar, averaged rain map resolution, [km].
    Output:
        averaged rain - 2D array.
    Example:
        >>> rain = array([[ 0.51,  1.02,  0.19,  0.96],\
                          [ 0.87,  0.58,  0.55,  0.01],\
                          [ 0.07,  0.74,  1.4 ,  0.91],\
                          [ 1.68,  0.7 ,  2.03,  0.78]])
        >>> dx_in = 0.25; dx_out = 0.5   # input and output rain resolution [km]
        >>> new_rain = average_map(rain, 0.25, 0.5)
        >>> shape(new_rain)
        (2L, 2L)
        >>> new_rain
        array([[ 0.745 ,  0.4275],
               [ 0.7975,  1.28  ]])
    '''
    c_rain = copy.copy(rain)
    h, w = shape(c_rain)
    # ratio between dx_out and dx_in
    t = int(dx_out/dx_in)
    if t%2 <> 0 and t<>1:
        raise IOError('Ratio t= dx_out/dx_in = '+ str(dx_out)+'/' + str(dx_in)+
                                                    ' should be even number ')
    else:
        new_rain = [ [c_rain[ t*i:t*i + t, t*j:t*j + t].mean()
                              for j in range(w / t)] for i in range(h / t)]
    return array(new_rain)


def _test():
    import doctest, raindata
    doctest.testmod(raindata)


if __name__ == '__main__':
    _test()
