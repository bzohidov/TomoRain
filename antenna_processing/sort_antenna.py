
#-------------------------------------------------------------------------------
# Name:      sort_antenna.py
# Purpose: This module is for processing, sorting antenna data and saving it.
#          The main operations are "read", "sort" and "save".
#
# Description:
#     1. Function "read_anten":
#       This reads antenna file (only Excel file, *.xls)
#    2. Function "sort_anten":
#       Takes only chosen columns (details) of antenna information.
#       Then, it unifies connected antenna pairs(link) including each
#       antenna detail. Antenna details are given below in "Antenna details".
#    3. Function "save_anten":
#       Saves connected antenna (link) details
#
# Antenna Details:
#    Each antenna has the following details:
#    1. CAF         scalar;  Number which shows digital name of antenna.
#    2. Type:       string;  Transmitter(Trans) and Receiver (Rec).
#    3. Position:   list;    [x, y] Cartesian coordinate system.
#    4. Frequency:  scalar;  Operate frequency (unit is [MHz] for our file).
#    5. Azimuth:    scalar;  Azimuth angle.
#    6. Place:      string;  Place name where antenna is located.
#    7. Company:    string;  Company name which owns the antenna.
#
# Example for one sorted antenna pair (link) details:
#    Link = [ Caf, Frequency, Recx, Recy, Trx, Try, RecAzimuth, TrAzimuth
#             RecPlace, TrPlace, RecCompany, TrCompany ]
#
# NOTICE: In order to find connected Transmitter and Receiver pairs
#        sorting procedure are performed by comparing only
#        their "frequency" and "CAF". However, this sorting can include
#        considering "azimuth" angle by choosing "with_azim = True".
#-------------------------------------------------------------------------------


import numpy as np
from xlrd import open_workbook
from pyproj import Proj
from pandas import DataFrame, Series
import copy

# my lib
import antenna



TRANSMITTER = 'Trans'
RECEIVER = 'Rec'



def read_vars_from_cell(worksheet, row):
    '''
    Reads information from a given worksheet in Excel file.
    Columns: 'CAF': 4th; 'Num_ant': 5th; 'Azimuth': 6th;
    'Tr': 8th; 'Rec': 9th; 'Latitude': 12th; 'Longitude': 13th column.
    '''
    for i in [4, 5, 6, 10, 12, 13, 14]:
        if worksheet.cell_type(row, i) == 0:
            return
    caf = worksheet.cell_value(row, 4)
    num_anten = worksheet.cell_value(row, 5)
    azim = worksheet.cell_value(row, 6)
    place = worksheet.cell_value(row, 10)
    company = worksheet.cell_value(row, 14)
    if worksheet.cell_type(row, 8) == 2:
        freq = worksheet.cell_value(row, 8)
        anten_type = TRANSMITTER
    elif worksheet.cell_type(row, 9) == 2:
        freq = worksheet.cell_value(row, 9)
        anten_type = RECEIVER
    else:
        return
    lat = antenna.dms2dd(worksheet.cell_value(row, 12))
    lon = antenna.dms2dd(worksheet.cell_value(row, 13))
    return caf, num_anten, azim, anten_type, freq, lat, lon, place, company

def find_anten(antennas, freq, upper = 1000):
    '''
    Returns antenna information(antennas) where frequency (freq) value
    is not greater than 'frequency + upper' value.
    '''
    return [ant for ant in antennas
                  if ant.frequency >= freq and ant.frequency < freq + upper]

def read_anten(filename):
    '''
    Returns antenna information from a given Excel file (*.xls).
    '''
    workbook = open_workbook(filename)
    # Fetch Excel sheet name index
    worksheet = workbook.sheet_by_index(0)
    num_rows = worksheet.nrows - 1
    curr_row = 0
    ants = []
    while curr_row < num_rows:
            curr_row += 1
            row = worksheet.row(curr_row)
            data = read_vars_from_cell(worksheet, curr_row)
            if data is None:
                continue
            ants.append(antenna.Antenna(*data))
    return ants


def sort_anten(filename, freq=None, coor_system='cartesian', ant_type=TRANSMITTER):
    '''
    Returns connected microwave link coordinates in meter at a frequency
    (freq in MHz) using conversion in WGS84 projection system of 'UTM'.

    Inputs:
        filename - string. Antenna file name in Excel.
    Optional:
        freq     - scalar. Frequency of antenna [Hz]. If 'None' then all
                   frequency ranges, otherwise certain frequencies in the
                   file are returned.
        coor_system - string. Coordinate system. This can be 'cartesian'
                   or 'geodetic' coordinate system. By default, 'cartesian'.
        ant_type - string. Antenna type whether 'Trans' or 'Rec' which are
                   TRANSMITTER or RECEIVER, respectively.
                   antennas If ant_type = TRANSMITTER
    Outputs:
        link     - 2D array. Information about antenna pairs(link)

    **NOTICE**
    Frequency unit in the given file is in MHz format!

    '''
    ants = read_anten(filename)
    # get connected antenna pairs
    connected = antenna.get_connected(ants)
    if freq==None:
        pair_ants = connected
    else:
        freq1 = freq / 1e6  # Hz -> MHz
        pair_ants = find_anten(connected, freq1)
    trans = [ant for ant in pair_ants if ant.anten_type == ant_type]
    link = []
    if coor_system == 'cartesian':
        # Set WGS84 coordinate projection system to convert (lon,lat)-->(x,y)
        p = Proj(proj='utm', ellps='WGS84')
        for tr in trans:
            tr_x, tr_y = p(tr.longitude, tr.latitude)
            for rec in tr.pairs:
                rec_x, rec_y = p(rec.longitude, rec.latitude)
                link.append([ int(tr.caf), freq, rec_x, rec_y, tr_x, tr_y,
                                rec.azimuth, tr.azimuth, rec.place, tr.place,
                                rec.company, tr.company ])
    elif coor_system == 'geodetic':
        for tr in trans:
            for rec in tr.pairs:
               link.append([ int(tr.caf), freq,  rec.longitude, rec.latitude,
                             tr.longitude, tr.latitude, rec.azimuth, tr.azimuth,
                             rec.place, tr.place, rec.company, tr.company ])
    return link


#TODO create file_columns names for latitude and longitude in 'to_file'.
def save_anten( from_file, to_file, frequency,
                coor_system = 'cartesian', connected_to = TRANSMITTER):
    '''
    Saves information about link to the pathname at given frequency
    in *.csv format.
    Inputs:
        from_file   - string. Unsorted links file path which must exist
                              E.g.: "data/raw_microwave_links.xls"
        to_file     - string. Sorted links file destination path.
                              E.g.: "data/CARTESIAN_sorted_microwave_links.csv"
        frequency   - scalar or list.  Antenna frequency
                              E.g.: = [18*1e9, 23*1e9, 38*1e9] ==> (18, 23, 38) GHz
    Optional:
        coor_system - string. Coordinate system, [meter].
                              This can be 'cartesian' or 'geodetic' coordinate
                              system. By default, 'cartesian'.
        connected_to - string. It can be TRANSMITTER or RECEIVER to which
                               connected antennas are found.
    Output:
        info - string. Information about file that has been saved
                       in the given location.
    '''
    columns_CARTESIAN = [ "Caf", "Frequency", "Recx",       "Recy",
                          "Trx",      "Try",       "RecAzimuth", "TrAzimuth",
                          "RecPlace", "TrPlace",   "RecCompany", "TrCompany"]

    columns_GEODETIC = [ "Caf", "Frequency",    "Reclongitude",
                         "Reclatitude",  "Trlongitude",  "Trlatitude",
                         "RecAzimuth",   "TrAzimuth",    "RecPlace",
                         "TrPlace",      "RecCompany",   "TrCompany"]

    # Unsorted antenna data file *.xls
    if np.isscalar(frequency):
        frequency = [frequency]
    # Find antenna pairs and its corresponding details.
    link = sum([sort_anten( from_file, i, coor_system, connected_to)
                                                       for i in frequency], [])

    sorted_link = DataFrame(link, columns = columns_GEODETIC)
    # Save antenna details
    sorted_link.to_csv(to_file)
    info =  str(" File: ") + from_file + str(" sorted and saved as ") + to_file
    return info

# TODO. Example should be added.
def main():
    pass

if __name__ == '__main__':
    main()

