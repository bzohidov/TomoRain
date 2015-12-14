
#-------------------------------------------------------------------------------
# Name:        antenna.py
# Purpose:     Provide processing functions and class for antenna sorting out
#-------------------------------------------------------------------------------

from numpy import around
import unicodedata


TRANSMITTER = 'Trans'
RECEIVER = 'Rec'


def dms2dd(value_dms):
    '''
    Converts GEO coordinate degree minute second(dms)
    to degree decimal (dd).

    Input:
        value_dms - degree minute second GEO coordinate value.
    Output:
        Degree Decimal(dd).

    NOTICE:
        value_dms should be in 'Unicode' format.
        Here, data from Excel file are in Unicode format

    Unicode example:
        - Latitude: unicode('87\xb0 43\? 41\" N'); 'N' - North.
        - Longitude: unicode('1\xb0 35\? 35\" W'); 'W' - West.
    '''

    val = value_dms.split(u' ')
    dms_deg = float(val[0].split(u'\xb0')[0])
    dms_min = float(val[1].split(u'\'')[0])
    dms_sec = float(val[2].split(u'"')[0])
    val_dd = dms_deg + dms_min/60. + dms_sec/3600.
    # Change positive value to negative one if it is 'West'
    if value_dms[-1]=='W':
      val_dd = -1*val_dd

    return val_dd


def is_valid_azimuth(angle1, angle2):
    '''
    Returns True or False by comparing two angles whether both of them
    are in corresponding direction (connected with each other).

    Inputs:
        - angle1, angle2 : antenna azimuth angles in [degree].
    Outputs:
        - Boolean
    '''
    return ( angle2 <= 180 and angle1 == angle2 + 180 ) or\
           ( angle2 >  180 and angle1 == angle2 - 180 )

def get_connected(antennas, with_azim = False):
    '''
    Returns information of link if Transmitter and Receiver are connected.
    Information about link includes:
        caf, number, azimuth, anten_type, frequency, latitude, longitude.

    Inputs:
        antennas - 1D list. Information about link
    Optional:
        with_azim - Boolean. If 'True' then link is checked for connectivity
                             considering azimuth angle of Antenna pairs,
                             Otherwise, not.
    Output:
        connected - 1D list. Information about connected link

    '''

    connected = [] # set()
    transmitters = [ant for ant in antennas if ant.anten_type == TRANSMITTER]
    receivers = [ant for ant in antennas if ant.anten_type == RECEIVER]
    for trans in transmitters:
        for rec in receivers:
            if trans.caf == rec.caf and trans.frequency == rec.frequency:
                if not with_azim or is_valid_azimuth(trans.azimuth, rec.azimuth):
                    trans.add_pair(rec)
                    rec.add_pair(trans)
                    if trans not in connected:
                        connected.append(trans)
                    if rec not in connected:
                        connected.append(rec)
                    #connected |= set((trans, rec))
    return connected

class Antenna:
    def __init__(self, caf, number, azimuth, anten_type, frequency,
                       latitude, longitude, place, company):
        self.caf = caf
        self.number = number
        self.azimuth = azimuth
        self.anten_type = anten_type
        self.frequency = frequency
        self.latitude = latitude
        self.longitude = longitude
        self.place = unicodedata.normalize('NFKD', place).encode('ascii', 'ignore')
        self.company = unicodedata.normalize('NFKD', company).encode('ascii', 'ignore')
        self.pairs = []

    def add_pair(self, antenna):
        self.pairs.append(antenna)

    @property
    def address(self, address):
        return self.address

    @address.setter
    def address(self, address):
        self.address = address

    def __repr__(self):
        #TODO qolgan o'zgaruvchilar this
        return 'Caf: %s \t Number: %s \t Azimuth: %s \t Frequency: %s \
                ant_type: %s \t Longitude: %s \t Latitude: %s \n'\
                % (self.caf, self.number, self.azimuth, self.frequency, \
                   self.anten_type , self.longitude, self.latitude,
                   self.place, self.company)
