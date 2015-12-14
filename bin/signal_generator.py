#-------------------------------------------------------------------------------
# Name:        signal_generate.py
# Purpose:     Generate Electromagnetic Signal between Transmitter and Receiver
#
# Author:      Bahtiyor
#
# Created:     16-11-2014
# Copyright:   (c) Bahtiyor 2014
# Licence:     <MIT>
#-------------------------------------------------------------------------------

import pandas as pn
import numpy as np
import pylab as pl
from copy import copy

# my lib
import raindata
import simulation

# file path to rainfall maps
rain_file = str('data/rain_sorted_fullevent_025km')

# Power law coefficients.
a_mp = [0.05911087, 0.1080751, 0.37898495]
b_mp = [1.08693514, 1.05342886, 0.92876888]

#decor_dist = [0.5, 1,2,3,4,5,6,7,8,9,10,11,12,15,20,25,30]

def signal_data(savefile, alpha=0, delta=0, model_type = 'theoretical'):
    '''
    GENERATES RAIN ATTENUATION SIGNAL for a given noise percent
    *alpha* and quantization step *delta*
    Inputs:
        savefile   - string. File path to save the RSL data.
        alpha      - scalar. Noise percentage.
        delta      - scalar. Quantization step.
        model_type - string. Attenuation model type. 'theoretical' is set by default.
                             Avaiable models: 'theoretical', 'empirical'
    Output:
        None
    '''
    #*******************  DATA  *******************
    link = pn.read_csv(str('data/mobilenetwork.csv'))
    frequency = np.array([18, 23, 38])*1e9
    # Network area
    border = ((0.,40.),(0.,40.))
    ((x0, x1), (y0, y1)) = border
    # Monitor area
    box = border
    # Rainfall types
    rain_type = [ 'lightrain' ,'shower']


    #***************** START SIMULATION *****************

    # Set simulation framework for the links in "border" area
    experiment = simulation.MonitorArea( link, frequency, border)
    experiment.set_link()

    for i in rain_type:
          f_att_data = simulation.create_folder(savefile + '/RSL_025km/')
          maps, names =  raindata.read_rain(rain_file + '/'+i)
          att_measurement = []
          name = []
          for eachmap, k in zip(maps, names):
                eachmap = eachmap[4*y0:4*y1, 4*x0:4*x1]
                k = k.split('_')
                k.sort()
                event_time = k[1] + '_' +k[0]
                total_rsl = experiment.generate_rsl( eachmap, alpha, delta,
                                                     model = model_type)[0]
                att_measurement.append(total_rsl)
                name.append(event_time)
          data = np.array(att_measurement).T
          df05 = pn.DataFrame(data, columns = name)
          df05.to_csv( f_att_data +i+ '.csv')

    return None

# Example for generating signal
def example():
    alpha = 5  # noise [%]
    delta = 0  # quantization [dB]
    savepath = '../EXPERIMENT/'
    signal_data(savepath, alpha, delta)

if __name__ == '__main__':
    #example()
    pass

