#-------------------------------------------------------------------------------
# Name:        rain_retrieval.py
# Purpose:
#
# Author:      Bahtiyor
#
# Created:     16-11-2014
# Copyright:   (c) Bahtiyor 2014
# Licence:     <your licence>
#-------------------------------------------------------------------------------

from numpy import loadtxt

# my lib
import raindata
import plots.plot_rainmap as plotrain
import simulation

def main():

    #---------------------------------------------------------------------------
    #                                INPUTS
    #---------------------------------------------------------------------------
    # Antenna frequency, [Hz]
    freq = [18*1e9, 23*1e9, 38*1e9]  # [Hz]
    # Area size, [km]
    monitor_area = ((0, 40)(0, 40))
    # resolution of the area, [km]
    area_resol = 0.25
    # retrieval resolution, [km]
    retri_resol = 0.5
    truemap = loadtxt('data/raindata/example_shower.txt') # Map size: 160x160
    MWlink =  read_csv('data/networkdata/example_network.csv')

    #---------------------------------------------------------------------------
    #                         LAUNCH SIMULATION MODEL
    #---------------------------------------------------------------------------
    nantes = simulation.MonitorArea(link, freq, monitor_area)
    nantes.set_rsl_discretization()
    nantes.set_retrieval_discretization()

    #---------------------------------------------------------------------------
    #                       GENERATE SIGNAL ATTENUATION
    #---------------------------------------------------------------------------
    signal = nantes.measure_att(truemap)
    signal_noisy = nantes.add_noise2rsl(truemap)

    #---------------------------------------------------------------------------
    #                            RAINFALL RETRIEVAL
    #---------------------------------------------------------------------------
    retrieved = nantes.reconstruct(signal_noisy)[-1]  # [-1] - retrieved map at "retri_resol"
    ground_truth = raindata.average_map(truemap, area_resol, retri_resol)

    #---------------------------------------------------------------------------
    #                    COMPARE "RETRIEVED" VS "GROUND-TRUTH"
    #---------------------------------------------------------------------------
    # a. Evaluate
    s = nantes.evaluate(ground_truth, retrieved, monitor_area, retri_resol)
    print 'Metrics: ', s
    # b. Plot
    plotting = plotrain.Rainplot('shower', '09-06-2012')
    plotting.compare(ground_truth, retrieved, retri_resol, retri_resol)


if __name__ == '__main__':
    main()
