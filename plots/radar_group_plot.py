#-------------------------------------------------------------------------------
# Name:        radar_rainfall_plot.py
# Purpose:     Rainfall events plot and save them into a given folder
#-------------------------------------------------------------------------------
import os
import pylab as pl

# my lib
import simulation as simul
import raindata
import plots.plot_rainmap as plotrain


def plot_allevents(rain_type, from_folder, to_folder, showfig=True, savefig=False, figfmt='png'):
    '''
    Plots and saves all raster rain images.
    Inputs:
       rain_type   - String list. Rain types, e.g.: ['shower', ...]
       from_folder - String.  File folder path, e.g.: 'myrainfolders/'
       to_folder   - String.  File destination path in which files are saved.
       showfig     - Boolean. Status for displaying a plotted map. Figure is displayed
                              if True, otherwise not displayed.
       savefig     - Boolean. Option for saving plotted figures. Figure is saved
                              if True, otherwise not saved.
       figfmt      - String.  Figure format. By default "png".
                              Available formats: 'png', 'jpeg', 'tif', 'bmp'
    '''
    for each1 in rain_type:
    	dir1 = os.listdir(from_folder + each1)
    	for each2 in dir1:
            newfolder = simul.create_folder(to_folder + '/'+each1 +'/' +each2)
            maps, names = raindata.read_rain( from_folder+ '/'+each1 +'/' +each2)
            for eachmap, eachname in zip(maps, names):
                s = eachname.split('_')
                s.sort
                eventime = str(s[1] + '_'+s[2])
                Rplot = plotrain.Rainplot(each1, eventime, colormap='mine')
                title = 'Radar '+each1+' map'
                eachmap[eachmap==0] = -999
                im = Rplot.plotrain(eachmap, 0.25, 0.25, title=title, locator=5)
                if showfig:
                    pl.show()
                if savefig:
                    pl.savefig(newfolder+'//'+eachname+'.'+figfmt)
                pl.clf()