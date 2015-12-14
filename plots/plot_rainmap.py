#-------------------------------------------------------------------------------
# Name:        plot_rainmap.py
# Purpose:     Plot functions for rainfall map
#-------------------------------------------------------------------------------
import numpy as np
import pylab as pl
from numpy import shape, arange, loadtxt
import matplotlib as mpl


#-------------------------------------------------------------
#    No  | Color         | Code      | Rain rate (mm/hour)
#    -----------------------------------------------------
#    0   |  Light Grey   | -         | No rain area
#    1   |  White        | #FFFFFF   | 0 - 0.1
#    2   |  Sky-blue     | #87CEEB   | 0.1 - 0.4
#    3   |  Light Blue   | #ADD8E6   | 0.4 - 0.8
#    4   |  Blue         | #0000FF   | 0.8 - 1.5
#    5   |  Light Cyan   | #E0FFFF   | 1.5 - 2
#    6   |  Cyan         | #00FFFF   | 2 - 4
#    7   |  Dark Cyan    | #008B8B   | 4 - 7
#    8   |  Yellow       | #FFFF00   | 7 - 10
#    9   |  Yellow-orange| #FFC200   | 10 - 15
#    10  |  Orange       | #FFA500   | 15 - 25
#    11  |  Orange-red   | #FF4500   | 25 - 45
#    12  |  Red          | #FF0000   | 45 - 70
#    13  |  Dark Red     | #8B0000   | 70 - 120
#    14  |  Maroon       | #800000   | 120 - 133
#    15  |  Dark Brown   | #5C4033   | 133 - ~
#--------------------------------------------------------------


class Rainplot:
    '''
    This class provides two functions: `plotrain` and `compare`
    Inputs:
        rain_type  - String. Rain type, e.g.: str('light_rain')
        event_time - String. Rain event time, e.g.: str('26/07/2014, 13:30 UTC')
    Optional:
        colorbarvals- list.  Colorbar discrete values. Its length should be 15
        colormap   - String. Colormap type. This can be 'mine' or one of the
                             matplotlib colormap names, e.g.: 'jet'
                             By default,  colormap = 'mine' which is based on
                             specified color hex code.
        plot_type  - String. Type of plotting. This can be 'imshow' or 'contour'
    '''
    def __init__(self, rain_type, event_time, colorbarvals=None,
                                    colormap = 'jet', plot_type = 'imshow'):
        self.rain_type = rain_type      # string: 'lightrain'
        self.event_time = event_time    # string: '26072007 , 13:30 UTC'
        if colormap == 'mine':
            self.cmap = mpl.colors.ListedColormap([ '#F2F2F2', '#D8D8D8',
                                                    '#B2B2FF', '#6666FF',
                                                    '#0000FF', '#00FFFF',
                                                    '#009999', '#195E5E',
                                                    '#FFFF00', '#FFCC00',
                                                    '#FF6600', '#FF0000',
                                                    '#CC0000', '#CD0000',
                                                               '#800000'])
        else:
            self.cmap = pl.get_cmap(colormap)
        # choose plotting type
        self.plot_type = plot_type
        if colorbarvals==None:
            self.clevs = [0, 0.1, 0.4, 0.8, 1.5, 2, 4, 7, 10, 15, 25, 30, 45, 70, 120]
        else:
            if len(colorbarvals)<>15:
                raise IOError('len(colorbarvals) should be equal to 15')
            self.clevs = colorbarvals
        self.norm = mpl.colors.BoundaryNorm(self.clevs, self.cmap.N)

    def plotrain( self, rainmap,  dx, dy, title='Radar rainfall map', text=None,
                                  ax = None, locator = 5, xlabel = 'East, [km]',
                                  ylabel='North, [km]',
                                  cbar_label = 'Rainfall rate, [mm/hour]',
                                  shrink=1., pad=0.05):
        '''
        Plots rainmap at a given resolution
        Inputs:
            rainmap - 2D array; Rainfall map
            dx, dy  - scalar; pixel resolution in x,y axis, respectively.
            title   - string; Figure title, e.g: 'Comparison between...'
            text    - string; Image source, i.e.: this can be "Radar" or
                                                              "Microwave links"
        Optional:
            locator        - scalar.        Step in x and y axis.
            xlabel, ylabel - string.        Names of x, y axis.
            cbar_label     - string.        Colorbar name.
            shrink,pad     - scalar (each). Colorbar position feature.
        '''
        rainmap[rainmap==-999] = -999
        h, w = shape(rainmap)
        x_max, y_max = w*dx, h*dy
        extent = [0, x_max, 0, y_max]
        if self.plot_type == 'imshow':
            im = pl.imshow(  rainmap,
                             extent = extent,
                             cmap = self.cmap,
                             norm = self.norm,
                             vmin = 0.,
                             origin = 'upper',
                             interpolation = 'none')
        elif self.plot_type == 'contour':
            im = pl.contourf( rainmap,
                              levels = self.clevs,
                              cmap = self.cmap,
                              norm = self.norm,
                              origin = 'upper',
                              interpolation = 'none',
                              extend = "both")
        im.cmap.set_under('#FFFFFF')   # white color
        im.cmap.set_over('#110000')    # dark Maroon
        pl.title(title + ' $ dx = '+str(dx)+' km$', fontsize=14)
        if text<>None:
            pl.text(dx+0.5, dy+0.5, text, backgroundcolor='#D6D6D6', color='b',
                                                                 fontsize=20)
        pl.xticks(arange(0, x_max + locator, locator))
        pl.yticks(arange(0, y_max + locator, locator))
        pl.xlabel(xlabel, fontsize = 14)
        pl.ylabel(ylabel, fontsize = 14)
        cbar = pl.colorbar(im, ticks = self.clevs, cax = ax, extend='both',
                                       shrink=shrink, pad = pad)
        cbar.set_ticklabels(self.clevs)
        cbar.set_label(cbar_label)
        return im

    def compare(self, observed, estimated, dx, dy, text= None, showfig = False,
                                                               savefig = None):
        '''
        Plots two subplots
        Inputs:
            observed  - 2D array; Observed rainmap.
            estimated - 2D array;  Estimated rainmap.
        Optional:
            showfig   - boolean; if True, figure is displayed, otherwise not.
            savefig   - string; filename to save the figure.
        '''
        # We dont want to show un-estimated pixels
        observed[estimated==-999] = -999
        # Figure name
        figname = self.rain_type + '_' + self.event_time +'_'+ str(dx) + ' km'

        pl.clf()
        # Configure figure
        fig = pl.figure( figname, figsize=(14,6))
        fig.subplots_adjust(wspace=0.2, hspace=0.6)
        cbarshrink, cbarpad = 0.9, 0.06

        ax1 = fig.add_subplot(121)
        title1 = 'Observed  '
        text1 = 'Radar'
        im1 = self.plotrain(observed, dx, dy, title1, text1, shrink = cbarshrink,
                                                                pad = cbarpad)
        ax2 = fig.add_subplot(122)
        if text==None:
            title2 = 'Estimated '
        else:
            title2 = 'Estimated ' + text
        text2 = 'HF links'
        im2 = self.plotrain(estimated, dx, dy, title2, text2,shrink = cbarshrink,
                                                                pad = cbarpad)
        if showfig:
            pl.show()
        # Save plotted image
        if savefig <> None:
            fig.savefig(savefig + figname + str('.png'))
            print 'Saved: ' + savefig + figname + str('.png')
        return fig


# Example
def pl_example():
    R_obs = abs(np.random.randn(160,160))
    R_est = 10*abs(np.random.randn(160,160))
    rainplot = Rainplot('Shower', ' 06-11-2009, 18:25', colormap = 'jet')
    rainplot.plotrain(R_obs, dx=0.25, dy=0.25, title ='Shower', text='')
    pl.show()
    rainplot.compare(R_obs, R_est, 0.25, 0.25, showfig=True)





if __name__ == '__main__':
    pl_example()

