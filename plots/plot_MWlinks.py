#-------------------------------------------------------------------------------
# Name:     plot_MWlinks.py
# Purpose:  This module provides plotting functions for microwave links.
#-------------------------------------------------------------------------------

import pylab as pl
from numpy import array, unique, arange
import pandas as pn
from PIL import Image

def plot_intersct(mwlinks, pnts, show='on'):
    '''
    Plots given intersected points showing each crossed part with the grid
    Inputs:
        mwlinks - (N,10) Data Frame. Microwave link containing 10 columns
                  [ID, freq, Recx, Recy, Trx, Try, Rec_azimuth, Tr_azimuth,
                                                                place, company]
        pnts   - (N,) array. Intersected points between beginning and ending
                             coordinates of microwave links.
        show   - string. Figure is displayed at 'on', not displayed at 'off'.
                         By default, it is set at 'on'.
    Outputs:
        No output
    '''
    f = pl.figure(num=1, figsize=(8, 6), dpi=80, facecolor='w', edgecolor='k')
    pl.subplot(111)
    for index, each in mwlinks.iterrrows():
        freq, rx, ry, tx, ty = each["Frequency"], each["Recx"], each["Recy"],\
                                                  each["Trx"],  each["Try"]
        rect = pl.Rectangle( (tx, ty), rx - tx, ry - ty,
                              facecolor="#60ff60", alpha=0.2 )
        pl.gca().add_patch(rect)
        pl.plot(pnts[i][1], pnts[i][2], pnts[i][3], pnts[i][4], 'bo')
    pl.xlabel()
    if show=='on':
        pl.show()



def plot_links( mwlinks, border, dx, dy, title=None, box=None,
                backgrnd_map = None, showfig='on' ,txt='on', grid=False):
    '''
    Plots microwave links.
    Inputs:
        mwlinks - (N,10) Data Frame. Microwave link containing 10 columns
                  [ID, freq, Recx, Recy, Trx, Try, Rec_azimuth, Tr_azimuth,
                   Rec_place, Rec_company, Tr_place, Tr_company]
        border  - (2,2) list/tuple. min and max coordinates of interest area.
                  ((x_min,x_max),(y_min,y_max)). By default, box = None
        dx      - scalar. Grid size in x axis.
        dy      - scalar. Grid size in y axis.
    Optional:
        backgrnd_map - string, file path in format `.*jpg`, `.*png`.
                               This will be plotted as background picture
                               behind microwave links. By default, it is None.
        show  -  string. Figure is displayed at 'on', not displayed at 'off'.
                         By default, it is set at 'on'.
    '''

    pl.figure(num=1)
    i = 0
    j = 0
    k = 0
    ((x_min, x_max), (y_min, y_max)) = border
    for index, each in mwlinks.iterrows():
        freq, rx, ry, tx, ty = each["Frequency"], each["Recx"], each["Recy"],\
                                                  each["Trx"],  each["Try"]
        if freq == 18*1e9:
            i +=1
            line1, = pl.plot([ tx, rx ], [ ty, ry ], 'b.-', linewidth = 2 )
        elif freq == 23*1e9:
            j +=1
            line2, = pl.plot([ tx, rx ], [ ty, ry ], 'r.-', linewidth = 2 )
        elif freq == 38*1e9:
            k +=1
            line3, = pl.plot([ tx, rx ], [ ty, ry ], 'g.-', linewidth = 2 )
    if box<>None:
        rect = pl.Rectangle( ( box[0][0], box[1][0]),
                               box[0][1] - box[0][0],
                               box[1][1] - box[1][0],
                               facecolor='yellow', alpha=0.4 ) #"#60ff60"
        pl.gca().add_patch(rect)

    pl.xticks(arange(x_min, x_max+dx, dx))
    pl.yticks(arange(y_min, y_max+dy, dy))
    if grid==True:
        pl.grid(True)
    pl.axis([x_min, x_max, y_min, y_max])
    pl.xlabel('x axis [East], km', fontsize=18)
    pl.ylabel('y axis [North], km', fontsize=18)
    used_freq = str(unique(mwlinks['Frequency'].values)/1e9)
    if title==None:
        title = 'Mobile phone networks at '+ used_freq +' GHz in Nantes'
        pl.title(title)
    # Put background image (city) behind microwave links
    if backgrnd_map <> None:
        im = array(Image.open(backgrnd_map))
        pl.imshow(im, extent= [x_min, x_max, y_min, y_max ])
    if showfig =='on':
        pl.show()