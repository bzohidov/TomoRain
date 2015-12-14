#-------------------------------------------------------------------------------
# Name:        plot_statistics.py
# Purpose:     To plot statistical metrics (bias, correlation, nash, rmse)
#              for rainfall map (observed and estimated)
#-------------------------------------------------------------------------------

import pylab as pl
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import pandas as pn


#TODO documentation needed
def doplot(data, title, ylim=None, yaxis = 'Quantity', meanval=False,
                                                      color='b',
                                                      label='quantity',
                                                      showfig = False):

    pl.title(title, fontsize = 15, color='k')
    fig = pl.plot(range(len(data)), data, color, label=label, linewidth=1.5)
    pl.ylabel(yaxis)
    #pl.xticks( np.arange(0, len(data)+len(data)*0.1, len(data)*0.1 ) )
    x_mean = np.mean(data)
    if meanval==True:
        pl.axhline( x_mean, 0, len(data), color='r', linewidth=1.3, label='mean')
    if ylim<>None:
        pl.ylim(ylim)
    pl.xlabel('Number of maps')
    pl.xticks(np.arange(0,len(data),1))
    pl.yticks(np.arange(0,1.1,0.1))
    pl.grid(False)
    pl.gca().yaxis.grid(True)
    pl.legend(loc='upper left', numpoints = 1)

    if showfig:
        pl.show()


    return fig


#TODO documentation needed
def plot_stats(figtitle, bias, corr, nash, rmse, show_fig=True, save_fig=None):

    pl.clf()
    #fig = pl.figure(figtitle, figsize=(16,8), facecolor='#D8D8D8')
    pl.suptitle( figtitle, fontsize=20, color='r')
    pl.subplots_adjust(left=0.1, right=0.9, wspace=0.2, hspace=0.2)

    pl.subplot(221)
    ylim = (0,1)
    pl1 = doplot(nash, ylim=ylim,  title='$Nash-Sutcliffe$ $Efficiency $',
                       yaxis = 'NSE, [-]', label='NSE')
    pl.subplot(222)
    ylim =(0,1)
    pl2 = doplot(corr, ylim=ylim, title='$Correlation$ $coefficient$ ',
                       yaxis = 'correlation, [-]', label='correlation')
    pl.subplot(223)
    ylim = (round(min(bias),1), round(max(bias),1))
    pl3 = doplot(bias, ylim=ylim, title='$Bias $', yaxis='Bias, [mm/hour]',
                               label = 'bias')
    pl.xlabel('Number of maps')
    pl.subplot(224)
    ylim = (0, np.ceil(max(rmse)))
    pl4 = doplot(rmse, ylim=ylim,  title='$Root$ $Mean$ $Square$ $Error $',
                       yaxis = 'RMSE, [mm/hour]', label = 'rmse' )
    pl.xlabel('Number of maps')

    if show_fig == True:
        pl.show()
    if save_fig<>None:
        pl.savefig(save_fig + '.png')
    return None


def scatter_plot( obs, est, figname, showfig=True, savef=None,
                  title1='Linear scale,', title2= 'Log-log scale.',
                  xlabel1='Observed', xlabel2='Observed',
                  ylabel1='Estimated', ylabel2='Estimated'):
        '''
        Fit a line, y = mx + c, through some noisy data-points:
            x - observed
            y - estimated
            m - slope
            c - intercept
        '''
        x = obs
        y = est

        corr = pearsonr(x, y)[0]
        R2 = round(corr**2, 2)
        # Regression line
        A = np.vstack([ x, np.ones(len(x)) ]).T
        m, c = np.linalg.lstsq(A, y)[0]

        max_v = np.ceil( max(np.concatenate((x, y))) )
        t = np.arange(0, max_v)
        y_est = m*t + c
        y_obs = 1*t + 0
        tex_est = 'y = ' + str(round(m, 2)) + '*x + ' + str(round(c, 2))
        tex_obs = ' 1:1'

        pl.clf()
        fig = pl.figure(figname, figsize=(8, 6))
        pl.subplots_adjust(left=0.1, right=0.9, wspace=0.3, hspace=0.3)

        pl.subplot(111)
        pl.plot(x, y, 'b.')
        pl.plot(t, y_est, 'r-', label = tex_est)
        pl.plot(t, y_obs, 'k-', label = tex_obs)
        pl.xlim((0, max_v))
        pl.ylim((0, max_v))
        pl.title(title1 + ' $R^{2}$ = ' + str(R2))
        pl.xlabel(xlabel1 + ', [mm/hour]')
        pl.ylabel(ylabel1 + ', [mm/hour]')
        pl.legend()
        if showfig:
            # show the plot
            pl.show()
        if savef<> None:
            pl.savefig(savef + '.png')

def regression_line(x, y, showfig=True, savef = None):

    # linfit.py - example of confidence limit calculation for linear regression fitting.

    # References:
    # - Statistics in Geography by David Ebdon (ISBN: 978-0631136880)
    # - Reliability Engineering Resource Website:
    # - http://www.weibull.com/DOEWeb/confidence_intervals_in_simple_linear_regression.htm
    # - University of Glascow, Department of Statistics:
    # - http://www.stats.gla.ac.uk/steps/glossary/confidence_intervals.html#conflim

    # fit a curve to the data using a least squares 1st order polynomial fit
    z = np.polyfit(x,y,1)
    p = np.poly1d(z)
    fit = p(x)

    # get the coordinates for the fit curve
    c_y = [np.min(fit),np.max(fit)]
    c_x = [np.min(x),np.max(x)]

    # predict y values of origional data using the fit
    p_y = z[0] * x + z[1]

    # calculate the y-error (residuals)
    y_err = y -p_y

    # create series of new test x-values to predict for
    p_x = np.arange(np.min(x),np.max(x)+1,1)

    # now calculate confidence intervals for new test x-series
    mean_x = np.mean(x)         # mean of x
    n = len(x)              # number of samples in origional fit
    t = 2.31                # appropriate t value (where n=9, two tailed 95%)
    s_err = np.sum(np.power(y_err,2))   # sum of the squares of the residuals

    confs = t * np.sqrt((s_err/(n-2))*(1.0/n + (np.power((p_x-mean_x),2)/
                ((np.sum(np.power(x,2)))-n*(np.power(mean_x,2))))))

    # now predict y based on test x-values
    p_y = z[0]*p_x+z[0]

    # get lower and upper confidence limits based on predicted y and confidence intervals
    lower = p_y - abs(confs)
    upper = p_y + abs(confs)

    # maximum value of all points
    max_v = np.ceil( max(np.concatenate((x, y))) )

    # set-up the plot
    pl.axes().set_aspect('equal')
    pl.xlabel('X values')
    pl.ylabel('Y values')
    ##pl.title('Linear regression and confidence limits')
    pl.title('Linear regression')

    # plot sample data
    pl.plot(x,y,'b.',label='Sample observations')

    # plot line of best fit
    pl.plot(c_x,c_y,'r-',label='Regression line')

    # plot confidence limits
    pl.plot(p_x,lower,'b--',label='Lower confidence limit (95%)')
    pl.plot(p_x,upper,'b--',label='Upper confidence limit (95%)')

    # set coordinate limits
    pl.xlim(0,max_v)
    pl.ylim(0,max_v)

    # configure legend
    pl.legend(loc=0)
    leg = pl.gca().get_legend()
    ltext = leg.get_texts()
    pl.setp(ltext, fontsize=10)

    if showfig:
        # show the plot
        pl.show()
    if savef<> None:
        pl.savefig(savef)


# Example for plotting
def example():

    # Given any data
    data1 = np.random.randn(50)

    # Plot data
    doplot(data = data1, title= 'Random data plot')

    # Plot statistics in different metrics
    data2 = np.random.sample(50)
    plot_stats('DATA' ,data1, data2, data2, data2)

    # Regression line plot
    obs = np.array([4.0,2.5,3.2,5.8,7.4,4.4,8.3,8.5])
    est = np.array([2.1,4.0,1.5,6.3,5.0,5.8,8.1,7.1])
    regression_line(obs, est)

if __name__ == '__main__':
    #main()
    example()

