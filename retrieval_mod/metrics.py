#-------------------------------------------------------------------------------
# Name:     metrics.py
#
# Purpose:  Compute statistical metrics for estimation purpose.
#           Formulas, i.e. nrmse, nbias were obtained from the paper below.
#
#           Paper name:
#           Zinevich et al. "Estimation of rainfall fields using commercial
#           microwave communication networks of variable density", 2008.
#-------------------------------------------------------------------------------


from numpy import sqrt, mean, array, std, exp, isscalar, nan
from scipy.stats import pearsonr


def average(val):
    '''
    Returns Mean value of a given array.

    Inputs:
        val  - 1D array. Quantity value
    Outputs:
        mean - scalar.
    '''
    return 1.0*sum(val) / len(val)


def stdev(val, BIAS=False):
    '''
    Returns Standard Deviation(STD) of a given array.

    Inputs:
        val   - 1D array. Quantity value
    Optional:
        BIAS  - boolean. If False, STD is computed without bias,
                         otherwise, with bias.
    Outputs:
        stdev - scalar.
    '''
    if isscalar(val) or len(val)==1:
        out = 0.0
    else:
        summa = 1.0*sum( abs(val - average( val ))**2 )
        if BIAS:
            out = sqrt( summa  / len(val) )
        else:
            out = sqrt( summa / (len(val)-1) )
    return out


def bias(obs, est):
    '''
    Returns Bias.
    Inputs:
        obs  - 1D array. Observed rain
        est  - 1D array. Estimated rain
    Outputs:
        bias - scalar.
    '''
    return 1.0*average(est - obs)


def pbias(obs, est):
    '''
    Returns Percent Bias between obs and est
    PBIAS = 100 * [ sum( sim - obs ) / sum( obs ) ]
    Percent bias (PBIAS) measures the average tendency of the simulated values
    to be larger or smaller than their observed ones.
    The optimal value of PBIAS is 0.0, with low-magnitude values indicating
    accurate model simulation. Positive values indicate overestimation bias,
    whereas negative values indicate model underestimation bias
    Formula:
        PBIAS = 100 * [ sum( sim - obs ) / sum( obs ) ]
    '''
    return 100 * ( 1.0 * sum(est - obs) / sum(obs) )


def mae(obs, est):
    '''
    Mean Absolute Error (MAE).
    '''
    return average( abs(est - obs) )


#TODO how to handle with "0" value in the denominator
def mape(obs, est):
    '''
    Mean Absolute Percentage Error [%].
    '''
    obs1 = obs[obs<>0]
    est1 = est[obs<>0]
    return 100 * average( abs( 1.0*(est1 - obs1) / obs1 ) )


def mse(obs, est, BIAS=False):
    '''
    Returns Mean Square Error(MSE).

    Inputs:
        obs  - 1D array. Observed rain.
        est  - 1D array. Estimated rain.
    Optional:
        BIAS - boolean. If "False", rmse is returned without bias,
                         If "True", rmse is returned with bias.
    Output:
        MSE  - scalar.
    '''
    if BIAS:
        out = (est - obs- bias(obs, est))**2
    else:
        out = (est - obs)**2
    return average(out)


#TODO how to handle with "0" value in the denominator
def mspe(obs, est, BIAS = False):
    '''
    Returns Mean Squared Percentage Error (MSPE) in [%].

    Inputs:
        obs  - 1D array. Observed rain.
        est  - 1D array. Estimated rain.
    Optional:
        BIAS - boolean. If "False", rmse is returned without bias,
                        If "True", rmse is returned with bias.
    Output:
        MSPE - scalar.
    '''
    obs1 = obs[obs<>0]
    est1 = est[obs<>0]

    if BIAS:
        out = (est1 - obs1 - bias(obs1, est1))**2
    else:
        out = (est1 - obs1)**2
    return 100 * average( 1.0 * out / obs1**2 )


def rmse(obs, est, BIAS=False):
    '''
    Returns Root Mean Square Error(rmse)

    Inputs:
        obs  - 1D array. Observed rain.
        est  - 1D array. Estimated rain.
    Optional:
        BIAS - boolean. If "False", rmse is returned without bias
                         If "True", rmse is returned with bias
    Outputs:
        rmse - scalar.
    '''
    return sqrt( mse( obs, est, BIAS ) )


def rmspe(obs, est, BIAS=False):
    '''
    Returns Root Mean Squared Percentage Error(RMSPE) in [%]

    Inputs:
        obs   - 1D array. Observed rain.
        est   - 1D array. Estimated rain.
    Optional:
        BIAS  - boolean. If "False", rmse is computed without bias
                        If "True", rmse is computed with bias
    Outputs:
        rmspe - scalar.
    '''
    return 100 * sqrt( mspe(obs, est, BIAS) / 100. )


def nbias(obs, est):
    '''
    Returns normalized bias (nbias).

    if nbias = 1, all rainfall is missed
    if nbias < 0,  algorithm overestimates the actual rainfall
    if nbias ~ 0, the correct reconstruction

    Inputs:
        obs   - (N, M) array. Observed rain. M - rain map number, N - its values
        est   - (N, M) array. Estimated rain. M - rain map number, N - its values
    Outputs:
        nbias - N dimensional vector.
    '''
    return average([  bias(i,j) / average(i) for i, j in zip(obs, est)])


def nrmse(obs, est, BIAS=False):
    '''
    Returns Normalized Root Mean Square Error(nrmse).

    Inputs:
        obs   - 2D array. Observed rain.
        est   - 2D array. Estimated rain.
    Optional:
        BIAS  - boolean. If "False", rmse is computed without bias
                        If "True", rmse is computed with bias
    Outputs:
        nrmse - N dimensional vector.
    '''
    return average([  sqrt( rmse(i, j, BIAS)**2 / average(( i - mean(j))**2))
                                                    for i, j in zip(obs, est)])

def ncorr(obs, est):
    '''
    The same document as in function "nbias"
    Inputs:
        obs   - (N, M) array. Observed rain. M - rain map number, N - its values
        est   - (N, M) array. Estimated rain. M - rain map number, N - its values
    Outputs:
        nnash - N dimensional vector.
    '''
    #avrg_est =
    nomi = [ ()]
    return [  corr(i,j) for i, j in zip(obs, est)]



def nash(obs, est):
    '''
    Nash-Sutcliffe (NS) criterion.

    Inputs:
        obs - 1D array. Observed rain.
        est - 1D array. Estimated rain.
    Outputs:
        NS  - scalar.
    '''
    denom = sum((obs - average(obs))**2)
    # To prevent zero division error
    if denom == 0:
        return nan
    return 1 - (1.0 * sum( (est - obs)**2) / denom)


def likelihood(obs, est, N=5.0):
    '''
    Returns Likelihood.

    Input:
        obs - 1D array. Observed rain.
        est - 1D array. Estimated rain.
    '''
    denom = sum((obs - average(obs))**2)
    # To prevent zero division error
    if denom == 0:
        return float("inf")
    return exp( - N * sum(( est - obs )**2) / denom)


def corr(obs, est):
    '''
    Pearsons Correlation coefficient.

    Inputs:
        obs  - 1D array. Observed rain.
        est  - 1D array. Estimated rain.
    Outputs:
        corr - scalar.
    '''
    return pearsonr(obs, est)[0]


#TODO This "rsquared" function should be corrected. It is wrong..!
def rsquared(obs, est):
    '''
    R squared effiency .
    Inputs:
        obs  - 1D array. Observed rain.
        est  - 1D array. Estimated rain.
    Outputs:
        R2 - scalar.
    Explanation:
        In statistics, the coefficient of determination, denoted R2 and
        pronounced R squared, indicates how well data points fit a statistical
        model- sometimes simply a line or curve. It is a statistic used
        in the context of statistical models whose main purpose is either
        the prediction of future outcomes or the testing of hypotheses,
        on the basis of other related information. It provides a measure of
        how well observed outcomes are replicated by the model, as the proportion
        of total variation of outcomes explained by the model
    '''
    o_mean = mean(obs)
    ss_tot = sum((i-o_mean)**2 for i in obs)
    ss_err = sum((i - j)**2 for i, j in zip(obs, est))
    r2 = 1 - (ss_err / ss_tot)
    return r2


def example_metric():
    import numpy as np
    frompath = 'C:/tomo_result/'
    obsname = 'imageInput_TS'
    estname = 'result_TS'
    maska_path = 'C:/ImageAreas_20IJ.txt'

    for each in ['lightrain', 'organized', 'unorganized','shower']:
        print each
        for i in range(1,21):

            oname = frompath + each+ '/' + obsname+str(i)+'_S20IJ.txt'
            ename = frompath + each+ '/' +  estname+str(i)+'_S20IJ.txt'

            obs = np.loadtxt(oname)
            est = np.loadtxt(ename)
            mask = np.loadtxt(maska_path)
            obs[mask == -999.] = -999.
            est[mask == -999.] = -999.

            obs = obs.reshape(400)
            est = est.reshape(400)

            obs = obs[obs<>-999]
            est = est[est<>-999]
            print round(rmse(obs, est),2)
    #        print ' nash: ', nash(obs, est)
    #        print rmse(obs, est)
     #       print ' bias: ', bias(obs, est)


    # Observed quantity
    obs = array([1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
                    16,17,18,30,31,32,33,34,35,36,37,38])
    # Estimated quantity
    est = array([0,1,2,3,4,5,6,7,8,8,9,10,11,12,13,14,
                    15,16,26,27,28,29,30,31,32,33,34])

    # Suppose, we have 3 samples with 9 parameters (3,9).
    obs1 = obs.reshape(3,9)
    est1 = est.reshape(3,9)
    #print 'NRMSE: ', nrmse(obs1, est1)
    #print 'NBIAS: ', nbias(obs1, est1)
    #for i, j in zip(obs1, est1):
    #    print ' corr: ', corr(i, j)
    #    print ' nash: ', nash(i, j)
    #    print ' rmse: ', rmse(i, j)
    #    print ' bias: ', bias(i, j)


if __name__ == '__main__':
    example_metric()