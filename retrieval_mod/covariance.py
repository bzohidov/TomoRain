#-------------------------------------------------------------------------------
# Name: covariance.py
# Purpose: This module provides apriori covariance matrix on parameter and data.
#
# 1) Uncertainties (variance) in data.
#    This covariance is assumed to be independent and the variance of
#    measurement is equally distributed.:
#                Cov_data = variance(d) * I,
#    where,
#         I - idendity matrix; d - data(attenuation) vector;
#
# 2) Variance in parameter represents confidence in parameter and as follows:
#         Cov_param(i,j) = ( sigmai*sigmaj ) * { exp(-0.5*( dmij/delta )^2)}
#    where,
#        sigmaij - scalar. Standard deviation of parameter at (i,j) location;
#        dij     - matrix. Distance which contains between ith and jth parameter.
#        decor   - scalar. Decorrelation distance [km]
#-------------------------------------------------------------------------------




from numpy import eye, around, exp




def covar_data(d_var):
    '''
    Returns data covariance matrix.
    Inputs:
        d_var - 1D array. Variance of data.
    Output:
        covar - 2D array. Data covariance matrix.
    '''
    covard = d_var * eye(len(d_var))
    return around(covard, 5)


def covar_param(p_var, decor, dij):
    '''
    Returns parameter covariance matrix.
    Inputs:
        p_var  - scalar. Parameter variance [-].
        decor  - scalar. Decorrelation distance [km].
        dij   - 2D array. Distance matrix [km].
    Output:
        covarp - 2D array. Parameter covariance.
    '''
    covarp = p_var * exp( -0.5 * (dij/ decor)**2 )
    return around(covarp, 5)
