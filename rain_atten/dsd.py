
#-------------------------------------------------------------------------------
#  This module contains two Drop Size Distribution (DSD) models,
#  namely, Gamma and Marshal-Palmer DSD for different rainfall types.
#
#  The number and size of raindrops within a unit volume is described by
#  the number concentration, N(D) [number m^-3 mm^-1], also called DSD,
#  where D is the spherical equivalent diameter of each raindrop [mm].
#
#  General formula for computing DSD is as follows:
#
#                    N(D) =  N0 * D^mu * exp(-lam * D)       (1)
#  where,
#        N(D), [m^-3 mm^(-1-mu)] -  the number of drops per unit volume
#                                   per drop diameter interval (dD);
#        N0,   [m^-3 mm^-1]      -  scaling parameter;
#        D,    [mm]              -  drop diameter;
#        mu,   [unitless]        -  shape of DSD or 'mu' parameter;
#        lam,  [mm^-1]           -  slope parameter;
#
#  Formula (1) indicates Gamma DSD. If mu=0 then MP DSD can be obtained.
#
#  Parameter 'lam' depends on rain rate (R) which is as follows:
#
#                        lam = alpha*R^beta                  (2)
#   where,
#        R,     [mm/hour]   - rain rate;
#        alpha, [unitless]  - coefficient;
#        beta,  [unitless]  - coefficient;
#  For example, for Marshall Palmer model: alpha = 4.1 and beta = -0.21
#  Coefficients (alpha and beta) are given depending on rainfall type.
#
#  More info.:
#    1. Ondrej Fiser (2010). The Role of DSD and Radio Wave Scattering in
#       Rain Attenuation, Geoscience and Remote Sensing New Achievements,
#       Pasquale Imperatore and Daniele Riccio (Ed.), ISBN: 978-953-7619-97-8,
#       InTech, DOI: 10.5772/9110.
#    2. C. R. Williams and K. S. Gage, 2009. "Raindrop size distribution
#       variability estimated using ensemble statistics".
#    3. Slope parameters of DSD is based on IEEE 802.16cc-99/24 (Nov 1, 1999).
#-------------------------------------------------------------------------------


import numpy as np
import scipy as sp


def mp(D, R, rain_type):
    '''
    Returns Marshall Palmer DSD for a given D and R depending on rain type.

    Inputs:
           D - scalar or 1D array; drop diameters ranges, [m];
           R - scalar; rainfall rate ranges, [mm/hour];
           rain_type - string; it can be 'average' or 'shower' or
                                         'widespread' or 'drizzle'

    Outputs:
           mp_dsd - 1D array; [m^-3]

    Following DSDs are available for 'rain_type' of mp function:
    ----------------------------------------------------------
      -  average,                 [ Marshal Palmer, 1948 ];
      -  shower or thunderstrom,  [ - ];
      -  widespread,              [ Joss et al., 1969 ];
      -  drizzle,                 [ Joss et al., 1969 ];
    ----------------------------------------------------------

    Example:
      >>> D = np.linspace(0, 0.007,10)  # [m]
      >>> R = 10  # [mm/hour]
      >>> mp(D, R, rain_type = 'average')
      array([  8.00000000e+06,   1.11984311e+06,   1.56756073e+05,
               2.19427759e+04,   3.07155829e+03,   4.29957922e+02,
               6.01856768e+01,   8.42481441e+00,   1.17930879e+00,
               1.65080102e-01])
    '''
    # to change single arguments to 1D array, otherwise it remains unchanged
    try:
        D = D
    except:
        D = np.array([D])

    denom = R**(-0.21)

    if rain_type == 'average':
                   N0 = 1e6*8.0         # [m^-4]
                   lam = 1e3*4.1*denom  # [m^-1]
    elif rain_type == 'shower' or\
         rain_type == 'thunderstrom':
                   N0 = 1e6*1.4
                   lam = 1e2*30*denom
    elif rain_type == 'widespread':
                   N0 = 1e6*7.0
                   lam = 1e3*4.1*denom
    elif rain_type == 'drizzle':
                   N0 = 1e6*30.
                   lam = 1e2*57*denom
    else:
         raise IOError('rain_type: `average`, `shower or thunderstrom`, `widespread`, `drizzle`.')

    return N0 * np.exp(- D * lam)


def gamma(D, R, rain_type):
    '''
    Returns Gamma DSD for a given D and R depending on rain type.

    Inputs:
           D   -  scalar or 1D array,  drop diameters ranges, [meter];
           R   -  scalar,  rainfall rate ranges,  [mm/hour];
           rain_type -  string, rainfall type,  [-];

    Outputs:
           gamma_dsd - 1D array; [m^-3 m^(-1-mu)];

    Following DSDs are available for 'rain_type' of gamma function:
    --------------------------------------------------------
      -  zhang_model,             [ Zhang DSD model, 1999 ];
      -  convective,              [ Iguchi  T.,  1999 ];
      -  stratiform,              [ Iguchi  T.,  1999 ];
    --------------------------------------------------------
    Example:
      >>> D = np.linspace(0, 0.007,10)  # [m]
      >>> R = 10   # [mm/hour]
      >>> gamma(D, R, rain_type='zhang_model')
      array([  0.00000000e+00,   1.96854013e+06,   3.71204798e+05,
               2.95302199e+04,   1.64991744e+03,   7.59576597e+01,
               3.09381710e+00,   1.15801538e-01,   4.07445468e-03,
               1.36743402e-04])
    '''
    # to change single arguments to 1D array, otherwise it remains unchanged
    try:
        D = D
    except:
        D = np.array([D])
    # mu parameter
    mu = 3
    if rain_type =='zhang_model':
                 N0 = 0.125*(1.42*1e10)  # [cm^-4/cm^3]
                 N0 = 1e8*N0  # [m^-4/m^3] or m^7
                 lam = 1e2*0.5*130*R**(-0.13)   #[1/m]
    elif rain_type =='convective':
                 N0 =1e6* 6.29e5*R**(-0.416)
                 lam =1e2*8.35*R**(-0.185)
    elif rain_type =='stratiform':
                 N0 = 1e7*2.57e4*R**(0.012)
                 lam = 1e2*5.5*R**(-0.129)
    else:
        raise IOError('rain_type: `zhang_model`, `convective`, `stratiform`.')

    return (D)**mu * N0 * np.exp(-D * lam)


def _test():
    import doctest, dsd
    doctest.testmod(dsd)

if __name__ == '__main__':
    _test()
