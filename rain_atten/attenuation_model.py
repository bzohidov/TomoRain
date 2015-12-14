#-------------------------------------------------------------------------------
#  Name: attenuation_model.py
#  Purpose: This module provides two rain attenuation models:
#           - Empirical model;
#           - Theoretical model.
#
# Rain Attenuation is the absorption of a microwave radio frequency (RF) signal
# by atmospheric rain, snow or ice, and losses which are especially prevalent
# at frequencies above 10 GHz.
#-------------------------------------------------------------------------------




import numpy as np
from scipy import linspace
from scipy.integrate import trapz

# from my lib.
import scattering as scat
import dsd
from powerlaw import get_coef_ab



def empirical(a,b, R, L):
    '''
    Returns attenuation along the link 'L' (between Transmitter and Receiver)
    using empirical model:
                          A = L * a * R^b
    Inputs:
          a,b  - scalar. Power law coefficients
          R    - scalar or 1D array. rain rate, [mm/hour];
          L    - scalar or 1D array. link length, [km];
    Output:
          A   -  scalar or 1D array, attenuation along the link, [dB];

    More info. on (a,b): https://www.itu.int/rec/R-REC-P.838-3-200503-I/en .

    # Example:
        #>>> # power law (a,b) at 18 GHz coef. by ITU-R, P.838-3-200503
        #>>> # a, b = 0.07393, 1.0404605978
        #>>> L = 1 # km
        #>>> R = np.array([1, 0, 15, 50 ])  # mm/hour
        #>>> freq = 18*1e9   # Hz
        #>>> empirical(freq, R, L)
        #array([ 0.05911087,  0.        ,  1.12202432,  4.15276687])
    '''
    A = np.dot(L, a * R**b)
    ##out =  (A / ( sum(L)*a) )**( 1./b )
    return A


def theoretical( freq, R, L, dsdtype = 'MP', D_min = 0, D_max = 7, D_num = 100,
                             temp = 20, rain_type = 'average'):
    '''
    Returns attenuation along the link 'L' (between Transmitter and Receiver)
    using theoretical model:

               A = L * log10(e) * integral(Q_ext(D) * N(D, R) dD)

    Inputs:
           freq  -  scalar. microwave signal frequency, [Hz];
           R     -  scalar. rain rate, [mm/hour];
           L     -  scalar. link length, [km];
    Optional:
           D_min -  min size of drop diameter, [mm];
                    By default, D_min = 0
           D_max -  max size of drop diameter, [mm];
                    By default, D_max = 7
           D_num -  number of drops in the range D_min and D_max;
                    By default, D_num = 100
           temp  -  scalar. temperature, [Celcius].
                    By default, temp = 20.
       dsd_type  -  string. type of drop size distribution: 'mp' or 'gamma';
                    By default, dsd_type = 'mp'.
       rain_type -  string. By default, rain_type = 'average'.
                    rain_type is chosen depending on dsd_type:
                    For 'mp' dsd_type: ('average', 'shower' or 'thunderstrom',
                                        'widespread', 'drizzle');
                    For 'gamma' dsd_type: ('zhang_model', 'convective',
                                           'stratiform').
    Output:
           A     -  scalar, attenuation along the link, [dB];

    # Example:
       >>> freq = 18*1e9 # frequency [Hz]
       >>> R = np.array([1, 0, 15, 50 ])
       >>> L = 1  # km
       >>> answer = [ theoretical( freq, r, L, dsdtype='gamma', temp=20, \
                                   rain_type='zhang_model' ) for r in R]
       >>> print 'GAMMA at 18 GHz: A = {0}'.format(answer)
       GAMMA at 18 GHz: A = [0.0464407, 0, 1.0105822, 3.7793992]
       >>> answer = [ theoretical( freq, r, L, dsdtype='MP', temp=20, \
                                   rain_type='average' ) for r in R]
       >>> print 'MP at 18 GHz: A = {0}'.format(answer)
       MP at 18 GHz: A = [0.0554683, 0, 1.2544763, 4.5017681]
    '''
    if R==0 or L==0 or freq==0:
        answer = 0
    else:
        # drop size ranges
        D_min = 1e-3 * D_min # [mm] --> [m]
        D_max = 1e-3 * D_max # [mm] --> [m]
        D = np.linspace(D_min, D_max, D_num)   # [m]
        # Extinction cross section function
        sig_ext_func = scat.mie_crossection(freq, D, temp)
        # Drop size distribution function
        if dsdtype == 'MP':
                    dsd_func = dsd.mp(D, R, rain_type='average')
        elif dsdtype == 'gamma':
                     dsd_func = dsd.gamma(D, R, rain_type='zhang_model')*0.5
        else:
            raise IOError('dsdtype should either `MP` or `gamma` , not '+str(dsdtype))
        # 10^6*4.4343.  10^6 is due to measurement unit adaptation
        const = 1e6 * 10 * np.log10(np.e)
        # Integrate attenuation
        att = L * const * trapz( sig_ext_func * dsd_func, x = D)
        answer = att * 1e-3   # [dB/meter]  ==> [dB/km]
        # Rounding the value
        answer = round(answer, 7)
    return answer


# Test example
def example():
    # Frequency
    freq = 18*1e9
    # Rain rates
    R = np.array([1, 0, 15, 50 ])
    # Link length [km]
    L = 1
    # Empirical model
    output = empirical(freq, R, L, mod_type = 'MP')
    print 'GAMMA at 18 GHz: A = {0}'.format(answer)
    # Theoretical model
    output = [ theoretical( freq, r, L, dsdtype='gamma', temp=20, \
                            rain_type='zhang_model' ) for r in R]
    print 'GAMMA at 18 GHz: A = {0}'.format(answer)


def _test():
    import doctest, attenuation_model
    doctest.testmod(attenuation_model)


if __name__ == '__main__':
    _test()
    #example()

