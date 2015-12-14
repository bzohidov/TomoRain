#-------------------------------------------------------------------------------
# Name: scattering.py
# Purpose: This module provides electromagnetic scattering computation
#          depending on the drop diameter (water), temperature, frequency.
#
# Assumptions:
#    - scattering by electromagnetic wave is computed for water only;
#    - drop shape is assumed as spherical;
#    - mie scattering is used due to spherical shape of drop.
#
# The main utility here is to use 'mie_crossection' function which provides
# extinction cross section function for the calculation of rain attenuation.
#
# Functions were obtained fom scattering package by open source packages:
# https://github.com/dopplershift/Scattering.
#-------------------------------------------------------------------------------

import numpy as np
import scipy.special as ss
import scipy as sp
import scipy.constants as const


def water(temp):
    '''
    Calculate various parameters for the calculation of the dielectric constant
    of liquid water using the extended Debye formula. Temp is in Celsius.
    '''
    eps_s = 78.54*(1.0 - 4.579e-3 * (temp-25.0) + 1.19e-5 * (temp-25.0)**2 \
        - 2.8e-8 * (temp-25.0)**3)
    eps_inf = 5.27137 + 0.0216474*temp + 0.00131198*temp*temp
    alpha = -16.8129/(temp + 273) + 0.0609265
    lam_s = 0.00033836 * np.exp(2513.98/(temp + 273))
    sigma = 12.5664e8
    return eps_s, eps_inf, alpha, lam_s, sigma


def refractive_index(wavelength, temp):
    '''
    Calculates the complex refractive index using an expand Debye formula.
    The argument to the function gives another function which will return the
    necessary constants.  Temperature is in Celsius, Wavelength in m.

    Formula comes from Ray (1972).

    Example:
        >>> from scipy.constants import c
        >>> freq = 30*1e9  # 30*10^9 Hz
        >>> wavelength = c/freq
        >>> refractive_index(wavelength, 20)  # or arg. '20' can be exluded.
        (5.6204095522321467+2.787777223518217j)
    '''

    (eps_s, eps_inf, alpha, lam_s, sigma) = water(temp)
    wavelength /= const.centi
    lam_ratio = (lam_s / wavelength) ** (1 - alpha)
    sin_alpha = np.sin(np.pi * alpha / 2.0)
    denom = 1 + 2 * lam_ratio * sin_alpha + lam_ratio * lam_ratio
    eps_real = eps_inf + (eps_s - eps_inf) * ((1 + lam_ratio * sin_alpha)
        / denom)
    eps_imag = (eps_s - eps_inf) * lam_ratio * (np.cos(np.pi * alpha / 2.0)
        / denom) + sigma * wavelength / 18.8496e10

    return np.sqrt(eps_real + 1.0j * eps_imag)


def mie_abcd(m, x):
    '''Computes a matrix of Mie coefficients, a_n, b_n, c_n, d_n,
    of orders n=1 to nmax, complex refractive index m=m'+im",
    and size parameter x=k0*a, where k0= wave number
    in the ambient medium, a=sphere radius;
    p. 100, 477 in Bohren and Huffman (1983) BEWI:TDD122
    C. Matzler, June 2002

    Inputs:
        - m  - scalar, refractive-index;
        - x  - 1D array, size parameter;

    Outputs:
        - an - 1D array, Mie coefficient;
        - bn - 1D array, Mie coefficient.

    Example:
        >>> import scipy.special as ss
        >>> m = 5.62 + 2.78j     # refractive index of water at 20 [Celsius]
        >>> d = 0.003            # drop diameter - 3 [mm] --> 0.003 [m]
        >>> lam = 0.01           # wavelength -   10 [mm] --> 0.01  [m]
        >>> x = np.pi * d / lam
        >>> an, bn = mie_abcd(m, x)
        >>> an
        array([  3.28284411e-01 -3.36266733e-01j,
                 3.55370117e-03 -2.36426920e-02j,
                 3.90046294e-05 -5.15173202e-04j,
                 3.89223349e-07 -6.78235070e-06j,
                 2.94035836e-09 -5.85429913e-08j,
                 1.66462864e-11 -3.54797147e-10j,   7.17772878e-14 -1.58903782e-12j])
        >>> bn
        array([  7.69466132e-02 +1.28116779e-01j,
                 6.75408598e-03 +6.94984812e-03j,
                 1.99953419e-04 +8.05215919e-05j,
                 2.23969208e-06 +8.10161445e-08j,
                 1.39633882e-08 -2.58644268e-09j,
                 6.06975879e-11 -1.97461249e-11j,   2.01460986e-13 -8.42856046e-14j])
    '''

    nmax = np.round(2 + x + 4*x**( 1.0 / 3.0))
    mx = m * x

    # Get the spherical bessel functions of the first (j) and second (y) kind,
    # and their derivatives evaluated at x at order up to nmax
    j_x,jd_x,y_x,yd_x = ss.sph_jnyn(nmax, x)

    # The above function includes the 0 order bessel functions, which aren't used
    j_x = j_x[1:]
    jd_x = jd_x[1:]
    y_x = y_x[1:]
    yd_x = yd_x[1:]

    # Get the spherical Hankel function of the first type (and it's derivative)
    # from the combination of the bessel functions
    h1_x = j_x + 1.0j*y_x
    h1d_x = jd_x + 1.0j*yd_x

    # Get the spherical bessel function of the first kind and it's derivative
    # evaluated at mx
    j_mx,jd_mx = ss.sph_jn(nmax, mx)
    j_mx = j_mx[1:]
    jd_mx = jd_mx[1:]

    # Get primes (d/dx [x*f(x)]) using derivative product rule
    j_xp = j_x + x*jd_x
    j_mxp = j_mx + mx*jd_mx
    h1_xp = h1_x + x*h1d_x

    m2 = m * m
    an = (m2 * j_mx * j_xp - j_x * j_mxp)/(m2 * j_mx * h1_xp - h1_x * j_mxp)
    bn = (j_mx * j_xp - j_x * j_mxp)/(j_mx * h1_xp - h1_x * j_mxp)
    return an, bn


def mie_crossection(freq, d, temp):
    '''
    Computation of Mie Cross section for given frequency, diameter and temperature
    (which should have the same units), using complex Mie Coefficients
    'an' and 'bn' for n=1 to nmax, calculated using 'mie_abcd'.
    s. Bohren and Huffman (1983) BEWI:TDD122, p. 103,119-122,477.
    C. Matzler, May 2002.

    Inputs:
        - freq - scalar, microwave signal frequency, [Hz]
        - d    - 1D array, drop diameter,  [m]
        - temp - scalar, temperature,      [Celcius].

    Outputs:
        - sig_ext  - 1D array, scattering cross section;

    # Manually disabled outputs:
        - sig_sca  -  Extinction cross section;
        - sig_abs  -  Absolute cross section;

    # Optional outputs:
        - qsca  -  Scattering efficiency coefficient;
        - qext  -  Extinction efficiency coefficient;
        - qabs  -  Absolute efficiency coefficient.

    Example:
        >>> freq = 30*1e9     # 30*10^9  [Hz]
        >>> d = np.linspace(0, 0.007, 10) # 10 number of drop sizes [m]
        >>> temp = 20  # [Celsius]
        >>> mie_cros = mie_crossection(freq, d, temp)
        >>> print mie_cros
        [  0.00000000e+00   7.03810246e-08   1.81717219e-06   8.79068113e-06
           2.22328296e-05   3.47690997e-05   4.70248999e-05   6.39949298e-05
           8.42451110e-05   1.04395505e-04]
    '''

    lam = const.c/freq # in meter
    xs = np.pi * d / lam

    m = refractive_index( lam, temp)   # m - complex number
    if(m.imag < 0):
            m = np.conj(m)

    # Mie Cross sections
    sig_ext = np.zeros_like(xs)
    #sig_sca = np.zeros_like(xs)
    #sig_abs = np.zeros_like(xs)

##    # Mie Efficiencies
##    qsca = np.zeros_like(xs)
##    qext = np.zeros_like(xs)
##    qabs = np.zeros_like(xs)

    for i,x in enumerate(xs.flat):
        if float(x)==0.0:               # To avoid a singularity at x=0
##            qsca[i] = 0.0
##            qext[i] = 0.0
##            qabs[i] = 0.0
            sig_ext[i] = 0.0
            #sig_sca[i] = 0.0
            #sig_abs[i] = 0.0
        else:
            an,bn = mie_abcd(m,x)
            n = np.arange(1, an.size + 1)
            c = 2 * n + 1
            # Mie cross section computation
            lam2 = lam**2/(2*np.pi)
            sig_ext[i] = lam2 * (c * (an.real + bn.real)).sum()

            #sig_sca[i] = lam2 * (c * (np.abs(an)**2 + np.abs(bn)**2)).sum()
            #sig_abs[i] = sig_ext[i] - sig_sca[i]
##            # Mie efficiency computation
##            x2 = x * x
##            qsca[i] = (2.0 / x2) * (c * (np.abs(an)**2 + np.abs(bn)**2)).sum()
##            qext[i] = (2.0 / x2) * (c * (an.real + bn.real)).sum()
##            qabs[i] = qext[i] - qsca[i]
    return sig_ext #, sig_sca, sig_abs


def _test():
    import doctest, scattering
    doctest.testmod(scattering)

if __name__ == '__main__':
    _test()

