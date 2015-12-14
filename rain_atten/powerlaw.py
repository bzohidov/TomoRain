#-------------------------------------------------------------------------------
# Name: powerlaw.py
# Purpose: This is a set of power law coefficients(a,b) for the calculation of
#          empirical rain attenuation model A = a*R^b.
#-------------------------------------------------------------------------------




def get_coef_ab(freq, mod_type = 'MP'):
    '''
    Returns power law coefficients according to model type (mod_type)
    at a given frequency[Hz].
    Input:
        freq - frequency, [Hz].
    Optional:
        mod_type - model type for power law. This can be 'ITU_R2005',
                   'MP','GAMMA'. By default, mode_type = 'ITU_R2005'.
    Output:
        a,b - power law coefficients.
    '''
    if mod_type == 'ITU_R2005':
        a = {'18': 0.07393, '23': 0.1285, '38': 0.39225}
        b = {'18':1.0404605978, '23': 0.99222272,'38':0.8686641682}
    elif mod_type == 'MP':
        a = {'18': 0.05911087, '23': 0.1080751, '38': 0.37898495}
        b = {'18': 1.08693514, '23': 1.05342886, '38': 0.92876888}
    elif mod_type=='GAMMA':
        a = {'18': 0.04570854, '23': 0.08174184, '38': 0.28520923}
        b = {'18': 1.09211488, '23': 1.08105214, '38': 1.01426258}
    freq = int(freq/1e9)
    return a[str(freq)], b[str(freq)]