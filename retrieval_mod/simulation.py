#-------------------------------------------------------------------------------
# Name: simulation.py
# Purpose: This module provides simulation class for rainfall field retrieval.
#
# Information for one link:
#     one_link = [ "ID" , "Caf",  "Frequency",  "Recx", "Recy", "Trx",
#                  "Try",  "RecAzimuth", "TrAzimuth", "RecPlace", "TrPlace",
#                  "RecCompany", "TrCompany"]
#-------------------------------------------------------------------------------




import numpy as np
import copy
import os
from copy import deepcopy
from pandas import Series, DataFrame

# my lib
import rain_atten.attenuation_model as att_model
from antenna_processing.link_process import get_link_length
from inverse_technique.nestingmethod import inverse_nesting_grid
from forward_process import discretize
import metrics
import uncertainty




def create_folder(newpath):
    '''
    Creates new path.
    If "newpath"  already exists then it return just its name.
    If "newpath" does not exist then it is created.
    Inputs:
        newpath - string. E.g.: newpath = "data/my_new_path"
    Output:
        newpath name - string.
    '''
    if not os.path.exists(newpath):
            os.makedirs(newpath)
    return newpath


def add_powerlaw(link, frequency, a_coef, b_coef):
    '''
    Returns given *link* with associated a_coef, b_coef columns at frequency.
    Inputs:
        link       -  DataFrame. Link data.
        frequency  -  list.      Antenna frequency ranges, [Hz]
        a_coef     -  list.      Power law coefficient for given frequency(s).
                                 By default, a_coef = None.
        b_coef     -  list.      Power law coefficients for given frequency(s).
                                 By default, b_coef = None.
    '''
    # To prevent mess up in system
    clink = deepcopy(link)
    clink['a_coef'] = Series(np.zeros(len(clink)), index=clink.index)
    clink['b_coef'] = Series(np.zeros(len(clink)), index=clink.index)
    # Give a_coef, b_coef value to respective Frequency
    for f, i, j in zip(frequency, a_coef, b_coef):
        number = link[ link['Frequency'] == f ].index
        clink['a_coef'][number] = [i] * len(number)
        clink['b_coef'][number] = [j] * len(number)
    return clink['a_coef'], clink['b_coef']


#TODO documentation should be improved in each function of the class
class MonitorArea:

    def __init__(self, link, frequency, border, a_coef, b_coef):
        '''
        Initialization parameters.
        Inputs:
            link       - DataFrame.  Link information.
            frequency  - 1D list.    Range of antenna frequencies, [Hz].
            border     - (2,2)tuple. Area border: ((x0, x1),(y0, y1)), [km].
        '''
        self.link = link
        self.frequency = frequency
        self.border = border
        self.link['a_coef'], self.link['b_coef'] = add_powerlaw( link,
                                                            frequency,
                                                            a_coef, b_coef )
        self.link['Length'] = get_link_length(link)
        self.lij = None
        self.pij = None
        self.nesting_lij = []
        self.nesting_pij = []

        # Retrieval variables
        self.retrieval_rain = None
        self.retrieval_iter = None
        self.retrieval_phi = None
        self.retrieval_p0 = None

    def set_rsl_discretization( self, resolx=0.25, resoly=0.25 ):
        '''
        Sets discretization for RSL at a given discretization step.
        Inputs:
            resolx  - scalar. Discretization in x axis, [km].
                              By default, resolx = 0.25
            resoly  - scalar. Discretization in y axis, [km].
                              By default, resoly = 0.25
        '''
        # discretize the area with (resolx,resoly) sq.km resolution
        self.lij, self.pij = discretize( self.border, resolx, resoly,
                                         self.link['Trx'],  self.link['Try'],
                                         self.link['Recx'], self.link['Recy'])

    def measure_att( self, rainmap, model='empirical'):
        '''
        Returns Received Signal Level (RSL) along the links.
        Inputs:
            rainmap  - 2D array. Rainfall map.
        Optional:
            model    - string. Attenuation model to compute signal reduce along
                               the link. It can be 'empirical' or 'theoretical'.
                               By default, model = 'empirical'.
        Output:
            total_rsl  - 1D array. RSL for given link(s).
            std_alpha  - 1D array. Standard deviation of measurement noise.
        '''
        rvector = rainmap.flatten()
        # measurement without noise
        if model == 'empirical':
            truedata = [ att_model.empirical(a, b, rvector[p], l)
                            for a, b, p, l in zip( self.link['a_coef'],
                                                   self.link['b_coef'],
                                                   self.pij, self.lij )]
        elif model == 'theoretical':
            truedata = [ sum([ att_model.theoretical( freq, rvector[p], l)
                         for p, l in zip(pix_, len_)])
                             for freq, pix_, len_ in zip( self.link['Frequency'],
                                                          self.pij, self.lij )]
        return np.array(truedata)

    def add_noise2rsl(self, rsl, alpha=0, delta=0, zero_att=1e-3):
        '''
        Adds noise to a given rsl measurement.
        Inputs:
            rsl - 1D array. Attenuation measurement.
        Optional:
            alpha    - scalar. Uncertainty of total attenuation measurement.
                               *alpha* takes value at interval [0, 100), [%]
                               According to Guili et al, 1991. This is the
                               certain percentage of measured total RSL.
                               By default, alpha=0 which means no uncertainty.
            delta    - scalar. Quantization step [dB]. This can be 0, 0.1 or 1.
                               By default, delta=0 which means no quantization.
            zero_att - scalar. Minimal value to avoid negative attenuation value
                               in rsl measurement. By default, 1e-3.
        Output:
            rsl_noisy - 1D array. Noisy attenuation measurement.
        '''
        # noisy measurement
        sigma_meas = alpha * rsl
        noise = uncertainty.general(sigma_meas)
        noisydata = rsl + noise
        # quantize the signal
        total_rsl = uncertainty.quantize(noisydata, delta)
        # To prevent "0" attenuation
        total_rsl[total_rsl <= 0] = zero_att
        return total_rsl

    def set_retrieval_discretization(self, dx=[2., 1., 0.5], dy=[2., 1., 0.5]):
        l = []
        p = []
        for x, y in zip(dx, dy):
            d = discretize( self.border, x, y, self.link['Trx'],  self.link['Try'],
                                               self.link['Recx'], self.link['Recy'])
            l.append(d[0])
            p.append(d[1])
        self.nesting_lij = l
        self.nesting_pij = p

    #TODO documentation should be completed
    def reconstruct(self, rsl, dx, dy, apriori = ['closest_links_based', None],
                       decorrelation = 4, delta = 0., beta = 0.,
                       global_prior_percent = None, iter_max = 20, accuracy = 0.01 ):
        '''
        Rainfall retieval for a given RSL measurement.

        Inputs:
            rsl     - 1D array.    RSL measurement vector.
            dx      - list.        Resolutions for retrieval in x axis.
            dy      - list.        Resolutions for retrieval in y axis.

        Optional:
            apriori - (2L,) list. Method for computing apriori rain:
                                   apriori(name, parameter)
                                   First element *name* is string,
                                   second one *parameter* is scalar.
                                   By default, apriori = ('limit_weight_based', None)
                                  Available methods:
                                  1) ('closest_links_based', parameter)
                                  2) ('constant_radius_based', parameter)
                                  3) ('limit_weight_based', parameter)
                                  where, parameter is as follows, respectively:
                                  1) link_number- scalar. Number of links,
                                     It can be within interval [1, M],
                                     M - maximum number of links in *border*.
                                  2) radius - scalar. Circle radius,
                                     It can be within interval (0, R],
                                     R - diagonal length of area *border*
                                  3) min_weight - scalar. Minimum weight for
                                     crossed part of links. It can be within
                                     interval (0, 1].
            decorrelation - scalar.   Decorrelation distance, [km].
                                      By default, decorr = 4
            delta         - scalar.   Quantization step [dB].This is 0, 0.1 or 1.
                                      By default, delta=0 which means no quantization.
            beta          - scalar.   Model error caused by the power law (A~R)
                                      relation. It is the certain percentage
                                      of measured total RSL.
                                      By default, beta = 0.
            var_percent   - scalar.   Percentage of variance of global apriori
                                      is equal to variance of local apriori value.
                                      By default, var_percent = None
            second_prior  - boolean.  If True, second prior is enabled and its
                                      value is taken as in first prior.
                                      If False, second prior is disabled.
                                      By default, second_prior = False
            iter_max      - scalar.   Maximum iteration number for retrieval.
                                      By default, iter_max = 20
            accuracy      - scalar.   Relative error between two consequent
                                      solutions. Iteration stops if accuracy is
                                      reached.                                      .
                                      By default, accuracy = 0.01
        Outputs:
            solu          - (3L,) list. Retrieved rainfall map at (dx, dy) resolutions
            r0            - list.     Apriori rain vector at (dx,dy) resolutions
            var_rsl       - 1D array. RSL data variance.
            var_model     - 1D array. Model error variance.
            var_quant     - 1D array. Quantization error variance.
            total_var     - 1D array. Total error variance.
                                      total_var = var_rsl + var_model + var_quant
        '''
        # RSL measurement variance
        var_rsl = np.var(rsl) * np.ones(len(rsl))
        # Model error variance
        var_model = (beta * rsl)**2
        # Quantization error variance
        var_quant = delta**2/12. * np.ones(len(rsl))
        # Total variance var(rsl) + var(model) + var(quant)
        total_var = var_rsl + var_model + var_quant
        # Retrieval part
        #TODO "p0" apriori value ham qaytarilishi kerak
        #     redisual ni hisoblash uchun kerak buladi
        solu, ite,\
        phi, p0 = inverse_nesting_grid( rsl, total_var, apriori,
                                        decorrelation,
                                        self.border, dx, dy,
                                        self.nesting_lij, self.nesting_pij,
                                        self.link['Trx'], self.link['Try'],
                                        self.link['Recx'], self.link['Recy'],
                                        self.link['Length'], self.link['a_coef'],
                                        self.link['b_coef'],
                                        global_prior_percent,
                                        iter_max, accuracy )
        #variance = np.vstack([var_rsl, var_model, var_quant, total_var]).T
        # TODO variance ni yam output ga chiqarish kerak after SENSITIVITY
        self.retrieval_rain = solu
        self.retrieval_iter = ite
        self.retrieval_phi = phi
        self.retrieval_p0 = p0
        return solu

    #TODO "stats" variable should be extented by adding other statistics.
    def evaluate(self, observed, estimated, zone, res = 0.5):
        '''
        Returns evaluation between observed and estimated rainfall maps
        by computing following metrics:
            ['BIAS', 'CORREALTION', 'Nash-Sutcliffe', 'RMSE', 'MAE',
                                                      'MEAN_OBS','MEAN_EST']
        Inputs:
            observed  -  2D array.    Observed rainfall map.
            estimated -  2D array.    Estimated rainfall map.
            zone      -  (2,2) tuple. Evaluation study zone [km x km]
        Optional:
            res       -  scalar.      Resolution for comparison, [km].
        Outputs:
            metrics   -  dictionary.  Statistical metrics:
                                     ['bias', 'corr', 'nash', 'rmse', 'mae',
                                                      'mean_obs','mean_est']
        '''
        # We neglect area that is not estimated for comparison purpose.
        estimated[estimated <= -999] = -999
        observed[estimated <= -999] = -999

        ((x0, x1), (y0, y1)) = zone
        # Cut the zone for evaluation
        t1, t2, t3, t4 = int(y0/res), int(y1/res), int(x0/res), int(x1/res)
        observed = observed[t1:t2, t3:t4]
        estimated = estimated[t1:t2, t3:t4]

        est = estimated[estimated<>-999].flatten()
        obs = observed[observed<>-999].flatten()
        stats = dict()
        stats['bias'] = metrics.bias(obs, est)
        stats['corr'] = metrics.corr(obs, est)
        stats['nash'] = metrics.nash(obs, est)
        stats['rmse'] = metrics.rmse(obs, est)
        stats['mae'] = metrics.mae(obs, est)
        stats['mean_obs'] = metrics.average(obs)
        stats['mean_est'] = metrics.average(est)
        # additional metrics can be added
        ##stats['likelihood'] = metrics.likelihood(obs, est)
        ##stats['mape'] = metrics.mape(obs, est)
        ##stats['mse'] = metrics.mse(obs, est)
        ##stats['mspe'] = metrics.mspe(obs, est)
        ##stats['rmspe'] = metrics.rmspe(obs, est)
        return stats