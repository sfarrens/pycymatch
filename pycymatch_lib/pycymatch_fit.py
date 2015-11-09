#  @file pycymatch_fit.py
#
#  PYCYMATCH FIT METHODS
#
#  Functions for fitting
#  p_lambda distributions.
#
#  @author Samuel Farrens
#  @version 1.0
#  @date 2015
#

import numpy as np
from numpy.lib.recfunctions import append_fields
from scipy.optimize import curve_fit


##
#  Gaussian function.
#
#  @param[in] x: Input data.
#  @param[in] amp: Amplitude.
#  @param[in] mean: Mean of distribution.
#  @param[in] sigma: Standard deviation.
#
#  @return Gaussian probability.
#
def gauss(x, amp, mean, sigma):

    with np.errstate(divide='ignore'):
        return amp * np.exp(-(x - mean) ** 2 / (2 * sigma ** 2))


##
#  Function to find best Gaussian fit to curve.
#
#  @param[in] x: X-axis data.
#  @param[in] y: Y-axis data.
#
#  @return Best fit parameters.
#
def fit(x, y):

    if max(y) > 0.0:
        mean = x[y == max(y)]

    else:
        mean = [np.sum(x * y) / len(x)]

    sigma = np.sqrt(sum(y * (x - mean[0]) ** 2) / len(x))

    if sigma < 0.0:
        sigma = 0.0

    try:
        popt, pcov = curve_fit(gauss, x, y, p0=[0, mean[0], sigma])
        return popt

    except:
        return [0., 0., 0.]


##
#  Function to find Gaussian fits to
#  P_lambda distributions.
#
#  @param[in] matrix: mass-observable
#  matrix.
#  @param[in] x_range: X-axis data.
#
#  @return Best fit parameters.
#
def get_fits(matrix, x_range):

    fit_param = []

    for i in range(len(matrix)):
        fit_param.append(fit(x_range, matrix[i]))

    return np.array(fit_param)
