# -*- coding: utf-8 -*-

import numpy as np
from scipy.optimize import curve_fit

def gaussian(x, sigma, x0=0, a=1):
    """
    One-dimensional Gaussian function for fitting.
    
    This function defines a Gaussian curve, commonly used for fitting data
    where the Gaussian distribution is expected.
    
    :param x: The input array of x-values where the Gaussian function is evaluated. [np.ndarray]
    :param sigma: The standard deviation of the Gaussian distribution. [float]
    :param x0: The mean or center of the Gaussian distribution (default is 0). [float, optional]
    :param a: The amplitude of the Gaussian peak (default is 1). [float, optional]
    :return: The calculated Gaussian function values for the given x. [np.ndarray]
    """
    return a * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def fit_gaussian(axis, data, **kwargs):
    """
    Optimization for fitting a Gaussian function to data.
    
    This function uses the scipy.optimize.curve_fit method to fit a Gaussian
    function to the provided data. Additional keyword arguments can be passed
    to curve_fit for more control over the fitting process.
    
    :param axis: The x-values (independent variable) of the data. [np.ndarray]
    :param data: The y-values (dependent variable) of the data to be fit. [np.ndarray]
    :param kwargs: Additional keyword arguments for scipy.optimize.curve_fit.
        Without p0 the start point (sigma, x0, a) is estimated from the moments of data.
    :return: Tuple containing the optimal parameters for the Gaussian fit and the covariance of the parameters. [tuple (np.ndarray, np.ndarray)]
    """
    if 'p0' not in kwargs:
        # curve_fit starts from (1, 1, 1) and stops there for data in metres
        w = np.clip(np.asarray(data, dtype=float), 0, None)
        x0 = np.sum(axis * w) / np.sum(w)
        kwargs['p0'] = [np.sqrt(np.sum((axis - x0) ** 2 * w) / np.sum(w)), x0, np.max(data)]
    params, params_covariance = curve_fit(gaussian, axis, data, **kwargs)
    return params, params_covariance
