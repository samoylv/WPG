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
    :return: Tuple containing the optimal parameters for the Gaussian fit and the covariance of the parameters. [tuple (np.ndarray, np.ndarray)]
    """
    params, params_covariance = curve_fit(gaussian, axis, data, **kwargs)
    return params, params_covariance
