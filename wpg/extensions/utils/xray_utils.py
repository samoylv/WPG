# -*- coding: utf-8 -*-

import numpy as np

#%%

E = .18505 # keV corresponds to 6.7 nm

def E2wavelength(E):
    
    """
    Converts energy in eV to wavelength in m
    
    :param E: Energy (in eV)
    """
    
    h = 4.135667662e-15
    c = 299792458
    
    wavelength = (h*c)/E
    
    return wavelength

def wavelength2E(wavelength):
    
    """
    Converts wavelength in m to energy in eV
    
    :param wavelength: wavelength (in m)
    """
    
    h = 4.135667662e-15
    c = 299792458
    
    E = (h*c)/wavelength
    
    return E


def theta_fwhm(ekev,qnC):
    
    """
    Calculates the theta fwhm of a gaussian electron beam
    
    param: ekev: energy of beam in keV
    param qnC: charge of electron beam in nC
    """
    
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC))*1e-6/ekev**0.85
    return theta_fwhm


