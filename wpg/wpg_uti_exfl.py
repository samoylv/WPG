# -*- coding: utf-8 -*-

__author__ = 'A. Buzmakov, L. Samoylova'

import time
import numpy as np

def calculate_theta_fwhm_cdr_s1(ekev,qnC):
    """
    Calculate angular divergence using empiric formula from XFEL CDR2011 for SASE1 undulator
    
    :param ekev: Energy in keV
    :param qnC: e-bunch charge, [nC]
    :return: theta_fwhm [rad]
    """
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC))*1e-6/ekev**0.85
    return theta_fwhm

def calculate_theta_fwhm_cdr_s3(ekev,qnC):
    """
    Calculate angular divergence using empiric formula from XFEL CDR2011 for SASE3 undulator
    
    :param ekev: Energy in keV
    :param qnC: e-bunch charge, [nC]
    :return: theta_fwhm [rad]
    """
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC))*1e-6/ekev**0.85
    return theta_fwhm


