#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 18:00:03 2024

@author: twguest
"""

import numpy as np

from wpg.wavefront import Wavefront

from phenom.source import sase_pulse as sp

import scipy.constants

h = scipy.constants.physical_constants['Planck constant in eV s'][0]

def analytical_pulse_energy(q, photon_energy):
    """
    Estimate of analytical_pulse_energy from electron bunch charge and radiation energy

    :param q: electron bunch charge [nC]
    :param photon_energy: radiation energy [eV]

    :return P: pulse energy [J]
    """

    P = 19*q/photon_energy
    return P

def analytical_pulse_duration(q):
    """
    Estimate analytical_pulse_duration from electron bunch charge

    :param q: electron bunch charge [nC]

    :return t: Duration of pulse [s]
    """

    t = (q*1e3)/9.8
    return t*1e-15


def analytical_pulse_width(photon_energy):
    """
    Estimate analytical_pulse_width (FWHM) from radiation energy (assumes symmetrical beam)

    :param photon_energy: radiation energy [eV]

    :return sig: Radiation pulse width [m]
    """

    sig = np.log((7.4e03/(photon_energy/1e03)))*6
    return sig/1e6


def analytical_pulse_divergence(photon_energy):

    """
    Estimate of analytical_pulse_divergence (half-angle) from electron bunch charge and radiation energy

    :param q: electron bunch charge [nC]
    :param photon_energy: radiation energy [eV]

    :return dtheta: pulse divergence [rad]
    """
    return ((14.1)/((photon_energy/1e03)**0.75)) / 1e06

def sase_pulse(x=None, y=None, t=None, photon_energy=10e3, pulse_energy=1e-03,
               pulse_duration=15e-15, bandwidth=1e-12, sigma=None, div=None,
               x0=0.0, y0=0.0, t0=0.0, theta_x=0.0, theta_y=0.0, domain='freq',
               polarization='horizontal'):
    """
    Define the SASE pulse with specified parameters or defaults.

    :param x: X-coordinates array [np.ndarray, optional]
    :param y: Y-coordinates array [np.ndarray, optional]
    :param t: Time array [np.ndarray, optional]
    :param photon_energy: Photon energy in eV [float or list/array, default=10e3]
    :param pulse_energy: Pulse energy in J [float or list/array, default=1e-03]
    :param pulse_duration: Pulse duration in seconds [float or list/array, default=15e-15]
    :param bandwidth: Bandwidth in eV [float or list/array, default=1e-12]
    :param sigma: Beam size in meters [float or list/array, calculated if None]
    :param div: Beam divergence in radians [float or list/array, calculated if None]
    :param x0: Initial x position [float or list/array, default=0.0]
    :param y0: Initial y position [float or list/array, default=0.0]
    :param t0: Initial time position [float or list/array, default=0.0]
    :param theta_x: Beam tilt angle in x [float or list/array, default=0.0]
    :param theta_y: Beam tilt angle in y [float or list/array, default=0.0]
    :param domain: Domain type ('freq' or 'time') [str, default='freq']
    :param polarization: 'horizontal', 'vertical' or 'circular'
    :return: List of electric fields of the SASE pulses [list of np.ndarray]
    """
    
    # Determine the maximum length of the input lists/arrays
    arg_lengths = [len(arg) if isinstance(arg, (list, np.ndarray)) else 1 for arg in 
                   [photon_energy, pulse_energy, pulse_duration, bandwidth, sigma, div, x0, y0, t0, theta_x, theta_y]]
    max_length = max(arg_lengths)
    
    # Convert all scalar arguments to lists of the appropriate length
    def to_list(arg):
        if isinstance(arg, (list, np.ndarray)):
            return arg
        else:
            return [arg] * max_length
    
    photon_energy = to_list(photon_energy)
    pulse_energy = to_list(pulse_energy)
    pulse_duration = to_list(pulse_duration)
    bandwidth = to_list(bandwidth)
    sigma = to_list(sigma) if sigma is not None else [analytical_pulse_width(e) for e in photon_energy]
    div = to_list(div) if div is not None else [analytical_pulse_divergence(e) for e in photon_energy]
    x0 = to_list(x0)
    y0 = to_list(y0)
    t0 = to_list(t0)
    theta_x = to_list(theta_x)
    theta_y = to_list(theta_y)
    
    # Generate the wavefronts
    wavefronts = []
    for i in range(max_length):
        efield = sp(x=x,
                    y=y,
                    t=t,
                    photon_energy=photon_energy[i],
                    pulse_energy=pulse_energy[i],
                    pulse_duration=pulse_duration[i],
                    bandwidth=bandwidth[i],
                    sigma=sigma[i],
                    div=div[i],
                    x0=x0[i],
                    y0=y0[i],
                    t0=t0[i],
                    theta_x=theta_x[i],
                    theta_y=theta_y[i],
                    domain=domain)
            
        nx, ny, nt = efield.shape

        wfr = Wavefront()

        # Setup E-field.
        wfr.data.arrEhor = np.zeros(shape=(nx, ny, nt, 2))
        wfr.data.arrEver = np.zeros(shape=(nx, ny, nt, 2))

        wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
        wfr.params.photonEnergy = photon_energy[i]
        
        wfr.params.Mesh.nSlices = nt
        wfr.params.Mesh.nx = nx
        wfr.params.Mesh.ny = ny      
        
        wfr.params.Mesh.sliceMin = np.min(t)
        wfr.params.Mesh.sliceMax = np.max(t)
        
        wfr.params.wDomain = 'time'
        
        wfr.set_electric_field_representation('f')
        
        wfr.params.Mesh.xMin = np.min(x)
        wfr.params.Mesh.xMax = np.max(x)
        wfr.params.Mesh.yMin = np.min(y)
        wfr.params.Mesh.yMax = np.max(y)

        wfr.params.Rx = 1
        wfr.params.Ry = 1
        
        #convert complex wavefield into wpg style electric field array
        arrE = np.zeros([efield.shape[0], efield.shape[1], efield.shape[2], 2])
        
        arrE[:,:,:,0] = efield.real
        arrE[:,:,:,1] = efield.imag
        
        # polarization
        if polarization == 'horizontal':
            wfr.data.arrEhor = arrE
            wfr.data.arrEver = np.zeros_like(arrE)
        elif polarization == 'vertical':
            wfr.data.arrEhor = None
            wfr.data.arrEver = np.zeros_like(arrE)
        elif polarization == 'circular':
            wfr.data.arrEhor = arrE
            wfr.data.arrEver = arrE

        wavefronts.append(wfr)    
    
    if len(wavefronts) == 1:
        wavefronts = wavefronts[0]
        
    return wavefronts
