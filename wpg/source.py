#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 11 18:00:03 2024

@author: twguest
"""

import numpy as np

from wpg.wavefront import Wavefront

from phenom.source import sase_pulse as sp


def sase_pulse(x=None, y=None, t=None, photon_energy=10e3, pulse_energy=1e-03,
               pulse_duration=15e-15, bandwidth=1e-12, sigma=50e-06, div=2.5e-03,
               x0=0.0, y0=0.0, t0=0.0, theta_x=0.0, theta_y=0.0, domain='freq',
               polarization = 'horizontal'):
    """
    Define the SASE pulse with specified parameters or defaults.
    
    :param x: X-coordinates array [np.ndarray, optional]
    :param y: Y-coordinates array [np.ndarray, optional]
    :param t: Time array [np.ndarray, optional]
    :param photon_energy: Photon energy in eV [float, default=10e3]
    :param pulse_energy: Pulse energy in J [float, default=1e-03]
    :param pulse_duration: Pulse duration in seconds [float, default=15e-15]
    :param bandwidth: Bandwidth in Hz [float, default=1e-12]
    :param sigma: Beam size in meters [float, default=50e-06]
    :param div: Beam divergence in radians [float, default=2.5e-03]
    :param x0: Initial x position [float, default=0.0]
    :param y0: Initial y position [float, default=0.0]
    :param t0: Initial time position [float, default=0.0]
    :param theta_x: Beam tilt angle in x [float, default=0.0]
    :param theta_y: Beam tilt angle in y [float, default=0.0]
    :param domain: Domain type ('freq' or 'time') [str, default='freq']
    :param polarization: 'horizontal', 'vertical' or 'circular'
    :return: Electric field of the SASE pulse [np.ndarray]
    """
    efield = sp(x=x,
                y=y,
                t=t,
                photon_energy=photon_energy,
                pulse_energy=pulse_energy,
                pulse_duration=pulse_duration,
                bandwidth=bandwidth,
                sigma=sigma,
                div=div,
                x0=x0,
                y0=y0,
                t0=t0,
                theta_x=theta_x,
                theta_y=theta_y,
                domain = 'freq')
 
    x = x
    y = y
    t = t
    photon_energy = photon_energy
    
    
    nx, ny, nt = efield.shape

    wfr = Wavefront()


    # Setup E-field.
    wfr.data.arrEhor = np.zeros(shape=(nx, ny, nt, 2))
    wfr.data.arrEver = np.zeros(shape=(nx, ny, nt, 2))

    wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wfr.params.photonEnergy = photon_energy
    
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

    
    
    return wfr