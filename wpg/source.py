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
    :param bandwidth: relative bandwidth, fraction of photon_energy (as in phenom) [float or list/array, default=1e-12]
    :param sigma: Beam size in meters [float or list/array, calculated if None]
    :param div: Beam divergence in radians [float or list/array, calculated if None]
    :param x0: Initial x position [float or list/array, default=0.0]
    :param y0: Initial y position [float or list/array, default=0.0]
    :param t0: Initial time position [float or list/array, default=0.0]
    :param theta_x: Beam tilt angle in x [float or list/array, default=0.0]
    :param theta_y: Beam tilt angle in y [float or list/array, default=0.0]
    :param domain: Domain type ('freq' or 'time') [str, default='freq']
    :param polarization: 'horizontal', 'vertical' or 'circular' (Ey = i*Ex, each component
        with 1/2 of the pulse energy)
    :return: wpg.Wavefront, or list of them if any parameter is a list/array
    """
    if x is None or y is None or t is None:
        raise ValueError('sase_pulse needs the mesh: x, y and t coordinate arrays must be given')
    if domain not in ('time', 'freq'):
        raise ValueError("domain must be 'time' or 'freq', given: {!r}".format(domain))
    if polarization not in ('horizontal', 'vertical', 'circular'):
        raise ValueError("polarization must be 'horizontal', 'vertical' or 'circular', given: {!r}".format(
            polarization))

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
            
        # phenom builds the transverse field with np.meshgrid(x, y) (indexing 'xy'),
        # so the array is (ny, nx, nt) -- the layout of WPG arrEhor/arrEver
        ny, nx, nt = efield.shape
        if (nx, ny, nt) != (len(x), len(y), len(t)):
            raise ValueError('phenom field of shape {} does not match (len(y), len(x), len(t)) = {}'.format(
                efield.shape, (len(y), len(x), len(t))))

        # complex field -> WPG (Re, Im) arrays of both polarizations
        if polarization == 'horizontal':
            e_hor, e_ver = efield, np.zeros_like(efield)
        elif polarization == 'vertical':
            e_hor, e_ver = np.zeros_like(efield), efield
        else:  # circular: Ey = i*Ex, total energy as for a linear polarization
            e_hor, e_ver = efield / np.sqrt(2), 1j * efield / np.sqrt(2)

        wfr = Wavefront()
        # domain flag and slice mesh in the units of the phenom field (no FFT here)
        if domain == 'time':
            wfr.params.wDomain = 'time'
            slice_min, slice_max = np.min(t), np.max(t)
        else:
            # phenom (spectrum.linear_SASE_spectrum) puts the spectrum on
            # E = linspace(E0 - 3*bw, E0 + 3*bw, nt), bw = E0*bandwidth
            wfr.params.wDomain = 'frequency'
            bw = photon_energy[i] * bandwidth[i]
            slice_min, slice_max = photon_energy[i] - 3 * bw, photon_energy[i] + 3 * bw
        wfr.params.wEFieldUnit = 'sqrt(W/mm^2)'
        wfr.params.photonEnergy = photon_energy[i]

        wfr.params.Mesh.nSlices = nt
        wfr.params.Mesh.nx = nx
        wfr.params.Mesh.ny = ny
        wfr.params.Mesh.sliceMin = slice_min
        wfr.params.Mesh.sliceMax = slice_max
        wfr.params.Mesh.xMin = np.min(x)
        wfr.params.Mesh.xMax = np.max(x)
        wfr.params.Mesh.yMin = np.min(y)
        wfr.params.Mesh.yMax = np.max(y)

        wfr.params.Rx = 1
        wfr.params.Ry = 1

        wfr.data.arrEhor = np.stack([e_hor.real, e_hor.imag], axis=-1)
        wfr.data.arrEver = np.stack([e_ver.real, e_ver.imag], axis=-1)

        wavefronts.append(wfr)
    
    if len(wavefronts) == 1:
        wavefronts = wavefronts[0]
        
    return wavefronts
