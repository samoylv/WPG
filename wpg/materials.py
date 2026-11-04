# -*- coding: utf-8 -*-
"""
Optical constants of materials and reflectivity of X-ray mirrors.

Refractive index ``n = 1 - delta + i*beta`` is taken from xraylib (optional dependency,
``pip install wpg[materials]``; T. Schoonjans et al., Spectrochim. Acta B 66, 776 (2011)).
Mirror reflectivity of a thick single-layer coating is calculated with the Fresnel equations,
roughness is taken into account with the Nevot-Croce factor (L. Nevot, P. Croce,
Rev. Phys. Appl. 15, 761 (1980)) - the same model as the CXRO mirror calculator
(https://henke.lbl.gov/optical_constants/mirror2.html).

The resulting :class:`Reflectivity` is applied to SRW mirrors (``SRWLOptMir*.set_reflect``):
SRW multiplies the sigma (perpendicular to the plane of incidence) and pi components of the field
by the complex amplitude reflectivity at the local grazing angle and photon energy, taking the
nearest node of a uniform table (no interpolation).

.. module:: wpg.materials
   :platform: Linux, Mac OSX, Windows
"""
from array import array

import numpy as np

#: h*c [eV*m], CODATA 2018
HC_EV_M = 1.239841984e-6


def _xraylib():
    try:
        import xraylib
    except ImportError as e:
        raise ImportError('Optical constants of materials need xraylib: '
                          'pip install xraylib (or pip install wpg[materials])') from e
    return xraylib


def density(material):
    """
    Tabulated density of an element ('Pt') or of a NIST compound ('Water, Liquid') [g/cm^3].

    :raises ValueError: for a chemical formula of a compound ('B4C'): give its density explicitly
    """
    xrl = _xraylib()
    try:
        return xrl.ElementDensity(xrl.SymbolToAtomicNumber(material))
    except ValueError:
        pass
    try:
        return xrl.GetCompoundDataNISTByName(material)['density']
    except ValueError:
        raise ValueError('No tabulated density for "{}": give the density in g/cm^3 '
                         '(e.g. B4C 2.52)'.format(material)) from None


def refractive_index(material, ekev, rho=None):
    """
    Complex refractive index ``n = 1 - delta + i*beta`` (xraylib).

    :param material: chemical formula ('B4C', 'Pt', 'SiO2') or NIST compound name
    :param ekev: photon energy [keV], number or array
    :param rho: density [g/cm^3], None - tabulated density (elements and NIST compounds only)
    :return: complex, the shape of ekev
    """
    xrl = _xraylib()
    rho = density(material) if rho is None else float(rho)
    e = np.asarray(ekev, dtype=float)
    n = np.array([complex(xrl.Refractive_Index_Re(material, v, rho),
                          xrl.Refractive_Index_Im(material, v, rho)) for v in e.ravel()]).reshape(e.shape)
    return n[()] if n.ndim == 0 else n


def delta_beta(material, ekev, rho=None):
    """(delta, beta) of ``n = 1 - delta + i*beta``, see :func:`refractive_index`."""
    n = refractive_index(material, ekev, rho)
    return 1. - np.real(n), np.imag(n)


def attenuation_length(material, ekev, rho=None):
    """
    Intensity attenuation length ``lambda / (4 pi beta)`` [m], as used by CRL and transmission elements.

    xraylib takes beta from the total cross-section (photoabsorption, Rayleigh and Compton scattering).
    CXRO (Henke) tables give beta of photoabsorption only, so their attenuation length is longer for
    light materials: Be at 8 keV 4.8 mm here, 6.1 mm with photoabsorption only (both from xraylib
    cross-sections). The CRL examples of WPG take attenuation lengths from Henke data.
    """
    _, beta = delta_beta(material, ekev, rho)
    return HC_EV_M / (np.asarray(ekev, dtype=float) * 1e3) / (4 * np.pi * beta)


def fresnel_reflectivity(n, theta, ekev=None, roughness=0.):
    """
    Complex amplitude reflectivity of a thick mirror (single interface vacuum - material).

    :param n: refractive index of the material, complex (broadcasts with theta)
    :param theta: grazing angle [rad]
    :param ekev: photon energy [keV], needed only for roughness
    :param roughness: RMS roughness [m] (Nevot-Croce factor)
    :return: (r_sigma, r_pi) - for the field perpendicular to and in the plane of incidence
    """
    n = np.asarray(n, dtype=complex)
    theta = np.asarray(theta, dtype=float)
    sin_t = np.sin(theta)
    n2 = n * n
    kz2 = np.sqrt(n2 - np.cos(theta) ** 2)  # normal component of the wave vector in the material / k
    kz2 = np.where(np.imag(kz2) < 0, -kz2, kz2)  # wave decaying into the material
    r_sigma = (sin_t - kz2) / (sin_t + kz2)
    r_pi = (n2 * sin_t - kz2) / (n2 * sin_t + kz2)
    if roughness:
        if ekev is None:
            raise ValueError('Photon energy is needed for roughness')
        k = 2 * np.pi * np.asarray(ekev, dtype=float) * 1e3 / HC_EV_M
        factor = np.exp(-2 * k ** 2 * roughness ** 2 * sin_t * kz2)
        r_sigma, r_pi = r_sigma * factor, r_pi * factor
    return r_sigma, r_pi


def _uniform_grid(values, name):
    v = np.atleast_1d(np.asarray(values, dtype=float))
    if v.ndim != 1:
        raise ValueError('{} should be a number or a 1D array'.format(name))
    if v.size > 1:
        step = np.diff(v)
        if not np.all(step > 0) or not np.allclose(step, step[0], rtol=1e-6, atol=0):
            raise ValueError('{} should be increasing with a constant step: SRW takes the nearest '
                             'node of a uniform table'.format(name))
    return v


class Reflectivity(object):
    """
    Complex amplitude reflectivity of a mirror vs photon energy and grazing angle.

    SRW takes the value of the nearest node (no interpolation), so the grids should be fine
    compared with the variation of the reflectivity.

    :param amplitude: complex amplitude reflectivity, shape (n_energies, n_angles) - the same for
        sigma and pi - or (n_energies, n_angles, 2) - sigma, pi
    :param ekev: photon energies [keV], uniform grid or one number
    :param angles: grazing angles [rad], uniform grid or one number
    """

    def __init__(self, amplitude, ekev, angles):
        self.ekev = _uniform_grid(ekev, 'Photon energies')
        self.angles = _uniform_grid(angles, 'Grazing angles')
        amp = np.asarray(amplitude, dtype=complex)
        if amp.ndim < 3:
            amp = amp.reshape(self.ekev.size, self.angles.size, 1)
        if amp.shape[:2] != (self.ekev.size, self.angles.size) or amp.shape[2] not in (1, 2):
            raise ValueError('Amplitude shape {} does not match {} energies x {} angles (x 1 or 2 '
                             'components)'.format(amp.shape, self.ekev.size, self.angles.size))
        self.amplitude = amp
        #: material parameters for Reflectivity.from_material, None for a table given by the user
        self.material = None

    @classmethod
    def constant(cls, R):
        """Reflectivity R (fraction of the reflected intensity) for all energies, angles and polarizations."""
        if not 0 <= R <= 1:
            raise ValueError('Intensity reflectivity should be in [0, 1], received {}'.format(R))
        return cls(np.sqrt(R), 0., 0.)

    @classmethod
    def from_material(cls, material, ekev, angles, rho=None, roughness=0.):
        """
        Fresnel reflectivity of a thick coating with optical constants from xraylib.

        :param material: chemical formula ('B4C', 'Pt') or NIST compound name
        :param ekev: photon energies [keV], uniform grid (e.g. ``np.linspace(8.9, 9.1, 101)``) or one number
        :param angles: grazing angles [rad], uniform grid covering the local angles on the mirror
            (e.g. ``np.linspace(0.8 * theta, 1.2 * theta, 201)``) or one number
        :param rho: density [g/cm^3], None - tabulated (elements and NIST compounds)
        :param roughness: RMS roughness [m]
        """
        e = _uniform_grid(ekev, 'Photon energies')
        a = _uniform_grid(angles, 'Grazing angles')
        n = refractive_index(material, e, rho)[:, None]
        r_sigma, r_pi = fresnel_reflectivity(n, a[None, :], e[:, None], roughness)
        refl = cls(np.stack([r_sigma, r_pi], axis=-1), e, a)
        refl.material = dict(material=material, rho=rho, roughness=roughness)
        return refl

    def intensity(self):
        """|amplitude|^2, shape (n_energies, n_angles, n_components)."""
        return np.abs(self.amplitude) ** 2

    def apply(self, mirror):
        """Set this reflectivity to an SRW mirror (SRWLOptMir*); returns the mirror."""
        amp, ekev, angles = self.amplitude, self.ekev, self.angles
        # SRW computes the node index dividing by the step, which is 0 for one node:
        # one node is given as two equal nodes
        if ekev.size == 1:
            amp, ekev = np.repeat(amp, 2, axis=0), np.array([ekev[0], ekev[0] + 1.])
        if angles.size == 1:
            amp, angles = np.repeat(amp, 2, axis=1), np.array([angles[0], angles[0] + 1e-3])
        n_e, n_a, n_c = amp.shape
        # SRW layout: complex values, photon energy changes fastest, then angle, then component
        flat = np.empty(2 * amp.size)
        amp = np.transpose(amp, (2, 1, 0)).ravel()
        flat[0::2], flat[1::2] = np.real(amp), np.imag(amp)
        mirror.set_reflect(array('d', flat), n_e, n_a, n_c,
                           1e3 * ekev[0], 1e3 * ekev[-1], 'lin', angles[0], angles[-1], 'lin')
        return mirror

    def __repr__(self):
        R = self.intensity()
        return 'Reflectivity({} energies {:g}-{:g} keV, {} angles {:g}-{:g} rad, R {:.3g}-{:.3g})'.format(
            self.ekev.size, self.ekev[0], self.ekev[-1], self.angles.size, self.angles[0],
            self.angles[-1], R.min(), R.max())
