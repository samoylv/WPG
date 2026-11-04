# -*- coding: utf-8 -*-
"""
This module contains base wrapper for SRWLWfr (Wavefront). It's implement numpy inter operations to SRWLWfr structure, serialization to HDF5, visualization tools, etc.

.. module:: wpg.wavefront
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

import array
import copy
import warnings

import numpy as np
import h5py

from srwpy import srwlpy 
import srwpy.srwlib as srwlib

import wpg.utils as utils
import wpg.glossary as glossary

from wpg.utils import srw_obj2str

from wpg.wpg_uti_math import fit_gaussian

# wpg.glossary warns on every load/store about fields not defined for the current
# wavefront (e.g. Q-space mesh); these messages are informational, hide only them
warnings.filterwarnings('ignore', category=UserWarning, module='wpg.glossary')


class Wavefront(object):
    """
    This is base class for manipulation with wavefronts in wpg module.

    One of most important field is _srwl_wf (instance of srwlib.SRWLWfr). SEtting and getting this field allows to call all SRWLpy functions.

    Wavefront fields (see :mod:`wpg.glossary`) are a live view of ``_srwl_wf``:
    ``wf.params`` (:class:`wpg.glossary.WavefrontParams`), ``wf.params.Mesh``
    (:class:`wpg.glossary.WavefrontMesh`) and ``wf.data`` (:class:`wpg.glossary.WavefrontData`).
    """

    params: glossary.WavefrontParams
    data: glossary.WavefrontData
    version = glossary.version_field

    # read-only shortcuts: wf.nx is wf.params.Mesh.nx, wf.photonEnergy is wf.params.photonEnergy
    nval = glossary.shortcut(glossary.WavefrontParams.nval)
    wSpace = glossary.shortcut(glossary.WavefrontParams.wSpace)
    wDomain = glossary.shortcut(glossary.WavefrontParams.wDomain)
    wFloatType = glossary.shortcut(glossary.WavefrontParams.wFloatType)
    wEFieldUnit = glossary.shortcut(glossary.WavefrontParams.wEFieldUnit)
    photonEnergy = glossary.shortcut(glossary.WavefrontParams.photonEnergy)
    Rx = glossary.shortcut(glossary.WavefrontParams.Rx)
    Ry = glossary.shortcut(glossary.WavefrontParams.Ry)
    dRx = glossary.shortcut(glossary.WavefrontParams.dRx)
    dRy = glossary.shortcut(glossary.WavefrontParams.dRy)
    xCentre = glossary.shortcut(glossary.WavefrontParams.xCentre)
    yCentre = glossary.shortcut(glossary.WavefrontParams.yCentre)
    nx = glossary.shortcut(glossary.WavefrontMesh.nx)
    ny = glossary.shortcut(glossary.WavefrontMesh.ny)
    nSlices = glossary.shortcut(glossary.WavefrontMesh.nSlices)
    zCoord = glossary.shortcut(glossary.WavefrontMesh.zCoord)
    xMin = glossary.shortcut(glossary.WavefrontMesh.xMin)
    xMax = glossary.shortcut(glossary.WavefrontMesh.xMax)
    yMin = glossary.shortcut(glossary.WavefrontMesh.yMin)
    yMax = glossary.shortcut(glossary.WavefrontMesh.yMax)
    qxMin = glossary.shortcut(glossary.WavefrontMesh.qxMin)
    qxMax = glossary.shortcut(glossary.WavefrontMesh.qxMax)
    qyMin = glossary.shortcut(glossary.WavefrontMesh.qyMin)
    qyMax = glossary.shortcut(glossary.WavefrontMesh.qyMax)
    sliceMin = glossary.shortcut(glossary.WavefrontMesh.sliceMin)
    sliceMax = glossary.shortcut(glossary.WavefrontMesh.sliceMax)
    hvx = glossary.shortcut(glossary.WavefrontMesh.hvx)
    hvy = glossary.shortcut(glossary.WavefrontMesh.hvy)
    hvz = glossary.shortcut(glossary.WavefrontMesh.hvz)
    nvx = glossary.shortcut(glossary.WavefrontMesh.nvx)
    nvy = glossary.shortcut(glossary.WavefrontMesh.nvy)
    nvz = glossary.shortcut(glossary.WavefrontMesh.nvz)

    @property
    def Mesh(self) -> glossary.WavefrontMesh:
        """Read-only shortcut of ``wf.params.Mesh``."""
        return self.params.Mesh


    def __init__(self, srwl_wavefront=None):
        """
        Create wavefront instance.

        :param srwl_wavefront: if present, wavefront inits with it's parameters
        :type srwl_wavefront: srwlib.SRWLWfr
        :return: Wavefront instance.
        """
        if srwl_wavefront is None:
            self._srwl_wf = srwlib.SRWLWfr()
        else:
            self._srwl_wf = srwl_wavefront

        self._version = 0.1  # HDF5 format version
        self._nval = 2  # complex electric field: (Re, Im)
        self.custom_fields = {}
        self.params = glossary.WavefrontParams(self)
        self.data = glossary.WavefrontData(self)

    @property
    def _wf(self):
        """The wavefront itself: glossary fields of the top level (version) are accessed via node._wf."""
        return self

    def _get_total_elements(self):
        """
        Get total amount of points in wavefront.

        :return: total amount of points in wavefront
        """
        return self.params.Mesh.nx * self.params.Mesh.ny * self.params.Mesh.nSlices

    def _allocate_srw_moments(self):
        """Allocate memory for SRW structures."""
        self._srwl_wf.arMomX = array.array(
            str(u'd'), [0] * self.params.Mesh.nSlices * 11)
        self._srwl_wf.arMomY = array.array(
            str(u'd'), [0] * self.params.Mesh.nSlices * 11)

    def _to_dict(self):
        """
        Convert wavefront to dictionary. Used for saving in HDF5 file.

        :return: dictionary view of wavefront
        """
        res = {f.path: glossary.field_value(self, f) for f in glossary.get_wf_fields()}
        res.update(self.custom_fields)
        return res

    def _update_from_dict(self, in_dict):
        """
        Update wavefront from dictionary. Used for loading wavefront to HDF5 file.

        Glossary fields are set in the order of :func:`wpg.glossary.get_wf_fields`
        (mesh size first, field arrays last), other keys go to custom_fields.

        :param in_dict: input dictionary
        :type in_dict: dict
        """
        values = {}
        for (key, value) in in_dict.items():
            if isinstance(key, bytes):
                key = key.decode('utf-8')
            values[key] = value

        for f in glossary.get_wf_fields():
            if f.path in values:
                f.__set__(glossary.field_node(self, f), values.pop(f.path))

        for (key, value) in values.items():
            utils.update_dict_slash_string(self.custom_fields, key, value)

    def _store_attributes(self, file_name):
        """
        Store wavefront attributes to HDF5 file.

        Attribute of each field is values of field.attributes

        :param file_name: output HDF5 file name
        :type  file_name: string
        """
        with h5py.File(file_name, 'r+') as h5f:
            for f in glossary.get_wf_fields():
                if f.path in h5f:
                    for (k, v) in f.attributes.items():
                        h5f[f.path].attrs[k] = v

    def store_hdf5(self, file_name):
        """
        Store wavefront to HDF5 file (attributes and values).

        :param file_name: output HDF5 file name
        :type  file_name: string
        """
        utils.store_dict_hdf5(file_name, self._to_dict())
        self._store_attributes(file_name)

    def load_hdf5(self, file_name):
        """
        Load wavefront from HDF5 file.

        :param file_name: output HDF5 file name
        :type  file_name: string
        """
        self.custom_fields = {}  # fields of a previously loaded file must not stay
        self._update_from_dict(utils.load_dict_slash_hdf5(file_name))

    def get_intensity(self, slice_number=None, polarization=None):
        """
        Return intensity of wavefront

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of intensities
        """

        if polarization == 'total' or (polarization is None):
            pol = 6
        elif polarization == 'horizontal':
            pol = 0
        elif polarization == 'vertical':
            pol = 1
        else:
            raise ValueError(
                'unknown polarization value, should be "total" or "horizontal" or "vertical"')

        res = np.zeros(self._get_total_elements(), dtype='float32')

        if not res.flags['C_CONTIGUOUS']:
            res = np.ascontiguousarray(res)

        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 0, 6, self.params.photonEnergy, 0, 0)
        res = np.asarray(res, dtype='float32')
        res = res.reshape((
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices))
        if slice_number is not None:
            res = res[:, :, slice_number]
        return res

    def get_phase(self, slice_number=None, polarization='horizontal'):
        """
        Return phase of wavefront.

        :param polarization: 'horizontal' (default, also for None) or 'vertical'.
            'total' raises ValueError: the phase of the sum of two polarizations is not defined.
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of phases
        """
        _part_polarization(polarization)
        res = np.arctan2(self.get_imag_part(slice_number=slice_number,
                                            polarization=polarization),
                         self.get_real_part(slice_number=slice_number,
                                            polarization=polarization))
        return res

    def _get_field_part(self, int_type, slice_number, polarization):
        """Real (int_type=5) or imaginary (int_type=6) part of one polarization component."""
        pol = _part_polarization(polarization)
        res = np.zeros(self._get_total_elements(), dtype='float32')
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, int_type, 6, self.params.photonEnergy, 0, 0)
        res = np.asarray(res, dtype='float32')
        res = res.reshape((
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices))
        if slice_number is not None:
            res = res[:, :, slice_number]
        return res

    def get_real_part(self, slice_number=None, polarization='horizontal'):
        """
        Return real part of wavefront.

        :param polarization: 'horizontal' (default, also for None) or 'vertical'.
            'total' raises ValueError: SRW returns only the horizontal component for it.
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of real parts
        """
        return self._get_field_part(5, slice_number, polarization)

    def get_imag_part(self, slice_number=None, polarization='horizontal'):
        """
        Return imaginary part of wavefront.

        :param polarization: 'horizontal' (default, also for None) or 'vertical'.
            'total' raises ValueError: SRW returns only the horizontal component for it.
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of imaginary parts
        """
        return self._get_field_part(6, slice_number, polarization)

    def get_limits(self, axis='z'):
        """
        Get wavefront mesh limits [xmin, xmax, ....].

        Used in 2D visualization tools (as pylab.imshow(wfr_data, extends=wrf.get_limits()))

        :params axis: 'x','y' or 'z'
        :type axis: string

        :return: tuple of 4 floats: coordinates [m] in R-space, angles [rad] in Q-space
            (SRW spatial frequency [1/m] times the wavelength); time [s] or photon energy [eV]
            along the slices
        """
        sr = self.params.Mesh
        rep = self.params.wSpace
        if rep == 'R-space':
            if axis == 'z':
                return sr.xMin, sr.xMax, sr.yMax, sr.yMin
            elif axis == 'x':
                return sr.sliceMin, sr.sliceMax, sr.yMax, sr.yMin
            elif axis == 'y':
                return sr.sliceMin, sr.sliceMax, sr.xMax, sr.xMin
        elif rep == 'Q-space':
            wl = self.wavelength
            if axis == 'z':
                return sr.qxMin * wl, sr.qxMax * wl, sr.qyMax * wl, sr.qyMin * wl
            elif axis == 'x':
                return sr.sliceMin, sr.sliceMax, sr.qyMax * wl, sr.qyMin * wl
            elif axis == 'y':
                return sr.sliceMin, sr.sliceMax, sr.qxMax * wl, sr.qxMin * wl

    def __str__(self):
        """
        String representation to enable print function.

        :return: String representation
        """
        def node_str(node):
            return ''.join('{0} = {1}\n'.format(f.name, getattr(node, f.name))
                           for f in sorted(node.fields(), key=lambda f: f.name))

        mesh_str = 'Mesh:\n\t' + \
            '\n\t'.join(node_str(self.params.Mesh).split('\n')) + '\n'
        radiation_str = 'Radiation:\n\t' + \
            '\n\t'.join(node_str(self.params).split('\n')) + '\n'


        return radiation_str + mesh_str     

    def srw_info(self):
        """
        Print self._srwl_wf string representation. Used for debugging.

        :return: string
        """
        return srw_obj2str(self._srwl_wf)


    def _construct_properties(self):
        """
        Construct a dictionary of wavefront properties.
        """
        props = {
            'wavelength': self.wavelength
            }            # Add more properties as needed
             
        return props

    def get_property(self, property_name):
        """
        Retrieve a property by name.

        :param property_name: Name of the property to retrieve.
        :return: The property value.
        """
        value = self._construct_properties()
        for part in property_name.split('.'):
            value = value.get(part)
            if value is None:
                raise AttributeError(f"Property '{property_name}' not found.")
        return value
    
    def set_electric_field_representation(self, domain):
        """
        wrapper for srwlpy.SetRepresElecField
        
        sets the electric field representation
        
        :param domain: choice ofangular - 'a' <-> real-space 'c', frequency 'f' <--> time 't' [str]
        """
        srwlpy.SetRepresElecField(self._srwl_wf, domain)
        
        
    # axis -> (representation field, its value, code for SetRepresElecField, mesh coordinate,
    #          index of the axis in the intensity array (ny, nx, nSlices))
    _AXES = {
        'x': ('wSpace', 'R-space', 'c', 'x', 1),
        'y': ('wSpace', 'R-space', 'c', 'y', 0),
        'qx': ('wSpace', 'Q-space', 'a', 'x', 1),
        'qy': ('wSpace', 'Q-space', 'a', 'y', 0),
        't': ('wDomain', 'time', 't', 'slice', 2),
        'f': ('wDomain', 'frequency', 'f', 'slice', 2),
    }

    def _represented(self, axis):
        """
        Wavefront in the representation where the axis is defined.

        :param axis: 'x', 'y', 'qx', 'qy', 't' or 'f'
        :return: self if it is already in this representation, otherwise a converted copy
            (the field of self is not changed)
        """
        if axis not in self._AXES:
            raise ValueError(f"Unknown axis {axis!r}, should be one of {list(self._AXES)}")
        field, value, code = self._AXES[axis][:3]
        if getattr(self.params, field) == value:
            return self
        res = Wavefront(copy.deepcopy(self._srwl_wf))
        res.set_electric_field_representation(code)
        return res

    def _mesh_range(self, axis):
        """(min, max, number of points) of the axis; the wavefront must be in its representation.

        Q-space mesh of SRW is a spatial frequency [1/m], it is converted to angle [rad]
        by multiplication by the wavelength (as in get_limits).
        """
        mesh = self._srwl_wf.mesh
        coord = self._AXES[axis][3]
        if coord == 'x':
            lo, hi, n = mesh.xStart, mesh.xFin, mesh.nx
        elif coord == 'y':
            lo, hi, n = mesh.yStart, mesh.yFin, mesh.ny
        else:
            lo, hi, n = mesh.eStart, mesh.eFin, mesh.ne
        if axis in ('qx', 'qy'):
            lo, hi = lo * self.wavelength, hi * self.wavelength
        return lo, hi, int(n)

    def _line_profile(self, axis, method):
        """Profile of intensity along the axis; the wavefront must be in its representation."""
        intensity = self.get_intensity()
        index = self._AXES[axis][4]
        if method == 'center':
            item = [n // 2 for n in intensity.shape]
            item[index] = slice(None)
            return intensity[tuple(item)]
        other_axes = tuple(i for i in range(3) if i != index)
        if method == 'mean':
            return intensity.mean(axis=other_axes)
        return intensity.sum(axis=other_axes)

    def get_pixel_size(self, axis):
        """
        Calculate and return the pixel size (mesh step) along the specified axis.

        The step is (max - min)/(n - 1), as in SRW and in get_axis.
        If the wavefront is not in the representation of the axis (real or reciprocal
        space, time or frequency), the step is taken from a converted copy:
        the wavefront itself is not changed.

        :param axis: The axis for which to calculate the pixel size.
                     Options are 'qx' (angle x, [rad]), 'qy' (angle y, [rad]),
                     'x' (spatial x, [m]), 'y' (spatial y, [m]), 't' (time, [s]),
                     'f' (photon energy, [eV]). [str]
        :return: The pixel size along the specified axis. [float]
        :raises ValueError: If an unknown axis is specified or the mesh has one point along it.
        """
        lo, hi, n = self._represented(axis)._mesh_range(axis)
        if n < 2:
            raise ValueError(f"Pixel size along {axis!r} is not defined for a mesh of {n} point(s)")
        return (hi - lo) / (n - 1)

    def get_axis(self, axis):
        """
        Retrieve the values of the specified axis.

        If the wavefront is not in the representation of the axis (real or reciprocal
        space, time or frequency), the axis is taken from a converted copy:
        the wavefront itself is not changed.

        :param axis: The axis for which to retrieve the values.
                     Options are 'qx' (angle x, [rad]), 'qy' (angle y, [rad]),
                     'x' (spatial x, [m]), 'y' (spatial y, [m]), 't' (time, [s]),
                     'f' (photon energy, [eV]). [str]
        :return: A numpy array containing the values of the specified axis. [np.ndarray]
        :raises ValueError: If an unknown axis is specified.
        """
        return np.linspace(*self._represented(axis)._mesh_range(axis))

    def get_profile(self, axis, method='mean'):
        """
        Extract a one-dimensional line profile of intensity along the specified axis.

        The intensity is calculated in the representation of the axis (real or reciprocal
        space, time or frequency); if the wavefront is in another representation, a converted
        copy is used: the wavefront itself is not changed.

        :param axis: The axis along which to extract the profile.
                     Options are 'x' (spatial x), 'y' (spatial y), 't' (time),
                     'f' (frequency), 'qx' (angle x), 'qy' (angle y). [str]
        :param method: Method to use for extracting the line profile.
                       Options are 'mean' (average over two other dimensions),
                       'center' (line through the central point of two other dimensions), or
                       'sum' (sum over two other dimensions). [str]
                       Default is 'mean'.
        :return: One-dimensional numpy array representing the profile along the specified axis,
            its length is the number of points along the axis (see get_axis). [np.ndarray]

        :raises AssertionError: If an invalid axis or method is specified.
        """
        assert axis in ['x', 'y', 't', 'f', 'qx', 'qy'], \
            "Axis must be 'x', 'y', 't', 'f', 'qx', or 'qy'."
        assert method in ['mean', 'center', 'sum'], \
            "Method must be 'mean', 'center', or 'sum'."
        return self._represented(axis)._line_profile(axis, method)

    def get_beam_width(self, axis, method='mean'):
        """
        Calculate the gaussian width of the intensity distribution.

        This function computes the gaussian width of the beam profile along the specified axis by
        fitting a Gaussian function a*exp(-(x - x0)**2/(2*sigma**2)) to the profile (see get_profile).
        The fit is made on the axis normalized by the centre and the RMS width of the profile,
        starting from these moments. The wavefront itself is not changed.

        :param axis: The axis along which to calculate the beam width.
                     Options are 'x' (spatial x), 'y' (spatial y),
                     'qx' (angle x), 'qy' (angle y), 't' (time), 'f' (frequency). [str]
        :param method: Method used to extract the profile, can be 'mean', 'center', or 'sum'. [str, optional]
                       Default is 'mean'.
        :return: A tuple containing the parameters of the fitted Gaussian
                 in the units of get_axis (m, rad, s or eV):
                 - sigma (RMS width of intensity, positive) [float]
                 - x0 (mean position) [float]
                 - a (amplitude, units of intensity) [float]
        :rtype: tuple (float, float, float)

        :raises AssertionError: If an invalid axis is specified.
        :raises ValueError: If the profile has no positive intensity.
        """
        assert axis in ['x', 'y', 'qx', 'qy', 't', 'f'], \
            "Axis must be 'x', 'y', 'qx', 'qy', 't', or 'f'."

        wf = self._represented(axis)
        axis_values = np.linspace(*wf._mesh_range(axis))
        profile = np.asarray(wf._line_profile(axis, method), dtype='float64')

        weight = np.clip(profile, 0, None)
        if not weight.sum() > 0:
            raise ValueError(f"Profile along {axis!r} has no positive intensity")
        center = (axis_values * weight).sum() / weight.sum()
        rms = np.sqrt(((axis_values - center) ** 2 * weight).sum() / weight.sum())
        if not rms > 0:
            raise ValueError(f"Profile along {axis!r} has zero width")
        a_max = profile.max()

        # fit on the axis of order 1: curve_fit makes an absolute step ~1e-8 for a zero parameter
        params, _ = fit_gaussian((axis_values - center) / rms, profile / a_max, p0=[1., 0., 1.])

        sigma = float(abs(params[0]) * rms)
        x0 = float(center + params[1] * rms)
        a = float(params[2] * a_max)
        return sigma, x0, a

    def get_efield(self, polarization='horizontal'):
        """
        Convert electric field data to a complex representation.

        This function converts the electric field components from real and imaginary parts
        into complex numbers, based on the specified polarization. The available polarizations
        are 'horizontal', 'vertical', and 'circular'.

        - For 'horizontal' polarization, it uses `arrEhor`.
        - For 'vertical' polarization, it uses `arrEver`.
        - For 'circular' polarization, it combines both `arrEhor` and `arrEver`.

        The resulting complex array is suitable for further computations in wave optics, including outside of WPG.

        :param polarization: Specifies the polarization of the electric field.
                             Options are 'horizontal', 'vertical', or 'circular'. [str, optional]
                             Default is 'horizontal'.

        :return: A complex-valued numpy array representing the electric field.
                 - For 'horizontal' and 'vertical' polarizations: [3D array]
                 - For 'circular' polarization: [4D array, with the last dimension size 2]
        :rtype: numpy.ndarray

        :raises ValueError: If an unsupported polarization is specified.
        """

        if polarization == 'horizontal':
            arr = self.data.arrEhor
            efield = arr[:, :, :, 0].astype('complex128') + (arr[:, :, :, 1] * 1j)

        elif polarization == 'vertical':
            arr = self.data.arrEver
            efield = arr[:, :, :, 0].astype('complex128') + (arr[:, :, :, 1] * 1j)

        elif polarization == 'circular':
            efield = np.zeros_like(self.data.arrEhor, dtype='complex128')
            efield[:,:,:,0] = self.data.arrEhor[:,:,:,0].astype('complex128') + (self.data.arrEhor[:,:,:,1] * 1j)
            efield[:,:,:,1] = self.data.arrEver[:,:,:,0].astype('complex128') + (self.data.arrEver[:,:,:,1] * 1j)

        else:
            raise ValueError("Unsupported polarization. Choose 'horizontal', 'vertical', or 'circular'.")

        return efield

    @property
    def wavelength(self):
        r"""
        Calculate and return the wavelength based on photon energy.

        This property calculates the wavelength using the relationship:

        .. math::
            \lambda = \frac{hc}{E}

        where:
            - \( \lambda \) is the wavelength in meters.
            - \( h \) is Planck's constant (4.135667696 x 10^-15 eV·s).
            - \( c \) is the speed of light (299792458 m/s).
            - \( E \) is the photon energy in electron volts (eV).

        :return: The calculated wavelength in meters. [float]
        :raises ValueError: If the photon energy is not defined in `params`.
        """
        if self.params.photonEnergy:
            h = 4.135667696e-15  # Planck's constant in eV·s
            c = 299792458.0  # Speed of light in m/s
            return h * c / self.params.photonEnergy
        else:
            raise ValueError("Photon energy is not defined.")

    @property
    def divergence(self):
        """
        Calculate and return the beam divergence (RMS angular width of intensity).

        Gaussian widths (see get_beam_width) of the intensity in Q-space along qx and qy.
        The wavefront itself is not changed.

        :return: A tuple (sigma_qx, sigma_qy) of the divergence widths in radians. [(float, float)]
        """
        wf = self._represented('qx')
        return wf.get_beam_width(axis='qx')[0], wf.get_beam_width(axis='qy')[0]

    @property
    def pulse_duration(self):
        """
        Calculate and return the pulse duration.

        Gaussian width (see get_beam_width) of the intensity in the time domain.
        The wavefront itself is not changed.

        :return: RMS duration of the intensity (sigma of the gaussian) in seconds. [float]
        """
        return self.get_beam_width(axis='t')[0]


def _part_polarization(polarization):
    """SRW polarization code for the real and imaginary parts and the phase of the field.

    'total' is rejected: for it SRW returns the horizontal component only.
    """
    if polarization is None or polarization == 'horizontal':
        return 0
    elif polarization == 'vertical':
        return 1
    elif polarization == 'total':
        raise ValueError(
            "polarization 'total' is not defined for the real and imaginary parts and the phase "
            "of the field (SRW returns the horizontal component only), use 'horizontal' or 'vertical'")
    raise ValueError(
        'unknown polarization value, should be "horizontal" or "vertical"')
