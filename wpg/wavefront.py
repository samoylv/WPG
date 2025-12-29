# -*- coding: utf-8 -*-
"""
This module contains base wrapper for SRWLWfr (Wavefront). It's implement numpy inter operations to SRWLWfr structure, serialization to HDF5, visualization tools, etc.

.. module:: wpg.wavefront
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

import array
import warnings

import numpy as np
import h5py

from srwpy import srwlpy 
import srwpy.srwlib as srwlib

import wpg.utils as utils
import wpg.glossary as glossary

from wpg.utils import srw_obj2str

from wpg.wpg_uti_math import fit_gaussian

warnings.filterwarnings('ignore', category=Warning)


class Wavefront(object):
    """
    This is base class for manipulation with wavefronts in wpg module.

    One of most important field is _srwl_wf (instance of srwlib.SRWLWfr). SEtting and getting this field allows to call all SRWLpy functions.
    """

    def __init__(self, srwl_wavefront=None):
        """
        Create wavefront instance.

        The most important wavefront fields dynamically initialize from :mod:`wpg.glossry`

        :param srwl_wavefront: if present, wavefront inits with it's parameters  
        :type srwl_wavefront: srwlib.SRWLWfr
        :return: Wavefront instance.
        """
        if srwl_wavefront is None:
            self._srwl_wf = srwlib.SRWLWfr()
        else:
            self._srwl_wf = srwl_wavefront

        self._wf_fields = {}
        self.custom_fields = {}

        for wf_field in glossary.get_wf_fields():
            wf = wf_field(self)
            self._add_field(wf)

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

    def _add_field(self, wf_field):
        """
        Add field to wavefront structure and create field.

        :param wf_field: field instance
        :type wf_field: wpg.glossary.RadiationField
        """

        class glossary_folder(object):
            """Glossary folder. Empty class to build dictionary tree."""
            pass

        def get_value(self):
            """Get value stored in field."""
            return wf_field.value

        def set_value(self, value):
            """
            Get value stored in field.

            :param value: value to be stored
            """
            wf_field.value = value

        def get_doc():
            """Get field documentation string."""
            return wf_field.value.__doc__

        if not isinstance(wf_field, glossary.RadiationField):
            raise TypeError('wf_field must be RadiationField')

        self._wf_fields[wf_field.glossary_name] = wf_field

        node = self
        keys_chain = wf_field.keys_chain

        for key in keys_chain[:-1]:
            if key not in node.__dict__:
                node.__dict__[key] = glossary_folder()
            node = node.__dict__[key]

        setattr(node.__class__, keys_chain[-1], property(get_value,
                                                         set_value, doc=get_doc()))

    def _to_dict(self):
        """
        Convert wavefront to dictionary. Used for saving in HDF5 file.

        :return: dictionary view of wavefront
        """
        res = {}
        for (key, value) in self._wf_fields.items():
            res[key] = value.value

        res.update(self.custom_fields)
        return res

    def _update_from_dict(self, in_dict):
        """
        Update wavefront from dictionary. Used for loading wavefront to HDF5 file.

        :param in_dict: input dictionary
        :type in_dict: dict
        """
        for (key, value) in in_dict.items():
            # python3 hack
            if isinstance(key, bytes):
                key = key.decode('utf-8')

            if key in self._wf_fields:
                self._wf_fields[key].value = value
            else:
                utils.update_dict_slash_string(self.custom_fields, key, value)

    def _store_attributes(self, file_name):
        """
        Store wavefront attributes to HDF5 file.

        Attribute of each field is values of field.attributes

        :param file_name: output HDF5 file name
        :type  file_name: string
        """
        with h5py.File(file_name, 'r+') as h5f:
            for (key, wff) in self._wf_fields.items():
                try:
                    if wff.glossary_name in h5f:
                        for (k, v) in list(wff.attributes.items()):
                            h5f[wff.glossary_name].attrs[k] = v
                except KeyError:
                    pass

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
        res = np.array(res, dtype='float32', copy=False)
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if slice_number is not None:
            res = res[:, :, slice_number]
        return res

    def get_phase(self, slice_number=None, polarization=None):
        """
        Return phase of wavefront.

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of phases
        """
        # TODO: bug with freeze

        if polarization == 'total' or (polarization is None):
            pol = 6
            print(
                'Attention!!! The "total" polarization behavior sometimes strange. Use "horizontal" or "vertical".'
            )
        elif polarization == 'horizontal':
            pol = 0
        elif polarization == 'vertical':
            pol = 1
        else:
            raise ValueError(
                'unknown polarization value, should be "total" or "horizontal" or "vertical"')

        res = np.arctan2(self.get_imag_part(slice_number=slice_number,
                                            polarization=polarization),
                         self.get_real_part(slice_number=slice_number,
                                            polarization=polarization))

        # res = array.array('f',[0]*self.get_total_elements())
        # res = srwlib.srwl.CalcIntFromElecField(res, self._srwl_wf, pol, 0, 6, self.params.photonEnergy, 0, 0.)
        # res = np.array(res, dtype='float32')
        # res.shape = (self.params.Mesh.ny,self.params.Mesh.nx,self.params.Mesh.nSlices)
        # if not slice is None:
        #     res = res[:, :, slice]

        return res

    def get_real_part(self, slice_number=None, polarization=None):
        """
        Return real part of wavefront.

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of real parts
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
            res, self._srwl_wf, pol, 5, 6, self.params.photonEnergy, 0, 0)
        res = np.array(res, dtype='float32', copy=False)
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if slice_number is not None:
            res = res[:, :, slice_number]
        return res

    def get_imag_part(self, slice_number=None, polarization=None):
        """
        Return imaginary part of wavefront.

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :type polarization: string
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :type slice_number: int or range
        :return: array of imaginary parts
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
            res, self._srwl_wf, pol, 6, 6, self.params.photonEnergy, 0, 0)
        res = np.array(res, dtype='float32', copy=False)
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if slice_number is not None:
            res = res[:, :, slice_number]
        return res

    def get_limits(self, axis='z'):
        """
        Get wavefront mesh limits [xmin, xmax, ....].

        Used in 2D visualization tools (as pylab.imshow(wfr_data, extends=wrf.get_limits()))

        :params axis: 'x','y' or 'z'
        :type axis: string

        :return: list of integers
        """
        sr = self.params.Mesh
        rep = self.params.wSpace
        if rep == 'R-space':
            print(rep)
            if axis == 'z':
                return sr.xMin, sr.xMax, sr.yMax, sr.yMin
            elif axis == 'x':
                return sr.sliceMin, sr.sliceMax, sr.yMax, sr.yMin
            elif axis == 'y':
                return sr.sliceMin, sr.sliceMax, sr.xMax, sr.xMin
        elif rep == 'Q-space':
            print(rep)
            wl = 12.398 * 1e-10 / (self.params.photonEnergy * 1e-3)  # WaveLength
            # wv = 2.*np.pi/wl
            # #WaveVector
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
        mesh_str = 'Mesh:\n\t' + \
            '\n\t'.join(srw_obj2str(self.params.Mesh).split('\n')) + '\n'
        radiation_str = 'Radiation:\n\t' + \
            '\n\t'.join(srw_obj2str(self.params).split('\n')) + '\n'


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
        value = self.properties
        for part in property_name.split('.'):
            value = value.get(part)
            if value is None:
                raise AttributeError(f"Property '{property_name}' not found.")
        return value
    
    def __getattr__(self, name):
        """
        Map higher level attributes to their nested counterparts
        """
        if hasattr(self.params, name):
            return getattr(self.params, name)
        elif hasattr(self.params.Mesh, name):
            return getattr(self.params.Mesh, name)
        raise AttributeError(f"'{self.__class__.__name__}' object has no attribute '{name}'")
        
    def set_electric_field_representation(self, domain):
        """
        wrapper for srwlpy.SetRepresElecField
        
        sets the electric field representation
        
        :param domain: choice ofangular - 'a' <-> real-space 'c', frequency 'f' <--> time 't' [str]
        """
        srwlpy.SetRepresElecField(self._srwl_wf, domain)
        
        
    def get_pixel_size(self, axis):
        """
        Calculate and return the pixel size along the specified axis.
        
        This function computes the pixel size for the given axis based on the 
        mesh grid parameters. It adjusts the electric field representation to 
        the appropriate domain (real or reciprocal space, time, or frequency) 
        and then calculates the pixel size.
    
        :param axis: The axis for which to calculate the pixel size. 
                     Options are 'qx' (reciprocal x), 'qy' (reciprocal y), 
                     'x' (spatial x), 'y' (spatial y), 't' (time), 'f' (frequency). [str]
        :return: The pixel size along the specified axis. [float]
        :raises ValueError: If an unknown axis is specified.
        """
        current_wSpace = self.params.wSpace
        current_wDomain = self.params.wDomain
    
        if axis == 'qx':
            self.set_electric_field_representation('a')
            pixel_size = (self.params.Mesh.qxMax - self.params.Mesh.qxMin) / self.params.Mesh.nx
        elif axis == 'qy':
            self.set_electric_field_representation('a')
            pixel_size = (self.params.Mesh.qyMax - self.params.Mesh.qyMin) / self.params.Mesh.ny
        elif axis == 'x':
            self.set_electric_field_representation('c')
            pixel_size = (self.params.Mesh.xMax - self.params.Mesh.xMin) / self.params.Mesh.nx
        elif axis == 'y':
            self.set_electric_field_representation('c')
            pixel_size = (self.params.Mesh.yMax - self.params.Mesh.yMin) / self.params.Mesh.ny
        elif axis == 't':
            self.set_electric_field_representation('t')
            pixel_size = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
        elif axis == 'f':
            self.set_electric_field_representation('f')
            pixel_size = (self.params.Mesh.sliceMax - self.params.Mesh.sliceMin) / self.params.Mesh.nSlices
        else:
            raise ValueError(f"Unknown axis '{axis}'")
    
        # Revert to original representation
        if current_wSpace in ['R-space', 'Q-space']:
            self.set_electric_field_representation('a' if current_wSpace == 'Q-space' else 'c')
        elif current_wDomain in ['time', 'frequency']:
            self.set_electric_field_representation('t' if current_wDomain == 'time' else 'f')
    
        return pixel_size
    
    def get_axis(self, axis):
        """
        Retrieve the values of the specified axis.
        
        This function generates an array of values for the given axis based on 
        the mesh grid parameters. It adjusts the electric field representation 
        to the appropriate domain (real or reciprocal space, time, or frequency) 
        and then constructs the axis values.
    
        :param axis: The axis for which to retrieve the values. 
                     Options are 'qx' (reciprocal x), 'qy' (reciprocal y), 
                     'x' (spatial x), 'y' (spatial y), 't' (time), 'f' (frequency). [str]
        :return: A numpy array containing the values of the specified axis. [np.ndarray]
        :raises ValueError: If an unknown axis is specified.
        """
        current_wSpace = self.params.wSpace
        current_wDomain = self.params.wDomain
    
        if axis == 'qx':
            self.set_electric_field_representation('a')
            axis_values = np.linspace(self.params.Mesh.qxMin, self.params.Mesh.qxMax, self.params.Mesh.nx)
        elif axis == 'qy':
            self.set_electric_field_representation('a')
            axis_values = np.linspace(self.params.Mesh.qyMin, self.params.Mesh.qyMax, self.params.Mesh.ny)
        elif axis == 'x':
            self.set_electric_field_representation('c')
            axis_values = np.linspace(self.params.Mesh.xMin, self.params.Mesh.xMax, self.params.Mesh.nx)
        elif axis == 'y':
            self.set_electric_field_representation('c')
            axis_values = np.linspace(self.params.Mesh.yMin, self.params.Mesh.yMax, self.params.Mesh.ny)
        elif axis == 't':
            self.set_electric_field_representation('t')
            axis_values = np.linspace(self.params.Mesh.sliceMin, self.params.Mesh.sliceMax, self.params.Mesh.nSlices)
        elif axis == 'f':
            self.set_electric_field_representation('f')
            axis_values = np.linspace(self.params.Mesh.sliceMin, self.params.Mesh.sliceMax, self.params.Mesh.nSlices)
        else:
            raise ValueError(f"Unknown axis '{axis}'")
    
        # Revert to original representation
        if current_wSpace in ['R-space', 'Q-space']:
            self.set_electric_field_representation('a' if current_wSpace == 'Q-space' else 'c')
        elif current_wDomain in ['time', 'frequency']:
            self.set_electric_field_representation('t' if current_wDomain == 'time' else 'f')
    
        return axis_values


    def get_profile(self, axis, method='mean'):
        """ 
        Extract a one-dimensional line profile along the specified axis.
        
        This function extracts a line profile from the intensity distribution
        of the wavefront along the specified axis using the specified method.
        The profile provides a representation of the intensity along a single dimension,
        useful for analyzing the spatial or temporal characteristics of the beam.
    
        :param axis: The axis along which to extract the profile. 
                     Options are 'x' (spatial x), 'y' (spatial y), 't' (time), 
                     'f' (frequency), 'qx' (reciprocal x), 'qy' (reciprocal y). [str]
        :param method: Method to use for extracting the line profile. 
                       Options are 'mean' (average over orthogonal dimensions), 
                       'center' (central slice), or 'sum' (sum over orthogonal dimensions). [str]
                       Default is 'mean'.
        :return: One-dimensional numpy array representing the profile along the specified axis. [np.ndarray]
        
        :raises AssertionError: If an invalid axis or method is specified.
        """
        assert axis in ['x', 'y', 't', 'f', 'qx', 'qy'], \
            "Axis must be 'x', 'y', 't', 'f', 'qx', or 'qy'."
        assert method in ['mean', 'center', 'sum'], \
            "Method must be 'mean', 'center', or 'sum'."
        
        intensity = self.get_intensity()
        axis_dict = {
            'x': (0, 'R-space', 'c', self.nx // 2, (0, 2)),
            'y': (1, 'R-space', 'c', self.ny // 2, (1, 2)),
            't': (2, 'time', 't', self.nSlices // 2, (0, 1)),
            'f': (2, 'frequency', 'f', self.nSlices // 2, (0, 1)),
            'qx': (0, 'Q-space', 'a', self.nx // 2, (1, 2)),
            'qy': (1, 'Q-space', 'a', self.ny // 2, (0, 2))
        }
        
        index, space, representation, center, mean_axes = axis_dict[axis]
    
        current_wSpace = self.params.wSpace
        current_wDomain = self.params.wDomain
        
        # Set to the desired representation
        self.set_electric_field_representation(representation)
    
        if method == 'mean':
            profile = intensity.mean(axis=mean_axes)
        elif method == 'center':
            profile = np.take(intensity, center, axis=index)
        elif method == 'sum':
            profile = intensity.sum(axis=mean_axes)
    
        # Revert back to original representation
        if current_wSpace in ['R-space', 'Q-space']:
            self.set_electric_field_representation('a' if current_wSpace == 'Q-space' else 'c')
        elif current_wDomain in ['time', 'frequency']:
            self.set_electric_field_representation('t' if current_wDomain == 'time' else 'f')
        
        return profile
    
    
    def get_beam_width(self, axis, method='mean'):
        """
        Calculate the gaussian width of the intensity distribution.
        
        This function computes the gaussian width of the beam profile along the specified axis by 
        fitting a Gaussian function to the extracted profile.
        
        :param axis: The axis along which to calculate the beam width. 
                     Options are 'x' (spatial x), 'y' (spatial y), 
                     'qx' (reciprocal x), 'qy' (reciprocal y), 't' (time), 'f' (frequency). [str]
        :param method: Method used to extract the profile, can be 'mean', 'center', or 'sum'. [str, optional]
                       Default is 'mean'.
        :return: A tuple containing the parameters of the fitted Gaussian:
                 - sigma (standard deviation) [float]
                 - x0 (mean position) [float]
                 - a (amplitude) [float]
        :rtype: tuple (float, float, float)
        
        :raises AssertionError: If an invalid axis is specified.
        """
        assert axis in ['x', 'y', 'qx', 'qy', 't', 'f'], \
            "Axis must be 'x', 'y', 'qx', 'qy', 't', or 'f'."
        
        # Get the appropriate axis values
        axis_values = self.get_axis(axis)
        
        # Extract the profile along the specified axis
        profile = self.get_profile(axis, method=method)
        
        # Initial guesses for Gaussian fitting
        sigma_guess = (axis_values.max() - axis_values.min()) / 2
        x0_guess = axis_values[np.argmax(profile)]
        a_guess = profile.max()
        
        initial_guesses = [sigma_guess, x0_guess, a_guess]
    
        # Fit the profile to a Gaussian function
        params, _ = fit_gaussian(axis_values, profile, p0=initial_guesses)
        
        sigma = params[0]
        x0 = params[1]
        a = params[2]
        
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
            efield = np.zeros_like(self.data.arrEhor[:,:,:,0], dtype='complex128')
            efield[:,:,:] = self.data.arrEhor[:,:,:,0].astype('complex128') + (self.data.arrEhor[:,:,:,1] * 1j)
        
        elif polarization == 'vertical':
            efield = np.zeros_like(self.data.arrEver[:,:,:,:-2], dtype='complex128')
            efield[:,:,:] = self.data.arrEver[:,:,:,0].astype('complex128') + (self.data.arrEver[:,:,:,1] * 1j)
        
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
            - \( c \) is the speed of light (3.0 x 10^8 m/s).
            - \( E \) is the photon energy in electron volts (eV).
    
        :return: The calculated wavelength in meters. [float]
        :raises ValueError: If the photon energy is not defined in `params`.
        """
        if self.params.get('photonEnergy'):
            h = 4.135667696e-15  # Planck's constant in eV·s
            c = 3.0e8  # Speed of light in m/s
            return h * c / self.params.photonEnergy
        else:
            raise ValueError("Photon energy is not defined.")
    
    @property
    def divergence(self):
        """ 
        Calculate and return the beam divergence width in reciprocal space.
        
        This property computes the divergence width for both the qx and qy axes, 
        typically in reciprocal space. Divergence is a measure of the spread of the 
        beam and is expressed in radians.
    
        :return: A tuple containing the divergence width in radians for the qx and qy axes. [(float, float)]
        """
        return self.get_beam_width(axis='qx'), self.get_beam_width(axis='qy')
    
    @property
    def pulse_duration(self):
        """
        Calculate and return the pulse duration.
        
        This property computes the gaussian width of the pulse 
        duration in the time domain. The pulse duration is an important parameter 
        characterizing the temporal width of the pulse.
    
        :return: The gaussian width of the pulse duration in seconds. [float]
        """
        return self.get_beam_width(axis='t')

