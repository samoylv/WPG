# -*- coding: utf-8 -*-
"""
This module contains base wrapper for SRWLWfr (Wavefront). It's implement numpy inter operations to SRWLWfr structure, serialization to HDF5, visualization tools, etc.

.. module:: wpg.wavefront
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import array
import warnings

import numpy as np
import h5py

import wpg.srwlib as srwlib

import wpg.utils as utils
import wpg.glossary as glossary

from wpg.utils import srw_obj2str

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
        with h5py.File(file_name) as h5f:
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

        res = array.array(str(u'f'), np.zeros(self._get_total_elements(), 'float32').tobytes())
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

        res = array.array(str(u'f'), np.zeros(self._get_total_elements(), 'float32').tobytes())
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 5, 6, self.params.photonEnergy, 0, 0)
        res = np.array(res, dtype='float32')
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

        res = array.array(str(u'f'), np.zeros(self._get_total_elements(), 'float32').tobytes())
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 6, 6, self.params.photonEnergy, 0, 0)
        res = np.array(res, dtype='float32')
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
            wl = 12.39 * 1e-10 / (self.params.photonEnergy * 1e-3)  # WaveLength
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
        mesh_str = 'Mesh:\n\t\t' + \
            '\n\t\t'.join(srw_obj2str(self.params.Mesh).split('\n')) + '\n\t'
        radiation_str = mesh_str + \
            '\n\t'.join(srw_obj2str(self.params).split('\n'))

        radiation_str = 'Radiation:\n\t' + radiation_str + '\n'

        data_ehor = '\tarrEhor = array of shape ' + \
            str(self.data.arrEhor.shape) + \
            ' // the 2-nd dimension is (re,im)\n'
        data_ever = '\tarrEver = array of shape ' + \
            str(self.data.arrEver.shape) + \
            ' // the 2-nd dimension is (re,im)\n'

        data_str = data_ehor + data_ever

        return radiation_str + data_str

    def srw_info(self):
        """
        Print self._srwl_wf string representation. Used for debugging.

        :return: string
        """
        return srw_obj2str(self._srwl_wf)
