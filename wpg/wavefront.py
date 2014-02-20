# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov'

import array
import warnings

import numpy
import h5py

import srwlib

import utils
import glossary

from utils import srw_obj2str

warnings.filterwarnings('ignore', category=Warning)

__author__ = 'makov'


class Wavefront(object):

    def __init__(self, srwl_wavefront=None):

        if srwl_wavefront is None:
            self._srwl_wf = srwlib.SRWLWfr()
        else:
            self._srwl_wf = srwl_wavefront

        self._wf_fields = {}
        self.custom_fields = {}

        for wf_field in glossary.get_wf_fields():
            wf = wf_field(self)
            self.add_field(wf)

    def get_total_elements(self):
        return self.params.Mesh.nx * self.params.Mesh.ny * self.params.Mesh.nSlices

    def allocate_moments(self):
        self._srwl_wf.arMomX = array.array(
            'd', [0] * self.params.Mesh.nSlices * 11)
        self._srwl_wf.arMomY = array.array(
            'd', [0] * self.params.Mesh.nSlices * 11)

    def add_field(self, wf_field):
        """ Add field to wavefront and create field"""

        class glossary_folder(object):

            """Glossary folder. Empty class."""

            pass

        def get_value(self):
            return wf_field.value

        def set_value(self, value):
            wf_field.value = value

        def get_doc():
            return wf_field.value.__doc__

        if not isinstance(wf_field, glossary.RadiationField):
            raise TypeError('wf_field must be RadiationField')

        self._wf_fields[wf_field.glossary_name] = wf_field

        node = self
        keys_chain = wf_field.keys_chain

        for key in keys_chain[:-1]:
            if not key in node.__dict__:
                node.__dict__[key] = glossary_folder()
            node = node.__dict__[key]

        setattr(node.__class__, keys_chain[-1], property(get_value,
                                                         set_value, doc=get_doc()))

    def to_dict(self):
        res = {}
        for (key, value) in self._wf_fields.iteritems():
            res[key] = value.value
        
        res.update(self.custom_fields)
        return res

    def update_from_dict(self, in_dict):
        for (key, value) in in_dict.iteritems():
            if key in self._wf_fields:
                self._wf_fields[key].value = value
            else:
                utils.update_dict_slash_string(self.custom_fields, key, value)

    def store_attributes(self, file_name):
        with h5py.File(file_name) as h5f:
            for (key, wff) in self._wf_fields.iteritems():
                try:
                    if wff.glossary_name in h5f:
                        for (k, v) in wff.attributes.items():
                            h5f[wff.glossary_name].attrs[k] = v
                except KeyError:
                    pass

    def store_hdf5(self, file_name):
        utils.store_dict_hdf5(file_name, self.to_dict())
        self.store_attributes(file_name)

    def load_hdf5(self, file_name):
        self.update_from_dict(utils.load_dict_slash_hdf5(file_name))

    def get_intensity(self, slice_number=None, polarization=None):
        """
        Return intensity of wavefront

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :return:
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

        res = array.array('f', [0] * self.get_total_elements())
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 0, 6, self.params.photonEnergy, 0, 0)
        res = numpy.array(res, dtype='float32')
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if not slice_number is None:
            res = res[:, :, slice_number]
        return res

    def get_phase(self, slice_number=None, polarization=None):
        """
        Return phase of wavefront

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :return:
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

        res = numpy.arctan2(self.get_imag_part(slice_number=slice_number,
                                               polarization=polarization),
                            self.get_real_part(slice_number=slice_number,
                                               polarization=polarization))

        # res = array.array('f',[0]*self.get_total_elements())
        # res = srwlib.srwl.CalcIntFromElecField(res, self._srwl_wf, pol, 0, 6, self.params.photonEnergy, 0, 0.)
        # res = numpy.array(res, dtype='float32')
        # res.shape = (self.params.Mesh.ny,self.params.Mesh.nx,self.params.Mesh.nSlices)
        # if not slice is None:
        #     res = res[:, :, slice]

        return res

    def get_real_part(self, slice_number=None, polarization=None):
        """
        Return real part of wavefront

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :param slice_number: slice number ti return, if None - get 3D array (all slices)
        :return:
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

        res = array.array('f', [0] * self.get_total_elements())
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 5, 6, self.params.photonEnergy, 0, 0)
        res = numpy.array(res, dtype='float32')
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if not slice_number is None:
            res = res[:, :, slice_number]
        return res

    def get_imag_part(self, slice_number=None, polarization=None):
        """
        Return imaginary part of wavefront

        :param polarization: 'total' or 'horizontal' or 'vertical'
        :param slice: slice number ti return, if None - get 3D array (all slices)
        :return:
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

        res = array.array('f', [0] * self.get_total_elements())
        res = srwlib.srwl.CalcIntFromElecField(
            res, self._srwl_wf, pol, 6, 6, self.params.photonEnergy, 0, 0)
        res = numpy.array(res, dtype='float32')
        res.shape = (
            self.params.Mesh.ny, self.params.Mesh.nx, self.params.Mesh.nSlices)
        if not slice_number is None:
            res = res[:, :, slice_number]
        return res

    def get_limits(self, axis='z'):
        sr = self.params.Mesh
        if axis == 'z':
            return sr.xMin, sr.xMax, sr.yMax, sr.yMin
        elif axis == 'x':
            return sr.sliceMin, sr.sliceMax, sr.yMax, sr.yMin
        elif axis == 'y':
            return sr.sliceMin, sr.sliceMax, sr.xMax, sr.xMin

    def __str__(self):
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
        return srw_obj2str(self._srwl_wf)
