#encoding: utf-8
from __future__ import print_function

__author__ = 'Hector Mauricio Castaneda Cortes'

import scipy.constants as const
import numpy as np
import h5py
from wpg import Wavefront
from copy import deepcopy

def open_slice_data(slice_dict, slice_index):
    '''
    Inputs:
      slice_dict ( h5 Object)
      slice_index (Slice Index)

    Output:
      slice_dict[slice_folder].value
      Contents of a particular entry of the Glossary (corresponding to the slice_index supplied as input)
    '''
    if (0 < slice_index < 10):
        sl_s = '00000' + str(slice_index)
    elif (10 <= slice_index <= 99):
        sl_s = '0000'+str(slice_index)
    elif(100<=slice_index< 1000):
        sl_s = '000'+str(slice_index)
    elif(1000<=slice_index<10000):
        sl_s = '00'+str(slice_index)
    elif(10000<=slice_index<100000):
        sl_s = '0'+str(slice_index)
    else:
        sl_s = str(slice_index)

    slice_folder = '/slice'+sl_s+'/field'
    return slice_dict[slice_folder].value


def vector_grid_conversion(_hf, _npoints, _nslices, _grid_size, _wv, _lambda_un):
    '''
    Inputs:
          _hf ( h5 Object)
          _npoints (Number of points in the grid, equivalent to ncar in GENESIS)
          _nslices (Number of slices)
          _grid_size (grid size extracted from h5)
          _wv radiation wavelength
          _lambda_un Undulator period
    Output:
          matrix_t   Matrix that contains the Electric field data along the grid
          (Real and imaginary parts)shape = (npoints,npoints, slice_count,2)))
    '''
    # Definition of constants
    e_charge = const.codata.value('elementary charge')
    vac_imp = const.codata.value('characteristic impedance of vacuum')
    eev = 1e6 * const.codata.value('electron mass energy equivalent in MeV')

    # Definition of internal variables
    h5f = _hf
    npt = _npoints
    nsl = _nslices
    mesh_size = 2.*_grid_size / (npt - 1)
    lmb = _wv
    lmb_u = _lambda_un
    matrix_t = np.zeros(shape=(npt, npt, nsl, 2))

    # Definition of parameters needed for the scaling factor to \sqrt{W} units
    xkw0 = 2. * np.pi / lmb_u
    xks = 2. * np.pi / lmb
    dxy =  xkw0 * mesh_size

    # Scaling factor in order to get the field in units of \sqrt{W}
    # (consistent with GENESIS v2)
    fact0 = dxy * eev * xkw0 / xks / np.sqrt(vac_imp)
    #fact = fact/ (xkw0 * xkw0)

    fact=float(np.sqrt(e_charge)*1e-2*fact0/(mesh_size))

    # Cycle over all slices, converting the 1D vector into a 3D Matrix and
    # dividing by the mesh size in order to get the field in units of \sqrt{W/mm^2}
    for islice in range(1, nsl):
        print('slice No ' + str(islice))
        tmp_data = open_slice_data(h5f, islice)
        for jc in range(npt):
            for lc in range(npt):
                ind = 0
                while ind < 2:
                    if ind == 0:
                        vect = tmp_data[0::2]
                    elif ind == 1:
                        vect = tmp_data[1::2]
                    matrix_t[jc, lc, islice - 1, ind] = fact * \
                        vect[(jc * npt) + lc]
                    ind = ind + 1
    return matrix_t


def read_genesis_file(gen_fname,lmb_und,pol='lin'):
    speed_of_light = const.codata.value('speed of light in vacuum')
    h_eV_s = const.codata.value('Planck constant in eV s')

    wf = Wavefront()

    ### Open the hdf5 output GENESIS file

    with h5py.File(gen_fname, 'r') as hf:
        grid_size = hf['gridsize'].value
        slice_count = hf['slicecount'].value
        print(slice_count)
        slice_spacing = hf['slicespacing'].value
        wavelength = hf['wavelength'].value
        data1 = open_slice_data(hf, 1)
        npoints = int(np.sqrt(len(data1) / 2))
        range_xy = (npoints - 1) * grid_size
         ### Fill in the fields of the wavefront object

        wf.params.wEFieldUnit = 'sqrt(W/mm^2)'
        wf.params.photonEnergy = h_eV_s * speed_of_light / wavelength
        print('Photon energy: ', wf.params.photonEnergy, 'eV')
        wf.params.wDomain = 'time'
        wf.params.Mesh.nSlices = slice_count
        wf.params.Mesh.nx = npoints
        wf.params.Mesh.ny = npoints
        pulse_length = (slice_count - 1) * slice_spacing / (speed_of_light)
        print('Pulse length: ', pulse_length)
        wf.params.Mesh.sliceMin = -pulse_length / 2.
        wf.params.Mesh.sliceMax = pulse_length / 2.

        wf.params.Mesh.xMin = -range_xy / 2.
        wf.params.Mesh.xMax = range_xy / 2.
        wf.params.Mesh.yMin = -range_xy / 2.
        wf.params.Mesh.yMax = range_xy / 2.
        wf.params.nval=2
        ### Extract the field data from the h5 file and fill in the data of the Electric field, by calling the
        ### vector_grid_conversion function

        if pol=='lin':
            wf.data.arrEhor = vector_grid_conversion(
                hf, npoints, slice_count, range_xy, wavelength, lmb_und)
            wf.data.arrEver = np.zeros_like(wf.data.arrEhor)
        else:
            matrix_E = vector_grid_conversion(hf, npoints, \
            slice_count, range_xy, wavelength, lmb_und)
            Re_Etotal = deepcopy(matrix_E[:,:,:,0])
            Im_Etotal = deepcopy(matrix_E[:,:,:,1])

            wf.data.arrEhor= deepcopy(matrix_E)
            wf.data.arrEver = deepcopy(wf.data.arrEhor)

            wf.data.arrEhor = np.multiply(np.sqrt(0.5),wf.data.arrEhor)
            wf.data.arrEver[:,:,:,0]= deepcopy(np.multiply(-1.0,wf.data.arrEhor[:,:,:,1]))
            wf.data.arrEver[:,:,:,1]= deepcopy(wf.data.arrEhor[:,:,:,0])
    return wf
