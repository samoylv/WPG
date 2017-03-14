from __future__ import print_function

__author__=’Hector Mauricio Castaneda Cortes’

import scipy.constants as const
import numpy as np
import h5py
from wpg import Wavefront


def open_slice_data(slice_dict, slice_index):
    '''
    Inputs:
      slice_dict ( h5 Object)
      slice_index (Slice Index)

    Output:
      slice_dict[slice_folder].value
      Contents of a particular entry of the Glossary (corresponding to the slice_index supplied as input)
    '''
    sl = [slice_index]
    if (0 < sl[0] < 10):
        sl_s = '00' + str(sl[0])
    elif (10 <= sl[0] <= 99):
        sl_s = '0' + str(sl[0])
    elif (sl[0] >= 100):
        sl_s = str(sl[0])
    slice_folder = '/slice000' + sl_s
    slice_folder = slice_folder + '/field'
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
    vac_imp = const.codata.value('characteristic impedance of vacuum')
    eev = 1e6 * const.codata.value('electron mass energy equivalent in MeV')

    # Definition of internal variables
    h5f = _hf
    npt = _npoints
    nsl = _nslices
    mesh_size = _grid_size / (npt - 1)
    lmb = _wv
    lmb_u = _lambda_un
    matrix_t = np.zeros(shape=(npt, npt, nsl, 2))

    # Definition of parameters needed for the scaling factor to \sqrt{W} units
    xkw0 = 2. * np.pi / lmb_u
    xks = 2. * np.pi / lmb
    dxy = xkw0 * mesh_size

    # Scaling factor in order to get the field in units of \sqrt{W}
    # (consistent with GENESIS v2)
    fact = dxy * eev * xkw0 / xks / np.sqrt(vac_imp)
    fact = fact / (xkw0 * xkw0)
    print(fact)

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
                        vect[(jc * npt) + lc] / (1000. * mesh_size)
                    ind = ind + 1
    return matrix_t

def read_genesis_file(gen_fname):
    speed_of_light = const.codata.value('speed of light in vacuum')
    h_eV_s = const.codata.value('Planck constant in eV s')
    lmb_und = 2.75e-2

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

    ### Definition of the Electric field field arrays where the electric field from the GENESIS output file
    ###  will be copied 

        wf.data.arrEhor = np.zeros(shape=(npoints, npoints, slice_count, 2))
        wf.data.arrEver = np.zeros(shape=(npoints, npoints, slice_count, 2))

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
        range_xy = (npoints - 1) * grid_size
        wf.params.Mesh.xMin = -range_xy / 2.
        wf.params.Mesh.xMax = range_xy / 2.
        wf.params.Mesh.yMin = -range_xy / 2.
        wf.params.Mesh.yMax = range_xy / 2.

    ### Extract the field data from the h5 file and fill in the data of the Electric field, by calling the
    ### vector_grid_conversion function   
        wf.data.arrEhor = vector_grid_conversion(
            hf, npoints, slice_count, range_xy, wavelength, lmb_und)

    return wf
