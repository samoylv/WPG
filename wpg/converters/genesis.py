#encoding: utf-8

__author__ = 'Hector Mauricio Castaneda Cortes'

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
      slice_dict[slice_folder][()]
      Contents of a particular entry of the Glossary (corresponding to the slice_index supplied as input)
    '''
    # GENESIS numbers the slice groups from 1 with six digits: /slice000001 ... /slice001000 ...
    if slice_index < 1:
        raise ValueError('GENESIS slice numbering starts at 1, got {}'.format(slice_index))
    slice_folder = '/slice%06d/field' % slice_index
    return slice_dict[slice_folder][()]


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
    vac_imp = const.value('characteristic impedance of vacuum')
    eev = 1e6 * const.value('electron mass energy equivalent in MeV')

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
    # GENESIS slices are numbered 1..nsl, slice k goes to the array index k - 1
    for islice in range(1, nsl + 1):
        print('slice No ' + str(islice))
        tmp_data = open_slice_data(h5f, islice)
        matrix_t[:, :, islice - 1, 0] = fact * tmp_data[0::2][:npt * npt].reshape(npt, npt) / (1000. * mesh_size)
        matrix_t[:, :, islice - 1, 1] = fact * tmp_data[1::2][:npt * npt].reshape(npt, npt) / (1000. * mesh_size)
    return matrix_t

def read_genesis_file(gen_fname):
    speed_of_light = const.value('speed of light in vacuum')
    h_eV_s = const.value('Planck constant in eV/Hz')
    lmb_und = 2.75e-2

    wf = Wavefront()

    ### Open the hdf5 output GENESIS file
    with h5py.File(gen_fname, 'r') as hf:
        grid_size = hf['gridsize'][()]
        slice_count = hf['slicecount'][()]
        print(slice_count)
        slice_spacing = hf['slicespacing'][()]
        wavelength = hf['wavelength'][()]
        data1 = open_slice_data(hf, 1)
        npoints = int(np.sqrt(len(data1) / 2))

    ### Definition of the Electric field field arrays where the electric field from the GENESIS output file
    ###  will be copied

    ### Fill in the fields of the wavefront object (mesh first: the field arrays must match it)
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
        wf.data.arrEver = np.zeros_like(wf.data.arrEhor)

    return wf
