import scipy.constants as const
import numpy as np
import numpy
import h5py
from wpg import Wavefront
from ocelot.adaptors.genesis import read_out_file, read_dfl_file_out

def vector_grid_conversion(data, _npoints, _nslices, _grid_size, _wv, _lambda_un):
    '''
    Inputs:
          data ( radiation field data [nz, ny, nx])
          _npoints (Number of points in the grid, equivalent to ncar in GENESIS)
          _nslices (Number of slices)
          _grid_size (grid size extracted from h5)
          _wv radiation wavelength
          _lambda_un Undulator period

    Output:
          matrix_t   Matrix that contains the Electric field data along the grid
          (Real and imaginary parts)shape = (npoints,npoints, slice_count,2)))

    Based on WPG/wpg/converters/genesis.py
    '''
    # The dfl field of GENESIS 2 is already in sqrt(W) per grid cell: sum |E|^2 over the grid
    # is the radiation power (checked with ocelot on tests/utest/lcls: 4.6728e9 W in both).
    # Hence no conversion factor from the dimensionless GENESIS field is applied here
    # (_wv and _lambda_un are kept in the signature for compatibility, they are not used);
    # the field is only divided by the mesh step in mm to get sqrt(W/mm^2).
    npt = _npoints
    nsl = _nslices
    mesh_size = _grid_size / (npt - 1)
    matrix_t = numpy.zeros(shape=(npt, npt, nsl, 2))

    # Slice out real and imaginary parts and store in wavefront matrix.
    # radiation field axes (ocelot) are z, y, x; wavefront axes are y, x, z, (re/im)
    matrix_t[:, :, :, 0] = data.real.transpose(1, 2, 0)
    matrix_t[:, :, :, 1] = data.imag.transpose(1, 2, 0)

    matrix_t = matrix_t / (1000. * mesh_size)  # 1e3 converts from m to mm.

    ### ORIGINAL CODE
    #for z_slice in range(_nslices):
        ##print('slice No ' + str(islice))
        #tmp_data = data[z_slice]
        #for jc in range(npt):
            #for lc in range(npt):
                #ind = 0
                #while ind < 2:
                    #if ind == 0:
                        #vect = tmp_data[0::2]
                    #elif ind == 1:
                        #vect = tmp_data[1::2]
                    #matrix_t[jc, lc, islice - 1, ind] = fact * \
                        #vect[(jc * npt) + lc] / (1000. * mesh_size)
                    #ind = ind + 1
    ### END ORIGINAL CODE

    return matrix_t

def read_genesis_file(genesis_out, genesis_dfl):
    '''
    Based on WPG/wpg/converters/genesis.py
    '''

    # Needed constants.
    speed_of_light = const.codata.value('speed of light in vacuum')
    h_eV_s = const.codata.value('Planck constant in eV s')

    # Undulator period: not used any more (the dfl field is already in sqrt(W) per cell).
    lmb_und = 2.75e-2

    # Initialize empty wavefront.
    wf = Wavefront()

    # Get the radiation field from genesis output.
    genesis_out = read_out_file(genesis_out)

    genesis_radiation_field = read_dfl_file_out(out=genesis_out, filePath=genesis_dfl, debug=0)

    # Extract geometry.
    slice_count = genesis_radiation_field.Nz()
    slice_spacing = genesis_radiation_field.dz
    wavelength = genesis_radiation_field.xlamds
    npoints = genesis_radiation_field.Nx()

    # Fill in the fields of the wavefront object (mesh first: the field arrays must match it)
    wf.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wf.params.photonEnergy = h_eV_s * speed_of_light / wavelength
    wf.params.wDomain = 'time'
    wf.params.Mesh.nSlices = slice_count
    wf.params.Mesh.nx = npoints
    wf.params.Mesh.ny = npoints
    pulse_length = (slice_count - 1) * slice_spacing / (speed_of_light)
    wf.params.Mesh.sliceMin = -pulse_length / 2.
    wf.params.Mesh.sliceMax = pulse_length / 2.
    # ocelot Lx() is the full GENESIS grid width leng = 2*dgrid = (ncar - 1) * meshsize,
    # so the WPG step range_xy / (npoints - 1) equals the GENESIS meshsize
    # (ocelot's own dfl.dx = Lx()/Nx is smaller by 1/ncar and is not used).
    range_xy = genesis_radiation_field.Lx()
    wf.params.Mesh.xMin = -range_xy / 2.
    wf.params.Mesh.xMax = range_xy / 2.
    wf.params.Mesh.yMin = -range_xy / 2.
    wf.params.Mesh.yMax = range_xy / 2.


    print ('Photon energy: ', wf.params.photonEnergy, 'eV')
    print ('Pulse length: %4.3e' % ( pulse_length ))

    # Extract the field data from the h5 file and fill in the data of the Electric field, by calling the
    # vector_grid_conversion function
    wf.data.arrEhor = vector_grid_conversion(
        genesis_radiation_field.fld, npoints, slice_count, range_xy, wavelength, lmb_und)
    wf.data.arrEver = numpy.zeros_like(wf.data.arrEhor)

    return wf


