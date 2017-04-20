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
    # Definition of constants
    vac_imp = const.codata.value('characteristic impedance of vacuum')
    eev = 1e6 * const.codata.value('electron mass energy equivalent in MeV')

    # Definition of internal variables
    npt = _npoints
    nsl = _nslices
    mesh_size = _grid_size / (npt - 1)
    lmb = _wv
    lmb_u = _lambda_un
    matrix_t = numpy.zeros(shape=(npt, npt, nsl, 2))


    # Definition of parameters needed for the scaling factor to \sqrt{W} units
    xkw0 = 2. * numpy.pi / lmb_u
    xks = 2. * numpy.pi / lmb
    dxy = xkw0 * mesh_size

    # Scaling factor in order to get the field in units of \sqrt{W}
    # (consistent with GENESIS v2)
    fact = dxy * eev * xkw0 / xks / numpy.sqrt(vac_imp)
    fact = fact / (xkw0 * xkw0)
    print fact

    # Cycle over all slices, converting the 1D vector into a 3D Matrix and
    # dividing by the mesh size in order to get the field in units of \sqrt{W/mm^2}

    # Slice out real and imaginary parts and store in wavefront matrix.
    # radiation field axes are z,y,x. Wavefront axes are x,y,z, (re/im)
    matrix_t[:,:,:,0] = data.real.swapaxes(0,2)
    matrix_t[:,:,:,1] = data.imag.swapaxes(0,2)

    matrix_t  = matrix_t * fact / ( 1000. * mesh_size) # 1e3 converts from m to mm.

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

    # Hard-code undulator period, should be read from genesis out object.
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

    # Definition of the Electric field arrays where the electric field from the GENESIS output file
    #  will be copied

    # Setup E-field.
    wf.data.arrEhor = numpy.zeros(shape=(npoints, npoints, slice_count, 2))
    wf.data.arrEver = numpy.zeros(shape=(npoints, npoints, slice_count, 2))

    # Fill in the fields of the wavefront object
    wf.params.wEFieldUnit = 'sqrt(W/mm^2)'
    wf.params.photonEnergy = h_eV_s * speed_of_light / wavelength
    wf.params.wDomain = 'time'
    wf.params.Mesh.nSlices = slice_count
    wf.params.Mesh.nx = npoints
    wf.params.Mesh.ny = npoints
    pulse_length = (slice_count - 1) * slice_spacing / (speed_of_light)
    wf.params.Mesh.sliceMin = -pulse_length / 2.
    wf.params.Mesh.sliceMax = pulse_length / 2.
    range_xy = genesis_radiation_field.Lx()
    wf.params.Mesh.xMin = -range_xy / 2.
    wf.params.Mesh.xMax = range_xy / 2.
    wf.params.Mesh.yMin = -range_xy / 2.
    wf.params.Mesh.yMax = range_xy / 2.


    print 'Photon energy: ', wf.params.photonEnergy, 'eV'
    print 'Pulse length: %4.3e' % ( pulse_length )

    # Extract the field data from the h5 file and fill in the data of the Electric field, by calling the
    # vector_grid_conversion function
    wf.data.arrEhor = vector_grid_conversion(
        genesis_radiation_field.fld, npoints, slice_count, range_xy, wavelength, lmb_und)

    return wf


