def test_import_wpg():
    import sys
    sys.path.insert(0, '..')
    import wpg


def test_import_wpg_members():
    import sys
    sys.path.insert(0, '..')
    from wpg import Wavefront, Beamline
    import wpg.beamline
    import wpg.glossary
    import wpg.generators
    import wpg.optical_elements
    import wpg.srwlib
    import wpg.srwlpy
    import wpg.utils
    import wpg.wpg_uti_exfl
    import wpg.wpg_uti_oe
    import wpg.wpg_uti_wf


def test_simple_gauusina_propagation():
    # TODO: fix propagation for coerrect results
    import sys
    sys.path.insert(0, '..')
    import os

    import wpg
    from wpg.generators import build_gauss_wavefront
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP
    from wpg.srwlib import srwl

    import numpy as np

    d2waist = 270.
    # beam parameters:
    qnC = 0.1  # [nC] e-bunch charge
    thetaOM = 3.6e-3
    ekev = 5.0

    # calculate angular divergence:
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC)) * 1e-6 / ekev ** 0.85
    theta_rms = theta_fwhm / 2.35
    sigX = 12.4e-10 / (ekev * 4 * np.pi * theta_rms)

    # define limits
    xmax = theta_rms * d2waist * 3.5
    xmin = - xmax
    ymin = xmin
    ymax = xmax
    nx = 300
    ny = nx
    nz = 3
    tau = 0.12e-15

    srw_wf = build_gauss_wavefront(
        nx, ny, nz, ekev, xmin, xmax, ymin, ymax, tau, sigX, sigX, d2waist)
    wf = wpg.Wavefront(srw_wf)
    b = Beamline()
    b.append(Drift(5), Use_PP())
    srwl.SetRepresElecField(wf._srwl_wf, 'f')
    b.propagate(wf)
    srwl.SetRepresElecField(wf._srwl_wf, 'c')
    srwl.ResizeElecField(srw_wf, 'c', [0, 0.25, 1, 0.25, 1])

    ti = wf.get_intensity()
    
    out_folder = os.path.join(os.path.dirname(__file__), 'tests_data')
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    wf_hdf5_out_file_path = os.path.join(out_folder, 'my_gauss.h5')
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf_out = wpg.Wavefront()
    wf_out.load_hdf5(wf_hdf5_out_file_path)
    return wf


def test_hisotry():
    import sys
    sys.path.insert(0, '..')
    import os

    import wpg
    from wpg.generators import build_gauss_wavefront
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP
    # from wpg.srwlib import srwl

    import numpy as np

    d2waist = 270.
    # beam parameters:
    qnC = 0.1  # [nC] e-bunch charge
    thetaOM = 3.6e-3
    ekev = 5.0

    # calculate angular divergence:
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC)) * 1e-6 / ekev ** 0.85
    theta_rms = theta_fwhm / 2.35
    sigX = 12.4e-10 / (ekev * 4 * np.pi * theta_rms)

    # define limits
    xmax = theta_rms * d2waist * 3.5
    xmin = - xmax
    ymin = xmin
    ymax = xmax
    nx = 300
    ny = nx
    nz = 3
    tau = 0.12e-15

    srw_wf = build_gauss_wavefront(
        nx, ny, nz, ekev, xmin, xmax, ymin, ymax, tau, sigX, sigX, d2waist)
    wf = wpg.Wavefront(srw_wf)
    b = Beamline()
    # b.append(Drift(5), Use_PP())
    # b.propagate(wf)
    # srwl.ResizeElecField(srw_wf, 'c', [0, 0.25, 1, 0.25, 1])

    out_folder = os.path.join(os.path.dirname(__file__), 'tests_data')
    if not os.path.exists(out_folder):
        os.mkdir(out_folder)

    wf_hdf5_out_file_path = os.path.join(out_folder, 'my_gauss_history.h5')
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf.custom_fields['/history/params/beamline/printout'] = str(b)
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf_out = wpg.Wavefront()
    wf_out.load_hdf5(wf_hdf5_out_file_path)
    wf_out.custom_fields['/history/params/beamline/printout'] = str(b)
    wf_out.store_hdf5(wf_hdf5_out_file_path)
    return wf
