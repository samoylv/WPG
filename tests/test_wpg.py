import importlib
import sys

import numpy as np
import pytest


def test_import_wpg():
    import wpg  # noqa: F401


def test_import_wpg_members():
    from wpg import Wavefront, Beamline  # noqa: F401
    import wpg.beamline  # noqa: F401
    import wpg.glossary  # noqa: F401
    import wpg.generators  # noqa: F401
    import wpg.optical_elements  # noqa: F401
    import wpg.utils  # noqa: F401
    import wpg.wpg_uti_exfl  # noqa: F401
    import wpg.wpg_uti_oe  # noqa: F401
    import wpg.wpg_uti_wf  # noqa: F401
    import srwpy.srwlib  # noqa: F401


def test_srwlib_compat_module():
    """wpg.srwlib is kept as a deprecated alias of srwpy.srwlib."""
    import srwpy.srwlib
    sys.modules.pop('wpg.srwlib', None)
    with pytest.warns(DeprecationWarning, match='srwpy.srwlib'):
        compat = importlib.import_module('wpg.srwlib')
    assert compat.srwl is srwpy.srwlib.srwl
    assert compat.SRWLOptD is srwpy.srwlib.SRWLOptD


def _gauss_wavefront(nx=300, nz=3):
    import wpg
    from wpg.generators import build_gauss_wavefront

    d2waist = 270.
    # beam parameters:
    qnC = 0.1  # [nC] e-bunch charge
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
    ny = nx
    tau = 0.12e-15

    srw_wf = build_gauss_wavefront(
        nx, ny, nz, ekev, xmin, xmax, ymin, ymax, tau, sigX, sigX, d2waist)
    return wpg.Wavefront(srw_wf)


def test_simple_gaussian_propagation(tmp_path):
    import wpg
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP
    from srwpy.srwlib import srwl

    wf = _gauss_wavefront()
    b = Beamline()
    b.append(Drift(5), Use_PP())
    srwl.SetRepresElecField(wf._srwl_wf, 'f')
    b.propagate(wf)
    srwl.SetRepresElecField(wf._srwl_wf, 'c')
    srwl.ResizeElecField(wf._srwl_wf, 'c', [0, 0.25, 1, 0.25, 1])

    ti = wf.get_intensity()
    assert ti.shape == (wf.params.Mesh.ny, wf.params.Mesh.nx, wf.params.Mesh.nSlices)
    assert np.all(np.isfinite(ti))
    assert ti.max() > 0

    wf_hdf5_out_file_path = tmp_path / 'my_gauss.h5'
    wf.store_hdf5(str(wf_hdf5_out_file_path))

    wf_out = wpg.Wavefront()
    wf_out.load_hdf5(str(wf_hdf5_out_file_path))
    np.testing.assert_allclose(wf_out.get_intensity(), ti, rtol=1e-6)


def test_history(tmp_path):
    import wpg
    from wpg.beamline import Beamline

    wf = _gauss_wavefront()
    b = Beamline()

    wf_hdf5_out_file_path = str(tmp_path / 'my_gauss_history.h5')
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf.custom_fields['/history/params/beamline/printout'] = str(b)
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf_out = wpg.Wavefront()
    wf_out.load_hdf5(wf_hdf5_out_file_path)
    wf_out.custom_fields['/history/params/beamline/printout'] = str(b)
    wf_out.store_hdf5(wf_hdf5_out_file_path)


def test_hdf5_tree(tmp_path, capsys):
    from wpg.utils import load_hdf5_tree, print_hdf5

    wf = _gauss_wavefront(nx=50, nz=1)
    path = str(tmp_path / 'tree.h5')
    wf.store_hdf5(path)

    tree = load_hdf5_tree(path)
    assert len(tree) > 0
    print_hdf5(path)
    assert capsys.readouterr().out
