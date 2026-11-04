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


# Known errors of SRW 4.2.1 (master 3c38a0d) on GPU, measured with a CUDA build on RTX 3090:
# without semi-analytical treatment the drift length is ignored for several slices,
# a single slice comes back with x and y transposed. Without CUDA SRW falls back to CPU and all cases pass.
GPU_CASES = [
    pytest.param(3, 1, id='3 slices, semi-analytical'),
    pytest.param(3, 0, id='3 slices, no semi-analytical',
                 marks=pytest.mark.xfail(reason='SRW GPU ignores drift length without semi-analytical treatment')),
    pytest.param(1, 1, id='1 slice',
                 marks=pytest.mark.xfail(reason='SRW GPU transposes x and y of a single-slice wavefront')),
]


@pytest.mark.parametrize('nz, semi', GPU_CASES)
def test_gpu_option(nz, semi):
    """gpu=1 gives the same intensity as CPU (with srwpy without CUDA SRW falls back to CPU)."""
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP
    from srwpy.srwlib import srwl

    res = []
    for gpu in (None, 1):
        wf = _gauss_wavefront(nx=200, nz=nz)
        b = Beamline()
        b.append(Drift(50), Use_PP(semi_analytical_treatment=semi))
        if nz > 1:
            srwl.SetRepresElecField(wf._srwl_wf, 'f')
        b.propagate(wf, gpu=gpu)
        if nz > 1:
            srwl.SetRepresElecField(wf._srwl_wf, 't')
        res.append(wf.get_intensity())
    np.testing.assert_allclose(res[1], res[0], rtol=0, atol=1e-4 * res[0].max())


def test_beamline_edit_by_name():
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP

    b = Beamline()
    names = ['d1', 'd2', 'd3']
    for name in names:
        d = Drift(1)
        d.name = name
        b.append(d, Use_PP())
    assert b.index == {'d1': 0, 'd2': 1, 'd3': 2}

    new_pp = Use_PP(zoom=2)
    b.edit_propagation_parameters('d2', new_pp)
    assert b.propagation_options[0]['propagation_parameters'][1] == new_pp.get_srw_pp()

    b.crop_beamline('unknown')
    assert len(b.propagation_options[0]['optical_elements']) == 3
    b.crop_beamline('d2')
    assert [el.name for el in b.propagation_options[0]['optical_elements']] == ['d1', 'd2']
    assert len(b.propagation_options[0]['propagation_parameters']) == 2


def test_glossary_fields_are_static():
    """Fields are declared in node classes: a typo raises instead of creating a new attribute."""
    import wpg
    from wpg.glossary import WavefrontMesh, get_wf_fields

    wf = wpg.Wavefront()
    assert isinstance(type(wf.params.Mesh).__dict__['nx'], type(WavefrontMesh.nx))
    with pytest.raises(AttributeError):
        wf.params.nx = 40  # nx belongs to params.Mesh
    with pytest.raises(AttributeError):
        wf.params.Mesh.nxx = 40
    paths = [f.path for f in get_wf_fields()]
    assert len(paths) == len(set(paths)) == 35
    assert paths[:2] == ['params/nval', 'params/wSpace'] and paths[-2:] == ['data/arrEhor', 'data/arrEver']


def test_field_array_sets_mesh(tmp_path):
    """Assigning an EM field array of another size updates the mesh (before: stray params.nx, mesh unchanged)."""
    wf = _gauss_wavefront(nx=50, nz=2)
    new = np.arange(30 * 40 * 2 * 2, dtype='float32').reshape(30, 40, 2, 2)
    with pytest.warns(UserWarning, match='Mesh'):
        wf.data.arrEhor = new  # also warns that arrEver is reset
    assert (wf.params.Mesh.ny, wf.params.Mesh.nx, wf.params.Mesh.nSlices) == (30, 40, 2)
    np.testing.assert_array_equal(wf.data.arrEhor, new)
    # the other polarization follows the new mesh (zeros), so the wavefront can be stored
    assert wf.data.arrEver.shape == new.shape and not wf.data.arrEver.any()
    wf.store_hdf5(str(tmp_path / 'resized.h5'))


def test_tilted_plane_roundtrip(tmp_path):
    """Unit vectors of the observation plane are floats (before: truncated to int on loading)."""
    import wpg

    wf = _gauss_wavefront(nx=20, nz=1)
    wf.params.Mesh.hvx, wf.params.Mesh.hvz = 0.8, 0.6
    wf.params.Mesh.nvx, wf.params.Mesh.nvz = -0.6, 0.8
    path = str(tmp_path / 'tilted.h5')
    wf.store_hdf5(path)
    wf2 = wpg.Wavefront()
    wf2.load_hdf5(path)
    m = wf2.params.Mesh
    assert (m.hvx, m.hvz, m.nvx, m.nvz) == pytest.approx((0.8, 0.6, -0.6, 0.8))


def test_wavelength():
    wf = _gauss_wavefront(nx=20, nz=1)
    assert wf.wavelength == pytest.approx(12.398419843320026e-7 / 5000.0, rel=1e-9)


def test_wavefront_shortcuts():
    """wf.nx, wf.photonEnergy, wf.Mesh are read-only shortcuts of wf.params(.Mesh) fields."""
    from wpg.glossary import WavefrontMesh, WavefrontParams

    wf = _gauss_wavefront(nx=20, nz=3)  # 3 slices: sliceMax of a single slice is NaN
    for f in WavefrontMesh.fields():
        assert getattr(wf, f.name) == getattr(wf.params.Mesh, f.name)
    for f in WavefrontParams.fields():
        assert getattr(wf, f.name) == getattr(wf.params, f.name)
    assert wf.Mesh is wf.params.Mesh
    with pytest.raises(AttributeError):
        wf.nx = 10
    with pytest.raises(AttributeError):
        wf.paramz  # noqa: B018  (no fallback for unknown names)


def test_field_array_on_fresh_wavefront():
    """Assigning field arrays to a fresh Wavefront (arEx, arEy are None in SRW) sets the mesh."""
    import wpg

    wf = wpg.Wavefront()
    ex = np.ones((8, 6, 2, 2), dtype='float32')
    wf.data.arrEhor = ex
    wf.data.arrEver = 2 * ex
    assert (wf.params.Mesh.ny, wf.params.Mesh.nx, wf.params.Mesh.nSlices) == (8, 6, 2)
    np.testing.assert_array_equal(wf.data.arrEver, 2 * ex)
