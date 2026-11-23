"""Regression tests for the fixes in the helper part of WPG
(wpg/useful_code, wpg/converters, s2e/prop), see WPG_notes/review/helpers.md.

Item numbers (pN) refer to the tables of helpers.md.
"""
import contextlib
import importlib
import io
import os
import pickle
import re
import sys
import types

import h5py
import numpy as np
import pytest

LCLS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'utest', 'lcls')
NCAR = 139                 # lcls.out: number of grid points
DX_GENESIS = 0.7481e-05    # lcls.out: meshsize [m]
LENG = (NCAR - 1) * DX_GENESIS
XLAMDS = 1.5e-10           # lcls.out: wavelength [m]
P_OUT = 4.6728e9           # lcls.out: last record of the power column [W]

EKEV = 5.0
LAM = 12.398419843320026e-10 / EKEV


def _sig_at(sig0, z):
    """RMS intensity size of a Gaussian beam with waist RMS size sig0 at distance z."""
    th = LAM / (4 * np.pi * sig0)
    return np.sqrt(sig0 ** 2 + (th * z) ** 2)


def _wf_power(wf):
    """sum of intensity * cell area in mm^2 (power for a field in sqrt(W/mm^2))."""
    m = wf.params.Mesh
    dx_mm = (m.xMax - m.xMin) / (m.nx - 1) * 1e3
    dy_mm = (m.yMax - m.yMin) / (m.ny - 1) * 1e3
    return wf.get_intensity(polarization='horizontal').sum() * dx_mm * dy_mm


def _rms_xy(wf):
    i = wf.get_intensity(polarization='horizontal').sum(axis=-1)  # (ny, nx)
    m = wf.params.Mesh
    x = np.linspace(m.xMin, m.xMax, m.nx)
    y = np.linspace(m.yMin, m.yMax, m.ny)
    px = i.sum(axis=0)
    py = i.sum(axis=1)
    return (np.sqrt((px * x ** 2).sum() / px.sum()), np.sqrt((py * y ** 2).sum() / py.sum()))


# ------------------------------------------------------------------ GENESIS 2 (ocelot dfl)

class _FakeDfl:
    """Minimal stand-in of ocelot RadiationField: fld[z, y, x], Lx() = leng as in ocelot."""

    def __init__(self, fld, leng):
        self.fld = fld
        self.dz = XLAMDS
        self.xlamds = XLAMDS
        self._leng = leng
        self.dx = leng / fld.shape[2]   # ocelot convention, not used by the converter

    def Nz(self):
        return self.fld.shape[0]

    def Nx(self):
        return self.fld.shape[2]

    def Lx(self):
        return self._leng


@pytest.fixture
def genesis_v2_stub(monkeypatch):
    """genesis_v2 imported against a stub ocelot reader; yields (module, state dict)."""
    state = {}
    ocelot = types.ModuleType('ocelot')
    adaptors = types.ModuleType('ocelot.adaptors')
    gen = types.ModuleType('ocelot.adaptors.genesis')
    gen.read_out_file = lambda path: path
    gen.read_dfl_file_out = lambda out, filePath, debug=0: state['dfl']
    monkeypatch.setitem(sys.modules, 'ocelot', ocelot)
    monkeypatch.setitem(sys.modules, 'ocelot.adaptors', adaptors)
    monkeypatch.setitem(sys.modules, 'ocelot.adaptors.genesis', gen)
    _forget_genesis_v2()
    mod = importlib.import_module('wpg.converters.genesis_v2')
    yield mod, state
    _forget_genesis_v2()


def _forget_genesis_v2():
    """Drop genesis_v2 from sys.modules and from the package namespace (it is bound to a reader)."""
    sys.modules.pop('wpg.converters.genesis_v2', None)
    import wpg.converters
    if hasattr(wpg.converters, 'genesis_v2'):
        delattr(wpg.converters, 'genesis_v2')


def _read_lcls_dfl():
    raw = np.fromfile(os.path.join(LCLS_DIR, 'lcls.out.dfl'), dtype=np.complex128)
    return raw.reshape(-1, NCAR, NCAR)  # (Nz, Ny, Nx), x fastest, as ocelot reads it


def test_genesis_v2_power_lcls(genesis_v2_stub):
    """p1: the dfl field is already in sqrt(W) per cell, the WPG power must equal lcls.out power.
    Before the fix: 1.03e-13 W instead of 4.67e9 W."""
    mod, state = genesis_v2_stub
    fld = _read_lcls_dfl()
    assert np.isclose((np.abs(fld) ** 2).sum(), P_OUT, rtol=1e-4)
    state['dfl'] = _FakeDfl(fld, LENG)
    with contextlib.redirect_stdout(io.StringIO()):
        wf = mod.read_genesis_file('lcls.out', 'lcls.out.dfl')
    m = wf.params.Mesh
    assert np.isclose((m.xMax - m.xMin) / (m.nx - 1), DX_GENESIS, rtol=1e-6)
    assert np.isclose(_wf_power(wf), P_OUT, rtol=1e-4)


def test_genesis_v2_axes_order(genesis_v2_stub):
    """p2: a beam wide in x must stay wide in x (WPG layout ny, nx, nSlices).
    Before the fix: rms_x = 37.4 um, rms_y = 143.2 um for sigma_x = 149.6 um, sigma_y = 37.4 um."""
    mod, state = genesis_v2_stub
    n = 101
    x = (np.arange(n) - n // 2) * DX_GENESIS
    X, Y = np.meshgrid(x, x)  # X varies along the last axis
    sx, sy = 20 * DX_GENESIS, 5 * DX_GENESIS
    f = np.exp(-X ** 2 / (4 * sx ** 2) - Y ** 2 / (4 * sy ** 2)).astype(np.complex128)
    state['dfl'] = _FakeDfl(f[None, :, :], (n - 1) * DX_GENESIS)
    with contextlib.redirect_stdout(io.StringIO()):
        wf = mod.read_genesis_file('x', 'y')
    rx, ry = _rms_xy(wf)
    assert np.isclose(ry, sy, rtol=0.01)
    assert 0.9 * sx < rx < sx   # rms_x is reduced by the grid window (+-2.5 sigma_x)


def test_genesis_v2_real_ocelot():
    """p1 with the real ocelot reader (skipped if ocelot-collab is not installed)."""
    pytest.importorskip('ocelot.adaptors.genesis')
    from ocelot.adaptors.genesis import read_out_file
    _forget_genesis_v2()
    genesis_v2 = importlib.import_module('wpg.converters.genesis_v2')
    out_path = os.path.join(LCLS_DIR, 'lcls.out')
    out = read_out_file(out_path)
    p_exit = float(np.ravel(out.p_int)[-1])  # out.power of ocelot is p_mid, not the power
    with contextlib.redirect_stdout(io.StringIO()):
        wf = genesis_v2.read_genesis_file(out_path, out_path + '.dfl')
    assert np.isclose(p_exit, P_OUT, rtol=1e-4)
    assert np.isclose(_wf_power(wf), p_exit, rtol=1e-4)


# ------------------------------------------------------------------ GENESIS HDF5 (genesis, genesis_v3)

def _make_genesis_h5(path, npt=8, nsl=4):
    with h5py.File(path, 'w') as f:
        f['gridsize'] = DX_GENESIS
        f['slicecount'] = nsl
        f['slicespacing'] = 32 * XLAMDS
        f['wavelength'] = XLAMDS
        for k in range(1, nsl + 1):
            v = np.zeros(2 * npt * npt)
            v[0::2] = k
            f['/slice%06d/field' % k] = v


@pytest.mark.parametrize('modname', ['genesis', 'genesis_v3'])
def test_genesis_h5_all_slices(tmp_path, modname):
    """p3: all GENESIS slices 1..nsl are read. Field k in slice k -> intensity k^2.
    Before the fix: 1 : 4 : 9 : 0."""
    mod = importlib.import_module('wpg.converters.' + modname)
    path = str(tmp_path / 'synthetic.h5')
    _make_genesis_h5(path)
    args = (path,) if modname == 'genesis' else (path, 0.03)
    with contextlib.redirect_stdout(io.StringIO()):
        wf = mod.read_genesis_file(*args)
    s = wf.get_intensity(polarization='horizontal').sum(axis=(0, 1))
    np.testing.assert_allclose(s / s[0], [1, 4, 9, 16], rtol=1e-5)


class _KeyRecorder(dict):
    def __getitem__(self, k):
        return {(): k}


@pytest.mark.parametrize('modname', ['genesis', 'genesis_v3'])
def test_open_slice_data_names(modname):
    """p14: six-digit slice group names; numbering from 1.
    Before the fix (genesis.py): '/slice0001000/field' for 1000, UnboundLocalError for 0."""
    mod = importlib.import_module('wpg.converters.' + modname)
    rec = _KeyRecorder()
    assert mod.open_slice_data(rec, 1) == '/slice000001/field'
    assert mod.open_slice_data(rec, 999) == '/slice000999/field'
    assert mod.open_slice_data(rec, 1000) == '/slice001000/field'
    with pytest.raises(ValueError):
        mod.open_slice_data(rec, 0)


def test_convert_from_phase(tmp_path):
    """p13: bytes/str mix. Before the fix: TypeError on the first dataset."""
    from wpg.converters.convert_from_phase import convert_from_phase
    src = str(tmp_path / 'in.h5')
    dst = str(tmp_path / 'out.h5')
    with h5py.File(src, 'w') as f:
        f['data/arrEhor'] = np.ones((2, 4, 3, 2), dtype='f4')
        f['data/arrEver'] = np.zeros((2, 4, 3, 2), dtype='f4')
        f['VERSION'] = 0.1
        f['nx'] = 3
    convert_from_phase(src, dst)
    with h5py.File(dst, 'r') as f:
        assert f['data/arrEx'].shape == (4, 6, 2)
        assert f['data/arrEy'].shape == (4, 6, 2)
        assert f['header/nx'][()] == 3
        assert np.isclose(f['version'][()], 0.1)


# ------------------------------------------------------------------ wfrutils, backpropagation

SIG0 = 10e-6     # waist RMS size along y [m], sigma_x0 = 2 SIG0
D = 10.0         # distance to the waist [m]
YOFF = 60e-6     # beam offset in y [m]
NXY = 201


@pytest.fixture(scope='module')
def asym_wf():
    """Single-slice Gaussian, offset in y, sigma_x = 2 sigma_y at the waist, y range 2x wider than x."""
    import wpg
    from wpg.generators import build_gauss_wavefront_xy
    sx0, sy0 = 2 * SIG0, SIG0
    sx = _sig_at(sx0, D)
    xmax, ymax = 5 * sx, 10 * sx
    srw = build_gauss_wavefront_xy(NXY, NXY, EKEV, -xmax, xmax, -ymax, ymax, sx0, sy0, D, yoff=YOFF)
    return wpg.Wavefront(srw)


def _asym_truth(wf):
    m = wf.params.Mesh
    dx = (m.xMax - m.xMin) / (m.nx - 1)
    dy = (m.yMax - m.yMin) / (m.ny - 1)
    sx, sy = _sig_at(2 * SIG0, D), _sig_at(SIG0, D)
    return dict(dx=dx, dy=dy, sx_px=sx / dx, sy_px=sy / dy, cx_px=(m.nx - 1) / 2,
                cy_px=(m.ny - 1) / 2 + YOFF / dy, sx=sx, sy=sy)


def test_calculate_peak_pos_and_fwhm_x(asym_wf):
    """p4: y0 from the y axis. Before the fix: y0 = 30.1 um instead of 60 um, calculate_fwhm_x = 0."""
    from wpg.useful_code import wfrutils
    t = _asym_truth(asym_wf)
    x0, y0 = wfrutils.calculate_peak_pos(asym_wf)
    assert abs(x0) <= t['dx']
    assert abs(y0 - YOFF) <= t['dy']
    fwhm_x = wfrutils.calculate_fwhm_x(asym_wf)
    assert abs(fwhm_x - 2.3548 * t['sx']) <= 2 * t['dx']


def test_fit_gaussian_axes(asym_wf):
    """p6, p18: x is the column (axis 1), y the row (axis 0).
    Before the fix: width_x = 9.92, width_y = 20.0 px; centre_x = 126.9, centre_y = 100.0 px;
    initial guess of moments 28.7 and 33.5 px."""
    from wpg.useful_code import backpropagation as bp
    t = _asym_truth(asym_wf)
    i2 = asym_wf.get_intensity(slice_number=0, polarization='horizontal')
    h, cx, cy, wx, wy = bp.fit_gaussian(i2)['params']
    assert np.isclose(abs(wx), t['sx_px'], rtol=0.01)
    assert np.isclose(abs(wy), t['sy_px'], rtol=0.01)
    assert abs(cx - t['cx_px']) < 0.1
    assert abs(cy - t['cy_px']) < 0.1
    # the fitted function evaluated on np.indices reproduces the data
    fit = bp.gaussian(h, cx, cy, wx, wy)(*np.indices(i2.shape))
    assert np.abs(fit - i2).max() < 1e-3 * i2.max()
    m = bp.moments(i2)
    assert np.isclose(m[3], t['sx_px'], rtol=0.05)
    assert np.isclose(m[4], t['sy_px'], rtol=0.05)


def test_fit_gaussian_rsquared_standard(asym_wf):
    """p17: R^2 = 1 - SS_res / SS_tot(data). Before the fix: 0.8286 instead of 0.8540 (5 % noise)."""
    from wpg.useful_code import backpropagation as bp
    i2 = asym_wf.get_intensity(slice_number=0, polarization='horizontal')
    noisy = i2 + np.random.default_rng(1).normal(0, 0.05 * i2.max(), i2.shape)
    r = bp.fit_gaussian(noisy)
    fit = bp.gaussian(*r['params'])(*np.indices(noisy.shape))
    r2_std = 1 - ((noisy - fit) ** 2).sum() / ((noisy - noisy.mean()) ** 2).sum()
    assert np.isclose(r['rsquared'], r2_std, rtol=1e-4)


def test_fit_gaussian_pulse_variance(asym_wf):
    """p6, p7: gaussian_parameters x/y on their places, pulse_variance_x/y are spatial variances in m^2.
    Before the fix: pulse_variance_x is an array of 201 numbers ~1e-42."""
    from wpg.useful_code import backpropagation as bp
    t = _asym_truth(asym_wf)
    bp.fit_gaussian_pulse(asym_wf)
    gp = asym_wf.custom_fields['misc']['gaussian_parameters']
    assert np.ndim(gp['pulse_variance_x']) == 0
    assert np.isclose(gp['pulse_variance_x'], t['sx'] ** 2, rtol=0.02)
    assert np.isclose(gp['pulse_variance_y'], t['sy'] ** 2, rtol=0.02)
    assert np.isclose(gp['pulse_variance'], t['sx'] ** 2 + t['sy'] ** 2, rtol=0.02)
    assert np.isclose(abs(gp['width_x'][0]), t['sx_px'], rtol=0.01)
    assert abs(gp['pulse_center_y'] - t['cy_px']) < 0.1


@pytest.fixture
def gauss3_wf():
    import wpg
    from wpg.generators import build_gauss_wavefront
    sx0 = 2 * SIG0
    xmax = 5 * _sig_at(sx0, D)
    return wpg.Wavefront(build_gauss_wavefront(64, 64, 3, EKEV, -xmax, xmax, -xmax, xmax, 0.12e-15,
                                               sx0, sx0, D))


def test_show_slices_hsv_runs(gauss3_wf):
    """p10: mesh steps without floor division. Before the fix: dt = 0, ValueError (NaN to integer)."""
    import matplotlib.pyplot as plt
    from wpg.useful_code import wfrutils
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            wfrutils.show_slices_hsv(gauss3_wf, slice_numbers=1)
    finally:
        plt.close('all')
    dt = float(re.search(r'^dt (\S+)', buf.getvalue(), re.M).group(1))
    m = gauss3_wf.params.Mesh
    assert np.isclose(dt, (m.sliceMax - m.sliceMin) / (m.nSlices - 1))


def _printed_total_power(wf):
    import matplotlib.pyplot as plt
    from wpg.useful_code import wfrutils
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            wfrutils.plot_wfront(wf, 'at 10 m', False, False, 1e-3, 1e-3, 'x', True)
    finally:
        plt.close('all')
    return float(re.search(r'Total power.*?:\s*(\S+) \[GW\]', buf.getvalue()).group(1))


def test_plot_wfront_power_units(gauss3_wf):
    """p5: the photon energy factor only for sqrt(Phot/s/0.1%bw/mm^2).
    Before the fix: 8.93e-16 GW printed for a field in sqrt(W/mm^2) with direct sum 1.115 GW."""
    ii = gauss3_wf.get_intensity(slice_number=0, polarization='horizontal')
    m = gauss3_wf.params.Mesh
    direct_gw = ii.sum() * (m.xMax - m.xMin) / (m.nx - 1) * (m.yMax - m.yMin) / (m.ny - 1) * 1e6 * 1e-9
    gauss3_wf.params.wEFieldUnit = 'sqrt(W/mm^2)'
    assert np.isclose(_printed_total_power(gauss3_wf), direct_gw, rtol=1e-4)
    gauss3_wf.params.wEFieldUnit = 'sqrt(Phot/s/0.1%bw/mm^2)'
    photon_j = gauss3_wf.params.photonEnergy / 6.24150934e18
    assert np.isclose(_printed_total_power(gauss3_wf), direct_gw * photon_j, rtol=1e-4)


# ------------------------------------------------------------------ s2e/prop

def test_add_history_external_link_relative(tmp_path, monkeypatch):
    """p8: the link to the parent /data is relative to the folder of the output file.
    Before the fix: link 'in/src.h5' (relative to the run folder), KeyError when read from out/."""
    from s2e.prop import propagate_s2e as ps
    monkeypatch.chdir(tmp_path)
    os.makedirs('in')
    os.makedirs('out')
    with h5py.File('in/src.h5', 'w') as f:
        f['data/arrEhor'] = np.arange(6.)
        f['misc/temporal_struct'] = np.ones((3, 2))
    with h5py.File('out/prop.h5', 'w') as f:
        f['data/arrEhor'] = np.zeros(2)
    ps.add_history('out/prop.h5', 'in/src.h5')
    with h5py.File('out/prop.h5', 'r') as f:
        assert f['history/parent/detail'].get('data', getlink=True).filename == '../in/src.h5'
    for cwd in (tmp_path, tmp_path / 'out', tmp_path / 'in'):
        monkeypatch.chdir(cwd)
        with h5py.File(str(tmp_path / 'out' / 'prop.h5'), 'r') as f:
            np.testing.assert_array_equal(f['history/parent/detail/data/arrEhor'][()], np.arange(6.))
            assert f['history/parent/detail/misc/temporal_struct'].shape == (3, 2)


BEAMLINE_SRC = '''
from wpg import Beamline
from wpg.optical_elements import Drift, Use_PP
def get_beamline():
    bl = Beamline()
    bl.append(Drift(1.0), Use_PP(semi_analytical_treatment=1))
    return bl
'''


@pytest.fixture(scope='module')
def s2e_run(tmp_path_factory):
    """One Gaussian pulse propagated through a 1 m drift by directory_process (worker pool)."""
    import wpg
    from wpg.generators import build_gauss_wavefront
    from s2e.prop import propagate_s2e as ps
    root = tmp_path_factory.mktemp('s2e')
    (root / 'in').mkdir()
    sx0 = 2 * SIG0
    xmax = 5 * _sig_at(sx0, D)
    wf0 = wpg.Wavefront(build_gauss_wavefront(64, 64, 3, EKEV, -xmax, xmax, -xmax, xmax, 0.12e-15,
                                              sx0, sx0, D))
    t = np.linspace(-1, 1, 5)
    wf0.custom_fields['misc'] = {'temporal_struct': np.stack([t, np.exp(-t ** 2)], axis=1)}
    wf0.store_hdf5(str(root / 'in' / 'FELsource_out_0000001.h5'))
    bl_file = root / 'beamline.py'
    bl_file.write_text(BEAMLINE_SRC)
    with contextlib.redirect_stdout(io.StringIO()):
        ps.directory_process(str(root / 'in'), str(root / 'out'), str(bl_file), 1)
    return root


def test_beamline_function_pickles(s2e_run):
    """p9: get_beamline loaded from a file can be pickled (registered in sys.modules).
    Before the fix: PicklingError "import of module 'custom_beamline' failed"."""
    from s2e.prop import propagate_s2e as ps
    get_beamline = ps.load_beamline_function(str(s2e_run / 'beamline.py'))
    assert pickle.loads(pickle.dumps(get_beamline)) is get_beamline


def test_directory_process_loaded_function(s2e_run, tmp_path):
    """The batch mode accepts get_beamline loaded by load_beamline_function: the 'spawn' workers
    load the same file. The pool used fork before, and with the OpenMP build of SRW the workers hung."""
    from s2e.prop import propagate_s2e as ps
    get_beamline = ps.load_beamline_function(str(s2e_run / 'beamline.py'))
    with contextlib.redirect_stdout(io.StringIO()):
        ps.directory_process(str(s2e_run / 'in'), str(tmp_path), get_beamline, 1)
    assert (tmp_path / 'prop_out_0000001.h5').exists()


def test_directory_process_rejects_notebook_function(monkeypatch):
    """A function of a notebook cell cannot be imported by a 'spawn' worker: clear error, not a hang."""
    from s2e.prop import propagate_s2e as ps

    def get_beamline():
        pass
    get_beamline.__module__ = '__main__'
    monkeypatch.delattr(sys.modules['__main__'], '__file__', raising=False)
    with pytest.raises(ValueError, match='notebook'):
        ps.directory_process('in', 'out', get_beamline, 1)


def test_directory_process_output(s2e_run):
    """p9: batch mode with a pool of workers produces the output with a readable parent link."""
    out = s2e_run / 'out' / 'prop_out_0000001.h5'
    assert out.exists()
    with h5py.File(str(out), 'r') as f:
        assert f['history/parent/detail/data/arrEhor'].size > 0


def test_calculate_fwhm_step(s2e_run):
    """p15: mesh step (xMax - xMin) / (nx - 1). Before the fix: step range / nx (63/64 of the true one)."""
    import wpg
    from s2e.prop import propagate_s2e as ps
    wf = wpg.Wavefront()
    wf.load_hdf5(str(s2e_run / 'out' / 'prop_out_0000001.h5'))
    m = wf.params.Mesh
    i_tot = wf.get_intensity(polarization='total').sum(axis=-1)
    iy, ix = np.unravel_index(np.argmax(i_tot), i_tot.shape)
    row = i_tot[iy]
    n_half = int((row > row.max() / 2).sum())
    fw = ps.calculate_fwhm(wf)
    assert np.isclose(fw['fwhm_x'], n_half * (m.xMax - m.xMin) / (m.nx - 1), rtol=1e-9)
    assert np.isclose(wf.custom_fields['misc']['xFWHM'], fw['fwhm_x'], rtol=1e-9)


def test_show_diagnostics(s2e_run, monkeypatch):
    """p11, p12: file name and history path. Before the fix: 'out/prop_out_0000001.h5.h5 not found',
    KeyError 'parent', then AttributeError plt.hold."""
    import matplotlib.pyplot as plt
    from s2e.prop import diagnostics
    assert diagnostics._prop_out_file_name('out/prop_out_0000001.h5') == 'out/prop_out_0000001.h5'
    assert diagnostics._prop_out_file_name('1') == 'prop_out_0000001.h5'
    assert diagnostics._prop_out_file_name('out/prop_out_0000001') == 'out/prop_out_0000001.h5'
    monkeypatch.chdir(s2e_run)
    plt.close('all')
    buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf):
            diagnostics.show_diagnostics('out/prop_out_0000001.h5')
        assert 'not found' not in buf.getvalue()
        assert 'No history' not in buf.getvalue()
        assert len(plt.get_fignums()) == 3
        ax = plt.figure(2).axes[0]
        assert [ln.get_label() for ln in ax.get_lines()] == ['original', 'propag']
    finally:
        plt.close('all')
