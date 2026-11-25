"""Regression tests for the fixes of review/optics.md (optical elements, sources, utilities).

Each test reproduces one finding of the review; the number of the finding is given in the docstring.
"""
import copy
import os
import warnings

import numpy as np
import pytest
from scipy.special import erf

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

EKEV = 8.0
S0 = 10e-6


def _gauss_t(nx=40, ny=40, nz=9, R=5 * S0, Ry=None, tau=5e-15, pulseEn=1e-3):
    """Time-domain Gaussian pulse (build_gauss_wavefront) as wpg.Wavefront."""
    from wpg.wavefront import Wavefront
    from wpg.generators import build_gauss_wavefront
    Ry = R if Ry is None else Ry
    return Wavefront(build_gauss_wavefront(nx, ny, nz, EKEV, -R, R, -Ry, Ry, tau, S0, S0, 0., pulseEn=pulseEn))


def _gauss_xy(nx=64, R=5 * S0, d2waist=0., **kw):
    from wpg.wavefront import Wavefront
    from wpg.generators import build_gauss_wavefront_xy
    return Wavefront(build_gauss_wavefront_xy(nx, nx, EKEV, -R, R, -R, R, S0, S0, d2waist, **kw))


# ---------------------------------------------------------------- 1. Mirror_plane_2d layout

@pytest.mark.parametrize('orient, nx, ny', [('x', 40, 40), ('x', 50, 10), ('y', 40, 40), ('y', 10, 50)])
def test_mirror_plane_2d_layout(tmp_path, orient, nx, ny):
    """1: OPD of a profile h = a*x_long must vary only along the reflection axis of the element."""
    from wpg.optical_elements import Mirror_plane_2d
    from wpg.wpg_uti_oe import get_opd
    a, length, width, theta = 1e-6, 0.4, 0.02, 2e-3
    xl = np.linspace(-0.25, 0.25, 101)
    yt = np.linspace(-0.015, 0.015, 31)
    X, Y = np.meshgrid(xl, yt, indexing='ij')
    fn = tmp_path / 'h2d.dat'
    np.savetxt(fn, np.column_stack([X.ravel(), Y.ravel(), a * X.ravel()]))

    op = Mirror_plane_2d(orient, theta, length, width, str(fn), nx=nx, ny=ny)
    assert (op.mesh.nx, op.mesh.ny) == (nx, ny)
    opd = get_opd(op)  # (ny, nx), x fastest -- SRW convention
    pv_along_x = np.mean(opd.max(axis=1) - opd.min(axis=1))
    pv_along_y = np.mean(opd.max(axis=0) - opd.min(axis=0))
    expected = 2 * np.sin(theta) * a * length
    along, across = (pv_along_x, pv_along_y) if orient == 'x' else (pv_along_y, pv_along_x)
    assert along == pytest.approx(expected, rel=1e-3), f'PV along mirror {along:.3e}, expected {expected:.3e}'
    assert across < 1e-6 * expected, f'PV across mirror {across:.3e} m, expected 0'


# ---------------------------------------------------------------- 2-4, 15. source.sase_pulse

def _sase(seed=0, **kw):
    pytest.importorskip('phenom')
    from wpg import source
    args = dict(x=np.linspace(-100e-6, 100e-6, 64), y=np.linspace(-300e-6, 300e-6, 48),
                t=np.linspace(-40e-15, 40e-15, 32), photon_energy=9e3, pulse_energy=1e-3,
                pulse_duration=10e-15, bandwidth=1e-3, sigma=20e-6, div=1e-6, domain='time')
    args.update(kw)
    np.random.seed(seed)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        return source.sase_pulse(**args), args


def test_sase_pulse_mesh_and_axes():
    """3: phenom returns (ny, nx, nt); mesh nx, ny and the beam size along x and y must be right."""
    wf, a = _sase()
    m = wf.params.Mesh
    assert (m.nx, m.ny, m.nSlices) == (64, 48, 32)
    assert (m.xMin, m.xMax, m.yMin, m.yMax) == pytest.approx((-100e-6, 100e-6, -300e-6, 300e-6))
    intens = wf.get_intensity().sum(axis=-1)
    assert intens.shape == (48, 64)
    x = np.linspace(m.xMin, m.xMax, m.nx)
    y = np.linspace(m.yMin, m.yMax, m.ny)
    ix, iy = intens.sum(axis=0), intens.sum(axis=1)
    rms_x = np.sqrt((ix * x ** 2).sum() / ix.sum())
    rms_y = np.sqrt((iy * y ** 2).sum() / iy.sum())
    # phenom field ~ exp(-r^2/(2 sigma^2)), intensity rms = sigma/sqrt(2)
    assert rms_x == pytest.approx(a['sigma'] / np.sqrt(2), rel=0.02), f'rms_x = {rms_x*1e6:.1f} um'
    assert rms_y == pytest.approx(a['sigma'] / np.sqrt(2), rel=0.02), f'rms_y = {rms_y*1e6:.1f} um'


def test_sase_pulse_polarization():
    """2: 'vertical' keeps the field in Ey, 'circular' is (Ex, iEx)/sqrt(2) with the same total energy."""
    wf_h, _ = _sase(polarization='horizontal')
    wf_v, _ = _sase(polarization='vertical')
    wf_c, _ = _sase(polarization='circular')
    ih = wf_h.get_intensity(polarization='horizontal').sum()
    assert wf_h.get_intensity(polarization='vertical').sum() == 0
    assert wf_v.get_intensity(polarization='horizontal').sum() == 0
    assert wf_v.get_intensity(polarization='vertical').sum() == pytest.approx(ih, rel=1e-5)
    ich = wf_c.get_intensity(polarization='horizontal').sum()
    icv = wf_c.get_intensity(polarization='vertical').sum()
    assert ich + icv == pytest.approx(ih, rel=1e-5), f'circular: {ich + icv:.4e}, horizontal: {ih:.4e}'
    assert ich == pytest.approx(icv, rel=1e-5)
    eh = wf_c.data.arrEhor[..., 0] + 1j * wf_c.data.arrEhor[..., 1]
    ev = wf_c.data.arrEver[..., 0] + 1j * wf_c.data.arrEver[..., 1]
    np.testing.assert_allclose(ev, 1j * eh, atol=1e-6 * np.abs(eh).max())


def test_sase_pulse_unknown_polarization():
    with pytest.raises(ValueError):
        _sase(polarization='elliptic')


@pytest.mark.parametrize('domain', ['time', 'freq'])
def test_sase_pulse_domain(domain):
    """4: wDomain and the slice mesh follow the domain of the phenom field."""
    wf, a = _sase(domain=domain)
    m = wf.params.Mesh
    if domain == 'time':
        assert wf.params.wDomain == 'time'
        assert (m.sliceMin, m.sliceMax) == pytest.approx((a['t'].min(), a['t'].max()))
    else:
        assert wf.params.wDomain == 'frequency'
        bw = a['photon_energy'] * a['bandwidth']
        assert (m.sliceMin, m.sliceMax) == pytest.approx((a['photon_energy'] - 3 * bw, a['photon_energy'] + 3 * bw))
    assert wf.params.photonEnergy == a['photon_energy']
    assert np.abs(wf.data.arrEhor).max() > 0


def test_sase_pulse_requires_grids():
    """15: x, y, t have no usable defaults -- clear error instead of TypeError deep inside."""
    pytest.importorskip('phenom')
    from wpg import source
    with pytest.raises(ValueError, match='x, y and t'):
        source.sase_pulse()


# ---------------------------------------------------------------- 5. propagate_wavefront

def test_propagate_wavefront_keeps_input():
    """5: propagate_wavefront must not change the input wavefront."""
    from wpg.beamline import Beamline
    from wpg.optical_elements import Drift, Use_PP
    from wpg import wpg_uti_wf
    wf_in = _gauss_xy(nx=64)
    z0 = wf_in.params.Mesh.zCoord
    imax0 = wf_in.get_intensity().max()
    bl = Beamline()
    bl.append(Drift(10.), Use_PP(semi_analytical_treatment=1))
    wf_out = wpg_uti_wf.propagate_wavefront(wf_in, bl)
    assert wf_out._srwl_wf is not wf_in._srwl_wf
    assert wf_in.params.Mesh.zCoord == z0, f'input z changed to {wf_in.params.Mesh.zCoord} m'
    assert wf_in.get_intensity().max() == pytest.approx(imax0)
    assert wf_out.params.Mesh.zCoord == pytest.approx(z0 + 10.)


# ---------------------------------------------------------------- 6. Mirror_plane x0, delim

def _bump_profile(path, delimiter=' '):
    xl = np.linspace(-0.25, 0.25, 501)
    h = 1e-9 * np.exp(-xl ** 2 / (2 * 0.02 ** 2))
    np.savetxt(path, np.column_stack([xl, h]), delimiter=delimiter)
    return str(path)


def test_mirror_plane_x0_shift(tmp_path):
    """6: x0 shifts the profile along the mirror: projected shift x0*sin(theta)."""
    from wpg.optical_elements import Mirror_plane
    from wpg.wpg_uti_oe import get_opd
    theta, length = 2e-3, 0.4
    fn = _bump_profile(tmp_path / 'h1d.dat')
    op0 = Mirror_plane('x', theta, length, 1e-3, fn, x0=0.)
    op1 = Mirror_plane('x', theta, length, 1e-3, fn, x0=0.1)
    xs = np.linspace(op0.mesh.xStart, op0.mesh.xFin, op0.mesh.nx)
    o0, o1 = get_opd(op0)[0], get_opd(op1)[0]
    shift = xs[np.argmax(np.abs(o1))] - xs[np.argmax(np.abs(o0))]
    step = xs[1] - xs[0]
    assert abs(abs(shift) - 0.1 * np.sin(theta)) <= step, f'shift {shift*1e6:.1f} um, expected 200 um'


def test_mirror_plane_delimiter(tmp_path):
    """6: delim is passed to loadtxt; whitespace delimiters keep reading files with several spaces/tabs."""
    from wpg.optical_elements import Mirror_plane
    from wpg.wpg_uti_oe import get_opd
    fn_c = _bump_profile(tmp_path / 'h1d_comma.dat', delimiter=',')
    fn_s = _bump_profile(tmp_path / 'h1d_space.dat')
    op_c = Mirror_plane('x', 2e-3, 0.4, 1e-3, fn_c, delim=',')
    op_s = Mirror_plane('x', 2e-3, 0.4, 1e-3, fn_s)
    np.testing.assert_array_equal(get_opd(op_c), get_opd(op_s))
    # sample profile of WPG: columns separated by ' \t'
    sample = os.path.join(os.path.dirname(__file__), '..', 'samples', 'data_common', 'mirror1.dat')
    if os.path.isfile(sample):
        op = Mirror_plane('x', 2e-3, 0.8, 1e-3, sample)
        assert np.abs(get_opd(op)).max() > 0


# ---------------------------------------------------------------- 7, 17, nz=1. generators

@pytest.mark.parametrize('builder', ['t', 'xy'])
def test_gauss_mode_orders(builder):
    """7: _my alone gives the HG01 mode (zero on axis), _mx alone does not crash."""
    from wpg.wavefront import Wavefront
    from wpg.generators import build_gauss_wavefront, build_gauss_wavefront_xy
    R = 5 * S0

    def build(**kw):
        if builder == 't':
            w = build_gauss_wavefront(41, 41, 3, EKEV, -R, R, -R, R, 5e-15, S0, S0, 0., **kw)
        else:
            w = build_gauss_wavefront_xy(41, 41, EKEV, -R, R, -R, R, S0, S0, 0., **kw)
        return Wavefront(w).get_intensity().sum(axis=-1)

    for kw in [dict(_my=1), dict(_mx=1)]:
        i = build(**kw)
        assert i[20, 20] / i.max() < 1e-3, f'{kw}: I(centre)/I_max = {i[20, 20] / i.max():.3f}'
    i = build(_my=1)
    assert i[20, :].max() / i.max() < 1e-3  # zero line along x (y = 0) for HG01


def test_gauss_frequency_units():
    """17: build_gauss_wavefront field in sqrt(W/mm^2); after time->frequency the spectrum integral is the pulse energy."""
    from srwpy.srwlib import srwl
    from wpg.wavefront import Wavefront
    from wpg.generators import build_gauss_wavefront
    from wpg.wpg_uti_wf import calc_pulse_energy
    R, pulseEn = 5 * S0, 1e-3
    expect = pulseEn * erf(4 / np.sqrt(2)) * erf(5 / np.sqrt(2)) ** 2
    srw_wf = build_gauss_wavefront(64, 64, 64, EKEV, -R, R, -R, R, 5e-15, S0, S0, 0., pulseEn=pulseEn)
    wf = Wavefront(srw_wf)
    assert wf.params.wEFieldUnit == 'sqrt(W/mm^2)'
    assert calc_pulse_energy(wf) == pytest.approx(expect, rel=1e-3)
    w = copy.deepcopy(srw_wf)
    srwl.SetRepresElecField(w, 'f')
    wf_f = Wavefront(w)
    m = wf_f.params.Mesh
    dx = (m.xMax - m.xMin) / (m.nx - 1)
    dy = (m.yMax - m.yMin) / (m.ny - 1)
    dE = (m.sliceMax - m.sliceMin) / (m.nSlices - 1)
    energy = wf_f.get_intensity().sum() * dx * dy * 1e6 * dE
    assert wf_f.params.wEFieldUnit == 'sqrt(J/eV/mm^2)'
    assert energy == pytest.approx(expect, rel=1e-3), f'spectrum integral {energy:.4e}, expected {expect:.4e} J'


def test_gauss_single_slice():
    """build_gauss_wavefront(nz=1): SRW left eFin = NaN and the field at t = -4 sigT;
    now a monochromatic wavefront (as build_gauss_wavefront_xy) that propagates correctly."""
    from srwpy import srwlib
    from srwpy.srwlib import srwl
    from wpg.wavefront import Wavefront
    from wpg.generators import build_gauss_wavefront, build_gauss_wavefront_xy
    R = 10 * S0
    with pytest.warns(UserWarning, match='monochromatic'):
        w = build_gauss_wavefront(128, 128, 1, EKEV, -R, R, -R, R, 5e-15, S0, S0, 0., _my=1)
    ref = build_gauss_wavefront_xy(128, 128, EKEV, -R, R, -R, R, S0, S0, 0., pulseEn=1e-3, pulseTau=5e-15, _my=1)
    assert (w.presFT, w.mesh.ne) == (0, 1)
    assert w.mesh.eStart == w.mesh.eFin == EKEV * 1e3
    np.testing.assert_array_equal(np.array(w.arEx), np.array(ref.arEx))
    with pytest.warns(UserWarning):
        w0 = build_gauss_wavefront(128, 128, 1, EKEV, -R, R, -R, R, 5e-15, S0, S0, 0.)
    srwl.PropagElecField(w0, srwlib.SRWLOptC([srwlib.SRWLOptD(2.)], [[0, 0, 1., 0, 0, 1., 1., 1., 1., 0, 0, 0]]))
    wf = Wavefront(w0)
    ix = wf.get_intensity().sum(axis=(0, 2))
    x = np.linspace(wf.params.Mesh.xMin, wf.params.Mesh.xMax, wf.params.Mesh.nx)
    zr = 4 * np.pi * S0 ** 2 / (12.398e-10 / EKEV)
    assert np.sqrt((ix * x ** 2).sum() / ix.sum()) == pytest.approx(S0 * np.sqrt(1 + (2 / zr) ** 2), rel=1e-3)


# ---------------------------------------------------------------- 8, 9. Use_PP

def test_use_pp_semi_analytical_modes():
    """8: all SRW modes of pp[3] are kept, wrong values raise."""
    from wpg.optical_elements import Use_PP
    for v in [0, 1, 2, 3, 4, 5, 100]:
        assert Use_PP(semi_analytical_treatment=v).pp[3] == v
    for v in [-1, 6, 2.5, 'a']:
        with pytest.raises(ValueError):
            Use_PP(semi_analytical_treatment=v)
    with pytest.raises(ValueError):
        Use_PP(auto_resize_before=2)


def test_use_pp_srw_pp():
    """9: srw_pp tuple accepted, wrong length raises, caller's list is not changed."""
    from wpg.optical_elements import Use_PP
    pp = [0, 0, 1.0, 0, 0, 2.0, 1.0, 2.0, 1.0, 0, 0, 0]
    assert Use_PP(srw_pp=tuple(pp)).pp[5] == 2.0
    with pytest.raises(ValueError):
        Use_PP(srw_pp=pp[:-1] + [0, 0, 0, 0, 0])  # 16 elements
    user = list(pp)
    p = Use_PP(srw_pp=user, zoom=3)
    assert p.pp[5] == 3.0 and user[5] == 2.0 and p.pp is not user
    assert len(Use_PP(srw_pp=pp + [0, 1, 0, 1, 0]).pp) == 17


# ---------------------------------------------------------------- 31. Empty

@pytest.mark.parametrize('pp_kind', ['use_pp', 'flat', 'nested'])
def test_empty_propagate_resizes(pp_kind):
    """31: Empty.propagate resizes the wavefront (SRW 4.2.1 rejects an empty SRWLOptC)."""
    from wpg.optical_elements import Empty, Use_PP
    wf = _gauss_xy(nx=128)
    r0 = wf.params.Mesh.xMax - wf.params.Mesh.xMin
    pp = Use_PP(zoom=2, sampling=0.5)
    arg = {'use_pp': pp, 'flat': pp.get_srw_pp(), 'nested': [pp.get_srw_pp()]}[pp_kind]
    Empty().propagate(wf, arg)
    r1 = wf.params.Mesh.xMax - wf.params.Mesh.xMin
    assert r1 / r0 == pytest.approx(2, rel=0.02)
    assert wf.params.Mesh.zCoord == 0.


# ---------------------------------------------------------------- 10, 12, 13. vis_tools

def test_vis_get_limits2d():
    """12: get_limits2d used the non-existing wf.srwl_wfr."""
    from wpg import vis_tools
    wf = _gauss_t(nx=20, ny=10, nz=3, R=50e-6, Ry=200e-6)
    lim = vis_tools.get_limits2d(wf, 'z')
    assert lim == pytest.approx([-50e-6, 50e-6, -200e-6, 200e-6])


def test_vis_show_2d_data_ndarray():
    """13: 2D ndarray input."""
    from wpg import vis_tools
    plt.figure()
    vis_tools.show_2d_data(func=None, data=np.ones((4, 5)))
    assert plt.gca().images[0].get_array().shape == (4, 5)
    plt.close('all')


@pytest.mark.parametrize('axis', ['x', 'y'])
def test_vis_show_2d_data_extent(axis):
    """10: extent of a x-t (y-t) map is the range of the axis left in the slice; no extra empty axes."""
    from wpg import vis_tools
    wf = _gauss_t(nx=40, ny=20, nz=9, R=50e-6, Ry=200e-6)
    m = wf.params.Mesh
    fig = plt.figure()
    vis_tools.show_intensity(wf, slice_number='center', axis=axis)
    assert len(fig.axes) == 1
    img = fig.axes[0].images[0]
    rows = img.get_array().shape[0]
    left, right, bottom, top = img.get_extent()
    if axis == 'x':  # data[iy, :, :]: rows along x
        assert rows == m.nx
        assert (top, bottom) == pytest.approx((m.xMin, m.xMax))
    else:
        assert rows == m.ny
        assert (top, bottom) == pytest.approx((m.yMin, m.yMax))
    assert (left, right) == pytest.approx((m.sliceMin, m.sliceMax))
    plt.close('all')


# ---------------------------------------------------------------- 11, 14, 16. wpg_uti_wf, wpg_uti_oe

def test_check_sampling_roi():
    """11: ROI test compares the window size (not Rx) with 3*FWHM; the y line says Vertical."""
    from wpg import wpg_uti_wf
    wf = _gauss_xy(nx=128, R=2 * S0, d2waist=50.)
    fw = wpg_uti_wf.calculate_fwhm(wf)
    rep = wpg_uti_wf.check_sampling(wf)
    roi_x = wf.params.Mesh.xMax - wf.params.Mesh.xMin
    assert roi_x < 3 * fw['fwhm_x'] < wf.params.Rx
    line_x = [ln for ln in rep.splitlines() if 'ROI' in ln and 'FWHM(x)' in ln][0]
    line_y = [ln for ln in rep.splitlines() if 'ROI' in ln and 'FWHM(y)' in ln][0]
    assert 'Increase ROI' in line_x, line_x
    assert line_y.startswith('Vertical'), line_y


def test_propagation_parameters_modes():
    """14: unknown mode raises ValueError, names are case-insensitive."""
    from wpg.wpg_uti_oe import propagation_parameters
    assert propagation_parameters(1, 1, 1, 1, mode='Fresnel')[3] == 0
    assert propagation_parameters(1, 1, 1, 1, mode='converge')[3] == 4
    with pytest.raises(ValueError):
        propagation_parameters(1, 1, 1, 1, mode='frenel')


def test_integral_intensity_single_slice():
    """16: one slice in time domain -- no division by zero."""
    from wpg import wpg_uti_wf
    wf = _gauss_xy(nx=32, pulseEn=1e-3, pulseTau=5e-15)
    wf.params.wDomain = 'time'
    wf.params.Mesh.sliceMin = wf.params.Mesh.sliceMax = 0.
    res = wpg_uti_wf.integral_intensity(wf, bPlot=False)
    assert np.isfinite(res)


# ---------------------------------------------------------------- minor findings

def test_averaged_intensity_returns_value():
    """19: averaged_intensity passed bPlot as threshold and returned None."""
    from wpg import wpg_uti_wf
    wf = _gauss_t(nx=32, ny=32, nz=9)
    res = wpg_uti_wf.averaged_intensity(wf)
    assert res == pytest.approx(wpg_uti_wf.integral_intensity(wf, bPlot=False))
    assert len(plt.get_fignums()) == 0


def test_integral_intensity_keeps_last_slice(capsys):
    """20: int0[min(aw):max(aw)] dropped the last meaningful slice."""
    from wpg import wpg_uti_wf
    wf = _gauss_t(nx=32, ny=32, nz=9)
    int0 = wf.get_intensity().sum(axis=(0, 1))
    n_expected = int((int0 > int0.max() * 0.01).sum())
    wpg_uti_wf.integral_intensity(wf, bPlot=False)
    assert f'number of meaningful slices: {n_expected}' in capsys.readouterr().out


def test_intensity_on_axis_vertical():
    """24: on-axis spectrum of a vertically polarized beam is not zero."""
    from wpg import wpg_uti_wf
    wf = _gauss_t(nx=33, ny=33, nz=9)
    wf.data.arrEver, wf.data.arrEhor = wf.data.arrEhor, np.zeros_like(wf.data.arrEhor)
    sz = wpg_uti_wf.get_intensity_on_axis(wf)
    assert sz[:, 1].max() == pytest.approx(1.0)


def test_fit_gaussian_without_p0():
    """25: fit_gaussian without p0 returned the start point (1, 1, 1) for data in metres."""
    from wpg.wpg_uti_math import fit_gaussian, gaussian
    x = np.linspace(-1e-4, 1e-4, 201)
    data = gaussian(x, 2e-5, 5e-6, 3e14)
    p, _ = fit_gaussian(x, data)
    assert p == pytest.approx([2e-5, 5e-6, 3e14], rel=1e-4)


def test_screen_existing_file_message(tmp_path):
    """28: the error message names the file."""
    from wpg.optical_elements import Screen
    fn = tmp_path / 'screen.h5'
    fn.write_bytes(b'')
    with pytest.raises(IOError, match='screen.h5'):
        Screen(str(fn))


# ---------------------------------------------------------------------------- issues #92 and #83

def test_crl_numpy_delta():
    """Issue #92: CRL with delta and attenuation length vs photon energy as numpy arrays failed with TypeError
    (SRW recognizes only lists); now the same element as with lists."""
    from wpg.optical_elements import CRL
    delta, att = np.array([4.80e-6, 4.79e-6, 4.78e-6]), np.array([6.0e-3, 6.1e-3, 6.2e-3])
    args = dict(_foc_plane=3, _shape=1, _apert_h=1e-3, _apert_v=1e-3, _r_min=5e-4, _n=10, _wall_thick=50e-6,
                _xc=0, _yc=0, _e_start=8400., _e_fin=8460., _nx=21, _ny=21)
    a = CRL(_delta=delta, _atten_len=att, **args)
    b = CRL(_delta=list(delta), _atten_len=list(att), **args)
    assert a.mesh.ne == 3 and np.array_equal(np.asarray(a.arTr), np.asarray(b.arTr))


def _gauss_9kev(d2waist):
    import wpg
    from wpg.generators import build_gauss_wavefront_xy
    return wpg.Wavefront(build_gauss_wavefront_xy(128, 128, 9., -2e-4, 2e-4, -2e-4, 2e-4, 20e-6, 20e-6, d2waist))


def _through(w, *elements):
    import wpg
    from wpg.optical_elements import Use_PP
    w = copy.deepcopy(w)
    bl = wpg.Beamline()
    for el in elements:
        bl.append(el, Use_PP())
    bl.propagate(w)
    return w


def test_transmission_amplitude_and_phase():
    """Issue #83: amplitude 0.5 and OPD lambda/8 multiply the field by 0.5 exp(+i pi/4)."""
    from wpg.optical_elements import Transmission
    lam = 1.239841984e-6 / 9e3
    w = _gauss_9kev(1.)
    t = Transmission(np.full((32, 32), 0.5), np.full((32, 32), lam / 8), 8e-4, 8e-4)
    e0 = w.get_real_part()[64, 64, 0] + 1j * w.get_imag_part()[64, 64, 0]
    w1 = _through(w, t)
    e1 = w1.get_real_part()[64, 64, 0] + 1j * w1.get_imag_part()[64, 64, 0]
    assert abs(e1 / e0) == pytest.approx(0.5, rel=1e-4)
    assert np.angle(e1 / e0) == pytest.approx(np.pi / 4, abs=1e-4)


def test_transmission_prism_deflects():
    """Issue #83: a linear OPD s*x deflects the beam by +s: after 20 m the centroid moves by 40 um for s = 2 urad."""
    from wpg.optical_elements import Drift, Transmission
    xs = np.linspace(-4e-4, 4e-4, 32)
    prism = Transmission(1., np.broadcast_to(2e-6 * xs[None, :], (32, 32)), 8e-4, 8e-4, ext_tr=1)
    w = _gauss_9kev(1.)

    def cx(w):
        m = w.params.Mesh
        p = w.get_intensity()[..., 0].sum(0)
        return (p * np.linspace(m.xMin, m.xMax, m.nx)).sum() / p.sum()

    shift = cx(_through(w, prism, Drift(20.))) - cx(_through(w, Drift(20.)))
    assert shift == pytest.approx(40e-6, rel=1e-2)


def test_transmission_wrong_shape():
    from wpg.optical_elements import Transmission
    with pytest.raises(ValueError):
        Transmission(np.ones(5), np.zeros(5), 1e-3, 1e-3)
    with pytest.raises(ValueError):
        Transmission(np.ones((4, 5)), np.zeros((5, 4)), 1e-3, 1e-3)
