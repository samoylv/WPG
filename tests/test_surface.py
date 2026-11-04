# -*- coding: utf-8 -*-
"""Surface height errors (wpg.surface): agreement with SRW, physics of the phase, random surfaces.

Beam for propagation tests: 9 keV, grazing angle 3 mrad, horizontal deflection; checks are described
in WPG_notes/surface/*.py.
"""
import numpy as np
import pytest

TH, EKEV = 3e-3, 9.


def _tr(opt):
    a = np.asarray(opt.arTr)
    return a[0::2], a[1::2]


# ---------------------------------------------------------------------------- same model as SRW

@pytest.mark.parametrize('orient, theta_r', [('x', None), ('y', None), ('x', 4e-3)])
def test_fill_opd_1d_equals_srw(orient, theta_r):
    """Vectorized 1D profile = srwl_opt_setup_surf_height_1d (profile off the mesh centre)."""
    from srwpy import srwlib
    from wpg.surface import fill_opd_1d
    rng = np.random.default_rng(1)
    x = np.linspace(-0.3, 0.5, 401)
    h = 1e-10 * rng.normal(size=x.size).cumsum()
    nx, ny = (600, 20) if orient == 'x' else (20, 600)
    sx, sy = (0.7 * TH, 1e-3) if orient == 'x' else (1e-3, 0.7 * TH)
    ref = srwlib.srwl_opt_setup_surf_height_1d(np.array([x, h]), orient, TH, _ang_r=theta_r or 0, _amp_coef=2.,
                                               _nx=nx, _ny=ny, _size_x=sx, _size_y=sy)
    new = fill_opd_1d(srwlib.SRWLOptT(nx, ny, sx, sy), x, h, orient, TH, theta_r, scale=2.)
    (a_r, o_r), (a_n, o_n) = _tr(ref), _tr(new)
    np.testing.assert_array_equal(a_n, a_r)
    np.testing.assert_allclose(o_n, o_r, rtol=0, atol=1e-12 * np.abs(o_r).max())


@pytest.mark.parametrize('orient', ['x', 'y'])
def test_fill_opd_2d_equals_srw(orient):
    """Vectorized 2D map = srwl_opt_setup_surf_height_2d (bilinear interpolation)."""
    from srwpy import srwlib
    from wpg.surface import fill_opd_2d
    rng = np.random.default_rng(2)
    xl, yt = np.linspace(-0.4, 0.4, 41), np.linspace(-0.01, 0.01, 11)
    h = 1e-9 * rng.normal(size=(yt.size, xl.size))
    tab = np.zeros((yt.size + 1, xl.size + 1))
    tab[0, 1:], tab[1:, 0], tab[1:, 1:] = xl, yt, h
    nx, ny = (120, 30) if orient == 'x' else (30, 120)
    sx, sy = (0.6 * TH, 0.015) if orient == 'x' else (0.015, 0.6 * TH)
    ref = srwlib.srwl_opt_setup_surf_height_2d(tab, orient, TH, _nx=nx, _ny=ny, _size_x=sx, _size_y=sy)
    new = fill_opd_2d(srwlib.SRWLOptT(nx, ny, sx, sy), xl, yt, h, orient, TH)
    np.testing.assert_allclose(_tr(new)[1], _tr(ref)[1], rtol=0, atol=1e-12 * np.abs(_tr(ref)[1]).max())


@pytest.mark.parametrize('orient', ['x', 'y'])
def test_aux_transm_vectorized(orient):
    """srwutils.AuxTransmAddSurfHeightProfileScaled (used by Mirror_plane and the tutorials) gives the
    result of SRW with the same angle of incidence and reflection; a profile with decreasing
    positions gives the same (the former loops gave zeros)."""
    from srwpy import srwlib
    from wpg.useful_code.srwutils import AuxTransmAddSurfHeightProfileScaled
    x = np.linspace(-0.4, 0.4, 801)
    h = 1e-9 * np.sin(2 * np.pi * x / 0.13)
    nx, ny = (900, 10) if orient == 'x' else (10, 900)
    sx, sy = (0.8 * TH, 1e-3) if orient == 'x' else (1e-3, 0.8 * TH)
    ref = srwlib.srwl_opt_setup_surf_height_1d(np.array([x, h]), orient, TH, _amp_coef=1.5, _nx=nx, _ny=ny,
                                               _size_x=sx, _size_y=sy)
    for prof in (np.array([x, h]), np.array([x[::-1], h[::-1]])):
        new = srwlib.SRWLOptT(nx, ny, sx, sy)
        AuxTransmAddSurfHeightProfileScaled(new, prof, orient, TH, 1.5)
        np.testing.assert_allclose(_tr(new)[1], _tr(ref)[1], rtol=0, atol=1e-12 * np.abs(_tr(ref)[1]).max())


# ---------------------------------------------------------------------------- physics

def _gauss(sig, d2waist, half, nx, ny=64, xoff=0.):
    import wpg
    from wpg.generators import build_gauss_wavefront_xy
    return wpg.Wavefront(build_gauss_wavefront_xy(nx=nx, ny=ny, ekev=EKEV, xMin=-half, xMax=half, yMin=-half,
                                                  yMax=half, sigX=sig, sigY=sig, d2waist=d2waist, xoff=xoff))


def _centre_rms(w):
    I = w.get_intensity()[..., 0].sum(0)
    x = np.linspace(w.params.Mesh.xMin, w.params.Mesh.xMax, w.params.Mesh.nx)
    p = I / I.sum()
    c = (p * x).sum()
    return c, np.sqrt((p * (x - c) ** 2).sum())


def _run(w, elements, d=None):
    import wpg
    from wpg.optical_elements import Drift, Use_PP
    bl = wpg.Beamline()
    for el in elements:
        bl.append(el, Use_PP(semi_analytical_treatment=1))
    if d is not None:
        bl.append(Drift(d), Use_PP(semi_analytical_treatment=1))
    bl.propagate(w)
    return w


def test_sine_profile_diffraction_orders():
    """h = a sin(2 pi x / 10 mm): phase grating, the share of the order m is J_m(phi0)^2 with
    phi0 = 2 k a sin(theta) = 1 rad, the angle of the order m lambda / (Lambda sin(theta))."""
    from scipy.special import jv
    from srwpy.srwlib import srwl
    from wpg.surface import Surface_errors
    lam = 1.239841984e-6 / (EKEV * 1e3)
    period = 1e-2
    a = 1. / (2 * 2 * np.pi / lam * np.sin(TH))
    x = np.linspace(-0.4, 0.4, 8001)
    w = _run(_gauss(1e-4, 1e-3, 8e-4, 1024), [Surface_errors('x', TH, (x, a * np.sin(2 * np.pi * x / period)),
                                                             width=3e-3, after_mirror=False)])
    srwl.SetRepresElecField(w._srwl_wf, 'a')
    m = w._srwl_wf.mesh
    ang = np.linspace(m.xStart, m.xFin, m.nx) * lam
    p = w.get_intensity()[..., 0].sum(0)
    p /= p.sum()
    step = lam / (period * np.sin(TH))
    for order in (-2, -1, 0, 1, 2):
        assert p[np.abs(ang - order * step) < 0.3 * step].sum() == pytest.approx(jv(order, 1.) ** 2, abs=2e-3)


def test_profile_sign_after_srw_mirror():
    """After an elliptical mirror (p = 50 m, q = 5 m): a concave profile c x^2 adds the focal power
    4c / sin(theta) (f = 50 m: size at 2.5 m 49.38 um instead of 54.87 um), a profile rising
    downstream s x increases the deflection by 2s (centre +5 um at 2.5 m for s = 1 urad)."""
    from wpg.optical_elements import Mirror_elliptical
    from wpg.surface import Surface_errors

    def wf():
        return _gauss(5e-6, 50., 6e-4, 256, 256)

    x = np.linspace(-0.25, 0.25, 2001)
    mirror = Mirror_elliptical('x', 50., 5., TH, TH, 0.5)
    c0, s_m = _centre_rms(wf())
    q = 1 / (1 / 5. + 1 / 50.)
    c = np.sin(TH) / (4 * 50.)
    _, s = _centre_rms(_run(wf(), [mirror, Surface_errors('x', TH, (x, c * x ** 2), width=1.2e-3)], 2.5))
    assert s == pytest.approx(s_m * (1 - 2.5 / q), rel=2e-3)
    c_ref, _ = _centre_rms(_run(wf(), [Mirror_elliptical('x', 50., 5., TH, TH, 0.5)], 2.5))
    c_t, _ = _centre_rms(_run(wf(), [Mirror_elliptical('x', 50., 5., TH, TH, 0.5),
                                     Surface_errors('x', TH, (x, 1e-6 * x), width=1.2e-3)], 2.5))
    assert c_t - c_ref == pytest.approx(2 * 1e-6 * 2.5, rel=1e-2)


def test_srw_mirror_inverts_image():
    """The output frame of an SRW mirror is inverted in the plane of deflection: a beam +50 um off
    axis is at -50 um after the mirror (the reason for Surface_errors(after_mirror=True))."""
    from wpg.optical_elements import Mirror_elliptical
    w = _gauss(5e-6, 50., 6e-4, 256, 256, xoff=5e-5)
    assert _centre_rms(w)[0] == pytest.approx(5e-5, rel=1e-3)
    c, _ = _centre_rms(_run(w, [Mirror_elliptical('x', 50., 1e6, TH, TH, 0.5)], 0.01))
    assert c == pytest.approx(-5e-5, rel=1e-2)


# ---------------------------------------------------------------------------- random surfaces

@pytest.mark.parametrize('exponent', [-1.5, -2.5])
def test_random_profile_statistics(exponent):
    """RMS height exactly sigma, mean 0, slope of the mean periodogram (100 realizations) = exponent."""
    from wpg.surface import random_profile
    L, n = 0.8, 2001
    prof = np.array([random_profile(L, n, 2e-9, exponent, seed=s)[1] for s in range(100)])
    np.testing.assert_allclose(prof.std(axis=1), 2e-9, rtol=1e-10)
    assert np.abs(prof.mean(axis=1)).max() < 1e-20
    dx = L / (n - 1)
    f = np.fft.rfftfreq(n, dx)
    p = (np.abs(np.fft.rfft(prof, axis=1)) ** 2).mean(0)
    sel = (f > 5 / L) & (f < f.max() / 5)
    assert np.polyfit(np.log(f[sel]), np.log(p[sel]), 1)[0] == pytest.approx(exponent, abs=0.05)
    np.testing.assert_array_equal(random_profile(L, n, 2e-9, exponent, seed=3)[1],
                                  random_profile(L, n, 2e-9, exponent, seed=3)[1])


def test_random_surface_and_element():
    """2D map: RMS sigma, shape (n_tr, n_long); Surface_errors of it is finite and of the given size."""
    from wpg.surface import Surface_errors, random_surface
    x, y, h = random_surface(0.8, 0.02, 401, 21, 1e-9, seed=0)
    assert h.shape == (21, 401)
    assert h.std() == pytest.approx(1e-9, rel=1e-10)
    opt = Surface_errors('y', TH, (x, y, h), width=0.015)
    assert (opt.mesh.nx, opt.mesh.ny) == (21, 401)
    assert opt.mesh.yFin - opt.mesh.yStart == pytest.approx(0.8 * np.sin(TH))
    amp, opd = _tr(opt)
    assert np.all(amp == 1.) and np.all(np.isfinite(opd)) and np.abs(opd).max() > 0
