# -*- coding: utf-8 -*-
"""Mirrors of WPG (wpg.optical_elements.Mirror_*) and reflectivity (wpg.materials).

Test geometry: monochromatic gaussian beam 9 keV, RMS intensity size at the waist 5 um, waist 50 m
before the mirror (RMS size at the mirror 109.7 um), grazing angle 3 mrad, horizontal deflection.
"""
from array import array

import numpy as np
import pytest

EKEV, P, Q, TH = 9., 50., 5., 3e-3


def _wf(sig=5e-6, d2waist=P, half=6e-4, nx=256, ny=256, polar='h'):
    import wpg
    from wpg.generators import build_gauss_wavefront_xy
    w = wpg.Wavefront(build_gauss_wavefront_xy(nx=nx, ny=ny, ekev=EKEV, xMin=-half, xMax=half, yMin=-half,
                                               yMax=half, sigX=sig, sigY=sig, d2waist=d2waist))
    if polar == 'v':
        s = w._srwl_wf
        s.arEy, s.arEx = s.arEx, array('f', [0.] * len(s.arEx))
    return w


def _power(w, comp=None):
    """Integral of the intensity over the cross-section: total or of one field component ('x' or 'y')."""
    s = w._srwl_wf
    arrs = {'x': [s.arEx], 'y': [s.arEy], None: [s.arEx, s.arEy]}[comp]
    m = w.params.Mesh
    return sum(float((np.asarray(a, dtype=float) ** 2).sum()) for a in arrs) * \
        (m.xMax - m.xMin) / (m.nx - 1) * (m.yMax - m.yMin) / (m.ny - 1)


def _rms_x(w):
    I = w.get_intensity()[..., 0].sum(0)
    x = np.linspace(w.params.Mesh.xMin, w.params.Mesh.xMax, w.params.Mesh.nx)
    p = I / I.sum()
    return np.sqrt((p * (x - (p * x).sum()) ** 2).sum())


def _propagate(w, mirror, d=None):
    import wpg
    from wpg.optical_elements import Drift, Use_PP
    bl = wpg.Beamline()
    bl.append(mirror, Use_PP(semi_analytical_treatment=1))
    if d is not None:
        bl.append(Drift(d), Use_PP(semi_analytical_treatment=1))
    bl.propagate(w)
    return w


# ---------------------------------------------------------------------------- orientation

def test_orientation_normalized_and_compatible():
    """The normal is a unit vector with roll; without roll it is the same as in WPG before 2026.11."""
    from wpg.optical_elements import _mirror_orientation
    for orient in 'xy':
        o = _mirror_orientation(orient, TH, roll=0.02)
        assert np.hypot(np.hypot(o['_nvx'], o['_nvy']), o['_nvz']) == pytest.approx(1., abs=1e-12)
    o = _mirror_orientation('x', TH)
    assert (o['_nvx'], o['_nvy'], o['_nvz'], o['_tvx'], o['_tvy']) == \
        pytest.approx((np.cos(TH), 0., -np.sin(TH), -np.sin(TH), 0.))
    o = _mirror_orientation('y', TH, yaw=1e-3, distance=10.)
    assert (o['_nvx'], o['_nvy'], o['_nvz'], o['_tvx'], o['_tvy'], o['_y']) == \
        pytest.approx((0., np.cos(TH), -np.sin(TH), 0., -np.sin(TH), np.tan(1e-3) * 10.))
    with pytest.raises(TypeError):
        _mirror_orientation('z', TH)


# ---------------------------------------------------------------------------- focusing

def test_focusing_mirrors_geometry():
    """Elliptical, toroidal and spherical mirrors (p = 50 m, q = 5 m) at 2.5 m: RMS size = geometric
    109.7 um * (1 - 2.5 / 5); convex hyperbolic with q = 5 m: 109.7 um * (1 + 2.5 / 5)."""
    from wpg.optical_elements import Mirror_elliptical, Mirror_hyperbolic, Mirror_spherical, Mirror_toroidal
    s_m = _rms_x(_wf())
    r_t = 2 * P * Q / ((P + Q) * np.sin(TH))
    focusing = [Mirror_elliptical('x', P, Q, TH, TH, 0.5), Mirror_toroidal('x', r_t, 1e40, TH, 0.5),
                Mirror_spherical('x', r_t, TH, 0.5)]
    for m in focusing:
        w = _propagate(_wf(), m, 2.5)
        assert _rms_x(w) == pytest.approx(s_m * 0.5, rel=2e-3)
    w = _propagate(_wf(), Mirror_hyperbolic('x', P, Q, TH, TH, 0.5, convex=True), 2.5)
    assert _rms_x(w) == pytest.approx(s_m * 1.5, rel=2e-3)


@pytest.mark.parametrize('collimating', [False, True])
def test_parabolic_mirror_focal_length(collimating):
    """Paraboloid f = 50 m: a source at 50 m is collimated (size kept at 50 m), a nearly parallel beam
    (waist 50 um 1000 m upstream, wavefront radius R = 1052 m) gives at 25 m the geometric size of a lens
    f = 50 m. SRW 'c' works only with the tangential vector opposite to its examples (_mirror_orientation)."""
    from wpg.optical_elements import Mirror_parabolic
    s_m = _rms_x(_wf())
    w = _propagate(_wf(), Mirror_parabolic('x', 50., TH, TH, 2., collimating=collimating), 50.)
    assert _rms_x(w) == pytest.approx(s_m, rel=2e-3)

    lam = 1.239841984e-6 / (EKEV * 1e3)
    sig0, z = 50e-6, 1000.
    R = z * (1 + (4 * np.pi * sig0 ** 2 / lam / z) ** 2)
    q = 1 / (1 / 50. - 1 / R)
    w0 = _wf(sig=sig0, d2waist=z, half=1.5e-3, nx=512)
    s_m = _rms_x(w0)
    w = _propagate(w0, Mirror_parabolic('x', 50., TH, TH, 2., collimating=collimating), 25.)
    assert _rms_x(w) == pytest.approx(s_m * (1 - 25. / q), rel=2e-3)


@pytest.mark.parametrize('collimating', [False, True])
def test_parabolic_mirror_short_keeps_energy(collimating):
    """Paraboloid 0.3 m long with a footprint RMS 37 mm: with the SRW default _npt = 500 it lost 9 %
    of the energy (limits of the surface), the wrapper sets _npt = 5000."""
    from wpg.optical_elements import Mirror_parabolic
    w = _wf(nx=128, ny=128)
    p0 = _power(w)
    _propagate(w, Mirror_parabolic('x', 50., TH, TH, 0.3, collimating=collimating), 0.01)
    assert _power(w) / p0 == pytest.approx(1., abs=1e-3)


# ---------------------------------------------------------------------------- reflectivity in SRW

def test_constant_reflectivity():
    """A number is the intensity reflectivity R (SRW multiplies the field by the amplitude sqrt(R))."""
    from wpg.optical_elements import Mirror_elliptical
    w = _wf(nx=128, ny=128)
    p0 = _power(w)
    _propagate(w, Mirror_elliptical('x', P, Q, TH, TH, 0.5, reflectivity=0.81))
    assert _power(w) / p0 == pytest.approx(0.81, abs=1e-3)
    with pytest.raises(ValueError):
        Mirror_elliptical('x', P, Q, TH, TH, 0.5, reflectivity=1.2)


@pytest.mark.parametrize('orient, sigma_comp', [('x', 'y'), ('y', 'x')])
def test_sigma_pi_components(orient, sigma_comp):
    """sigma is perpendicular to the plane of deflection: Ey for a horizontal mirror, Ex for a vertical one."""
    from wpg.materials import Reflectivity
    from wpg.optical_elements import Mirror_elliptical
    refl = Reflectivity(np.array([0.9, 0.6]).reshape(1, 1, 2), EKEV, TH)
    for polar, comp in [('h', 'x'), ('v', 'y')]:
        w = _wf(nx=128, ny=128, polar=polar)
        p0 = _power(w)
        _propagate(w, Mirror_elliptical(orient, P, Q, TH, TH, 0.5, reflectivity=refl))
        expected = 0.81 if comp == sigma_comp else 0.36
        assert _power(w, comp) / p0 == pytest.approx(expected, abs=1e-3)


def test_reflectivity_vs_energy_per_slice():
    """Table vs photon energy: every slice of a pulse in the frequency domain gets its own value."""
    from srwpy.srwlib import srwl
    from wpg.generators import build_gauss_wavefront
    from wpg.materials import Reflectivity
    from wpg.optical_elements import Mirror_elliptical
    import wpg
    w = wpg.Wavefront(build_gauss_wavefront(64, 64, 5, EKEV, -6e-4, 6e-4, -6e-4, 6e-4, 0.12e-15, 5e-6, 5e-6, P))
    srwl.SetRepresElecField(w._srwl_wf, 'f')
    m = w.params.Mesh
    ekev = np.linspace(m.sliceMin, m.sliceMax, m.nSlices) / 1e3
    R = np.linspace(0.3, 0.9, m.nSlices)
    I0 = w.get_intensity().sum(axis=(0, 1))
    _propagate(w, Mirror_elliptical('x', P, Q, TH, TH, 0.5, reflectivity=Reflectivity(np.sqrt(R)[:, None], ekev, TH)))
    np.testing.assert_allclose(w.get_intensity().sum(axis=(0, 1)) / I0, R, atol=2e-3)


# ---------------------------------------------------------------------------- materials (xraylib)

# CXRO mirror calculator (https://henke.lbl.gov/optical_constants/mirror2.html, s-polarization,
# fetched 04.11.2026, WPG_notes/mirrors/fetch_cxro_reflectivity.py): grazing angle 1..6 mrad
CXRO = {
    ('B4C', 2.52, 9., 0.): [0.999507, 0.998814, 0.997202, 0.117910, 0.02735247, 0.01057971],
    ('B4C', 2.52, 9., 1e-9): [0.999495, 0.998784, 0.997131, 0.110489, 0.02355941, 0.008290933],
    ('Pt', 21.45, 9., 0.): [0.981860, 0.963385, 0.943826, 0.922179, 0.896872, 0.864963],
    ('Si', 2.33, 12., 0.): [0.991195, 0.974801, 0.112359, 0.01867395, 0.00620622, 0.002707002],
}


@pytest.mark.parametrize('case', list(CXRO))
def test_fresnel_reflectivity_vs_cxro(case):
    """xraylib + Fresnel + Nevot-Croce agrees with CXRO (Henke tables) within 0.003 + 5 % (critical edge)."""
    pytest.importorskip('xraylib')
    from wpg import materials
    material, rho, ekev, sigma = case
    theta = np.arange(1, 7) * 1e-3
    r_s, r_p = materials.fresnel_reflectivity(materials.refractive_index(material, ekev, rho), theta, ekev, sigma)
    np.testing.assert_allclose(np.abs(r_s) ** 2, CXRO[case], rtol=0.05, atol=3e-3)
    np.testing.assert_allclose(np.abs(r_p) ** 2, np.abs(r_s) ** 2, rtol=1e-3)  # grazing incidence


def test_material_reflectivity_in_srw():
    """Pt mirror 9 keV at 5 mrad: the reflected energy equals |r|^2 at the mirror angle."""
    pytest.importorskip('xraylib')
    from wpg import materials
    from wpg.materials import Reflectivity
    from wpg.optical_elements import Mirror_elliptical
    th = 5e-3
    refl = Reflectivity.from_material('Pt', EKEV, np.linspace(0.8 * th, 1.2 * th, 201))
    w = _wf(nx=128, ny=128)
    p0 = _power(w)
    _propagate(w, Mirror_elliptical('x', P, Q, th, th, 0.5, reflectivity=refl))
    _, r_pi = materials.fresnel_reflectivity(materials.refractive_index('Pt', EKEV), th)
    assert _power(w) / p0 == pytest.approx(abs(r_pi) ** 2, abs=2e-3)  # Ex is pi for a horizontal mirror
    assert 0.85 < abs(r_pi) ** 2 < 0.92


def test_materials_helpers():
    pytest.importorskip('xraylib')
    from wpg import materials
    assert materials.density('Pt') == pytest.approx(21.45, rel=1e-2)
    with pytest.raises(ValueError, match='density'):
        materials.density('B4C')
    import xraylib
    # free-electron estimate far from absorption edges: delta = r_e lambda^2 n_e / (2 pi),
    # Be (Z = 4, A = 9.012, rho = 1.848 g/cm^3) at 8 keV: 5.32e-6
    lam = materials.HC_EV_M / 8e3
    n_e = 1.848 / 9.012 * 6.02214076e23 * 4 * 1e6
    delta, beta = materials.delta_beta('Be', 8., 1.848)
    assert delta == pytest.approx(2.8179403e-15 * lam ** 2 * n_e / (2 * np.pi), rel=0.01)
    # attenuation length from beta equals 1 / (rho * total cross-section) of xraylib (with scattering)
    mu = 1.848 * xraylib.CS_Total(4, 8.) * 100  # 1/m
    assert materials.attenuation_length('Be', 8., 1.848) == pytest.approx(1 / mu, rel=1e-3)


def test_reflectivity_table_validation():
    from wpg.materials import Reflectivity
    with pytest.raises(ValueError, match='constant step'):
        Reflectivity(np.ones((3, 1)), [8., 9., 9.5], TH)
    with pytest.raises(ValueError, match='shape'):
        Reflectivity(np.ones((2, 3)), [8., 9.], TH)
    r = Reflectivity(np.ones((2, 3, 2)), [8., 9.], [1e-3, 2e-3, 3e-3])
    assert r.intensity().shape == (2, 3, 2)


# ---------------------------------------------------------------------------- output axis

def _centre(w, axis):
    I = w.get_intensity()[..., 0]
    m = w.params.Mesh
    if axis == 'x':
        c, p = np.linspace(m.xMin, m.xMax, m.nx), I.sum(0)
    else:
        c, p = np.linspace(m.yMin, m.yMax, m.ny), I.sum(1)
    return float((c * p).sum() / p.sum())


@pytest.mark.parametrize('orient', ['x', 'y'])
def test_fixed_output_axis_pitch_misalignment(orient):
    """fixed_axis=True (default): theta0 = thetaE + delta deflects the beam by 2 delta off the design axis
    (centre shift 2 delta d at d = 2.5 m); fixed_axis=False (SRW default): the beam stays on the axis;
    without misalignment both give the same field."""
    from wpg.optical_elements import Mirror_elliptical
    d, delta = 2.5, 1e-6
    ref = _centre(_propagate(_wf(), Mirror_elliptical(orient, P, Q, TH, TH, 0.5), d), orient)
    for fixed, expected in [(True, 2 * delta * d), (False, 0.)]:
        m = Mirror_elliptical(orient, P, Q, TH, TH + delta, 0.5, fixed_axis=fixed)
        c = _centre(_propagate(_wf(), m, d), orient)
        assert c - ref == pytest.approx(expected, abs=0.02e-6)
    w_fixed = _propagate(_wf(nx=128, ny=128), Mirror_elliptical(orient, P, Q, TH, TH, 0.5), d)
    w_free = _propagate(_wf(nx=128, ny=128), Mirror_elliptical(orient, P, Q, TH, TH, 0.5, fixed_axis=False), d)
    np.testing.assert_allclose(w_fixed.get_intensity(), w_free.get_intensity(), rtol=0,
                               atol=1e-6 * w_free.get_intensity().max())


def test_fixed_output_axis_parameters():
    """The axis is added to the 12 propagation parameters at propagation (stored ones unchanged);
    17 parameters given by the user are kept."""
    import wpg
    from wpg.beamline import _with_output_axis
    from wpg.optical_elements import Mirror_elliptical, Use_PP
    m = Mirror_elliptical('x', P, Q, TH, TH, 0.5)
    pp = Use_PP().get_srw_pp()
    full = _with_output_axis(m, pp)
    assert len(full) == 17 and full[12:] == pytest.approx([np.sin(2 * TH), 0., np.cos(2 * TH), np.cos(2 * TH), 0.])
    user = list(pp) + [0., 0., 1., 1., 0.]
    assert _with_output_axis(m, user) is user
    assert _with_output_axis(Mirror_elliptical('x', P, Q, TH, TH, 0.5, fixed_axis=False), pp) is pp
    bl = wpg.Beamline()
    bl.append(m, Use_PP())
    bl.propagate(_wf(nx=64, ny=64))
    assert len(bl.propagation_options[0]['propagation_parameters'][0]) == 12
