# -*- coding: utf-8 -*-
"""Regression tests for the fixes of WPG_notes/review/core.md (items 1-7, 10-12, 14-18).

Test wavefront as in tests/test_wpg.py::_gauss_wavefront: 5 keV, 270 m from the waist,
RMS intensity size at the waist sigX = 12.0 um, RMS intensity divergence
theta = lam / (4 pi sigX) = 1.644 urad, RMS pulse duration of intensity 0.120 fs.
RMS intensity size after a drift z from the waist: sigma(z) = sqrt(sigX^2 + (theta z)^2).
"""
import warnings

import numpy as np
import pytest

EKEV = 5.0
D2WAIST = 270.
QNC = 0.1
THETA = (17.2 - 6.4 * np.sqrt(QNC)) * 1e-6 / EKEV ** 0.85 / 2.35  # 1.644e-6 rad
SIGX0 = 12.4e-10 / (EKEV * 4 * np.pi * THETA)
TAU = 0.12e-15


def _wf(nx=200, nz=3, ny=None, sigy_factor=1.0):
    import wpg
    from wpg.generators import build_gauss_wavefront
    ny = nx if ny is None else ny
    xmax = THETA * D2WAIST * 3.5
    ymax = xmax * ny / nx
    srw_wf = build_gauss_wavefront(nx, ny, nz, EKEV, -xmax, xmax, -ymax, ymax, TAU,
                                   SIGX0, SIGX0 * sigy_factor, D2WAIST)
    return wpg.Wavefront(srw_wf)


def _sigma(sig0, z):
    lam = 12.398419843320026e-10 / EKEV
    return np.sqrt(sig0 ** 2 + (lam / (4 * np.pi * sig0) * z) ** 2)


def _rms(axis, profile):
    p = np.asarray(profile, dtype=float)
    c = (axis * p).sum() / p.sum()
    return np.sqrt(((axis - c) ** 2 * p).sum() / p.sum())


def _state(wf):
    m = wf._srwl_wf.mesh
    return (wf.params.wSpace, wf.params.wDomain, m.xStart, m.xFin, m.yStart, m.yFin, m.eStart, m.eFin,
            np.array(wf._srwl_wf.arEx, dtype='float32').copy(), np.array(wf._srwl_wf.arEy, dtype='float32').copy())


def _assert_same_state(s0, s1):
    assert s0[:8] == s1[:8]
    np.testing.assert_array_equal(s0[8], s1[8])
    np.testing.assert_array_equal(s0[9], s1[9])


# ------------------------------------------------------------- 1. Q-space profile and divergence

def test_divergence_in_radians():
    """core.md 1: was (sigma, x0, a) tuples with sigma_qx = -9.10e3 (1/m, R-space profile)."""
    wf = _wf(nx=200, nz=3)
    sx, sy = wf.divergence
    assert isinstance(sx, float) and isinstance(sy, float)
    assert sx == pytest.approx(THETA, rel=0.01)
    assert sy == pytest.approx(THETA, rel=0.01)


def test_q_profiles_not_swapped():
    """core.md 1: get_profile('qx') was the R-space y-profile; sigY = 2 sigX -> theta_y = theta / 2."""
    wf = _wf(nx=200, nz=3, sigy_factor=2.0)
    rms_qx = _rms(wf.get_axis('qx'), wf.get_profile('qx'))
    rms_qy = _rms(wf.get_axis('qy'), wf.get_profile('qy'))
    assert rms_qx == pytest.approx(THETA, rel=0.01)
    assert rms_qy == pytest.approx(THETA / 2, rel=0.01)


# ------------------------------------------------------------- 2, 3. getters keep the representation

@pytest.mark.parametrize('call', [
    lambda wf: wf.get_axis('qx'),
    lambda wf: wf.get_axis('qy'),
    lambda wf: wf.get_pixel_size('qx'),
    lambda wf: wf.get_profile('qx'),
    lambda wf: wf.get_beam_width('qy'),
    lambda wf: wf.divergence,
    lambda wf: wf.get_axis('f'),
    lambda wf: wf.get_pixel_size('f'),
    lambda wf: wf.get_profile('f'),
])
def test_getters_do_not_change_wavefront(call):
    """core.md 2, 3: getters made FFTs of the user wavefront (mesh shift -0.5 px, |dI|/I = 1.5 %)
    and left it in another domain."""
    wf = _wf(nx=200, nz=3)
    s0 = _state(wf)
    call(wf)
    _assert_same_state(s0, _state(wf))


def test_getters_q_space_frequency_wavefront():
    """core.md 3: Q-space, frequency domain -> get_axis('t') left the wavefront in time domain."""
    wf = _wf(nx=64, nz=3)
    wf.set_electric_field_representation('a')
    wf.set_electric_field_representation('f')
    s0 = _state(wf)
    wf.get_axis('t')
    wf.get_profile('x')
    wf.get_pixel_size('y')
    _assert_same_state(s0, _state(wf))
    assert (wf.params.wSpace, wf.params.wDomain) == ('Q-space', 'frequency')


def test_odd_nx_getter_keeps_field_finite():
    """core.md 2: nx = 101 -> one get_axis('qx') gave 10239 NaN in the user field."""
    wf = _wf(nx=101, nz=3)
    s0 = _state(wf)
    wf.get_axis('qx')
    _assert_same_state(s0, _state(wf))
    assert np.isfinite(np.asarray(wf._srwl_wf.arEx)).all()


def test_time_profile_of_frequency_domain_wavefront():
    """core.md 3: get_profile('t') took intensity before the change of domain (returned the spectrum)."""
    wf_t = _wf(nx=40, nz=21)
    wf_f = _wf(nx=40, nz=21)
    wf_f.set_electric_field_representation('f')
    p_t = wf_t.get_profile('t')
    p_f = wf_f.get_profile('t')
    # shapes are compared: SRW 4.2.1 round trip time -> frequency -> time multiplies the intensity
    # by (nSlices / (nSlices - 1))**2 = 1.1025 for 21 slices
    np.testing.assert_allclose(p_f / p_f.max(), p_t / p_t.max(), rtol=0, atol=1e-3)
    np.testing.assert_allclose(wf_f.get_axis('t'), wf_t.get_axis('t'), rtol=1e-6)


# ------------------------------------------------------------- 4, 5. gaussian fit

def test_pulse_duration():
    """core.md 4: 41 slices, TAU = 0.120 fs -> pulse_duration gave sigma = -0.0386 fs (tuple)."""
    wf = _wf(nx=60, nz=41)
    tau = wf.pulse_duration
    assert isinstance(tau, float)
    assert tau == pytest.approx(TAU, rel=0.01)


def test_beam_width_positive():
    """core.md 5: get_beam_width('x') gave sigma = -4.4405e-4 m, analytic +4.4405e-4 m."""
    wf = _wf(nx=200, nz=3)
    sigma, x0, a = wf.get_beam_width('x')
    assert sigma == pytest.approx(_sigma(SIGX0, D2WAIST), rel=1e-3)
    assert abs(x0) < 1e-3 * sigma
    assert a == pytest.approx(wf.get_profile('x').max(), rel=0.01)


# ------------------------------------------------------------- 6. pixel size

def test_pixel_size_is_mesh_step():
    """core.md 6: step was (max - min)/n, 1.0 % (nx = 100) and 5.0 % (21 slices) smaller."""
    wf = _wf(nx=100, nz=21)
    for axis in ['x', 'y', 't', 'qx', 'f']:
        values = wf.get_axis(axis)
        assert wf.get_pixel_size(axis) == pytest.approx(values[1] - values[0], rel=1e-9), axis


# ------------------------------------------------------------- 7. profile through the centre, axes order

def test_profile_center_and_axes():
    """core.md 7: get_profile('x', method='center') returned an array (200, 3)."""
    wf = _wf(nx=64, ny=96, nz=5, sigy_factor=2.0)
    intensity = wf.get_intensity()
    ny, nx, ns = intensity.shape
    np.testing.assert_array_equal(wf.get_profile('x', method='center'), intensity[ny // 2, :, ns // 2])
    np.testing.assert_array_equal(wf.get_profile('y', method='center'), intensity[:, nx // 2, ns // 2])
    np.testing.assert_array_equal(wf.get_profile('t', method='center'), intensity[ny // 2, nx // 2, :])
    np.testing.assert_allclose(wf.get_profile('x', method='sum'), intensity.sum(axis=(0, 2)))
    assert wf.get_profile('x').shape == (nx,)
    assert wf.get_profile('y').shape == (ny,)
    assert wf.get_profile('qx', method='center').shape == (nx,)
    # x and y are not swapped: sigma_y at the waist is 2 times larger -> 2 times smaller at 270 m
    assert _rms(wf.get_axis('x'), wf.get_profile('x')) == pytest.approx(_sigma(SIGX0, D2WAIST), rel=0.01)
    assert _rms(wf.get_axis('y'), wf.get_profile('y')) == pytest.approx(_sigma(2 * SIGX0, D2WAIST), rel=0.01)


# ------------------------------------------------------------- 10, 11. polarization

def _wf_two_pol():
    wf = _wf(nx=64, nz=3)
    ey = 0.5j * wf.get_efield('horizontal')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        wf.data.arrEver = np.stack([ey.real, ey.imag], axis=-1).astype('float32')
    return wf, ey


def test_parts_default_horizontal_and_total_rejected():
    """core.md 10: 'total' (default) returned only Ex silently: 52 % off Re Ex + Re Ey."""
    wf, ey = _wf_two_pol()
    ex = wf.get_efield('horizontal')
    np.testing.assert_array_equal(wf.get_real_part(), wf.get_real_part(polarization='horizontal'))
    np.testing.assert_array_equal(wf.get_imag_part(), wf.get_imag_part(polarization='horizontal'))
    np.testing.assert_allclose(wf.get_real_part(polarization='vertical'), ey.real,
                               atol=1e-6 * np.abs(ey).max())
    ph = wf.get_phase()
    assert np.abs(np.angle(np.exp(1j * (ph - np.angle(ex))))).max() < 1e-5
    for f in (wf.get_real_part, wf.get_imag_part, wf.get_phase):
        with pytest.raises(ValueError, match='total'):
            f(polarization='total')


def test_efield_vertical():
    """core.md 11: get_efield('vertical') -> ValueError: could not broadcast (200,200,3) into (200,200,3,0)."""
    wf, ey = _wf_two_pol()
    ev = wf.get_efield('vertical')
    assert ev.shape == (64, 64, 3)
    np.testing.assert_allclose(ev, ey, atol=1e-6 * np.abs(ey).max())


# ------------------------------------------------------------- 12. beamline: groups after Empty

def _named_drift(name, length):
    from wpg.optical_elements import Drift
    d = Drift(length)
    d.name = name
    return d


def test_beamline_names_in_all_groups():
    """core.md 12: append(d1), append(Empty()), append(d2) -> index {'d1': 0}, KeyError 'd2'."""
    from wpg.beamline import Beamline
    from wpg.optical_elements import Empty, Use_PP
    b = Beamline()
    d1, d2, d3 = _named_drift('d1', 1.), _named_drift('d2', 2.), _named_drift('d3', 3.)
    b.append(d1, Use_PP())
    b.append(Empty(), Use_PP(zoom=2.))
    b.append(d2, Use_PP())
    b.append(d3, Use_PP())
    assert len(b.propagation_options) == 2
    assert b.index == {'d1': 0, 'd2': 1, 'd3': 2}

    b.edit_element_property('d2', 'L', 3.)
    assert d2.L == 3.
    assert b.get_element_properties('d2')['L'] == 3.
    new_pp = Use_PP(zoom=1.5)
    b.edit_propagation_parameters('d3', new_pp)
    assert b.propagation_options[1]['propagation_parameters'][1] == new_pp.get_srw_pp()

    d4 = _named_drift('d4', 4.)
    b.replace_element('d3', d4)
    assert b.propagation_options[1]['optical_elements'][1] is d4

    b.remove_element('d2')
    assert b.index == {'d1': 0, 'd4': 1}
    assert len(b.propagation_options[1]['propagation_parameters']) == 1

    b.crop_beamline('d1')
    assert len(b.propagation_options) == 1
    assert b.index == {'d1': 0}


def test_beamline_empty_first():
    """Empty as the first (or only) element: SRW 4.2.1 rejects SRWLOptC([], [pp])."""
    from wpg.beamline import Beamline
    from wpg.optical_elements import Empty, Use_PP
    wf = _wf(nx=64, nz=3)
    wf.set_electric_field_representation('f')
    r0, nx0, z0 = wf.params.Mesh.xMax - wf.params.Mesh.xMin, wf.params.Mesh.nx, wf.params.Mesh.zCoord
    b = Beamline()
    b.append(Empty(), Use_PP(zoom=2., sampling=0.5))
    b.propagate(wf)
    m = wf.params.Mesh
    assert (m.xMax - m.xMin) / r0 == pytest.approx(2., rel=0.03)
    assert m.nx == pytest.approx(nx0, abs=2)
    assert m.zCoord == z0


def test_beamline_crop_and_drift_after_empty():
    """Drift after Empty is propagated and its length is added to zCoord."""
    from wpg.beamline import Beamline
    from wpg.optical_elements import Empty, Use_PP
    wf = _wf(nx=200, nz=3)
    wf.set_electric_field_representation('f')
    b = Beamline()
    b.append(Empty(), Use_PP())
    b.append(_named_drift('d1', 5.), Use_PP())
    b.propagate(wf)
    assert wf.params.Mesh.zCoord == pytest.approx(D2WAIST + 5.)
    assert wf.get_beam_width('x')[0] == pytest.approx(_sigma(SIGX0, D2WAIST + 5.), rel=1e-3)


# ------------------------------------------------------------- 14-16. HDF5

def test_store_list_of_strings(tmp_path):
    """core.md 14: ['a', 'bb'] -> TypeError: No conversion path for dtype('<U2'); 16: str came back as bytes."""
    wf = _wf(nx=20, nz=3)
    wf.custom_fields['/misc/list_str'] = ['a', 'bb']
    wf.custom_fields['/misc/text'] = 'hello'
    fn = str(tmp_path / 'wf.h5')
    wf.store_hdf5(fn)
    import wpg
    wf2 = wpg.Wavefront()
    wf2.load_hdf5(fn)
    assert list(wf2.custom_fields['misc']['list_str']) == ['a', 'bb']
    assert wf2.custom_fields['misc']['text'] == 'hello'


def test_load_clears_custom_fields(tmp_path):
    """core.md 15: custom fields of the previous file stayed after load_hdf5."""
    import wpg
    wf = _wf(nx=20, nz=3)
    fn = str(tmp_path / 'wf.h5')
    wf.store_hdf5(fn)
    wf2 = wpg.Wavefront()
    wf2.custom_fields['/misc/old'] = 'from the previous file'
    wf2.load_hdf5(fn)
    assert '/misc/old' not in wf2.custom_fields
    assert 'misc' not in wf2.custom_fields


# ------------------------------------------------------------- 17. get_limits

def test_get_limits_silent(capsys):
    """core.md 17: get_limits() printed 'R-space' on every call."""
    wf = _wf(nx=20, nz=3)
    wf.get_limits()
    assert capsys.readouterr().out == ''

