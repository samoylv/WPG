# -*- coding: utf-8 -*-
"""Sampling report and choice of propagation parameters (wpg.sampling, wpg_uti_wf.check_sampling).

Gaussian beams 9 keV, mesh 256 x 256: RMS intensity size after a drift sigma(z) = sigma0 sqrt(1 + (z/zR)^2),
zR = 4 pi sigma0^2 / lambda; radius of curvature R(z) = z (1 + (zR/z)^2). Checks: WPG_notes/sampling/*.py.
"""
import copy
import warnings

import numpy as np
import pytest

EKEV = 9.
LAM = 1.239841984e-6 / (EKEV * 1e3)


def _gauss(sig, d2waist, half, n=256, xoff=0.):
    import wpg
    from wpg.generators import build_gauss_wavefront_xy
    return wpg.Wavefront(build_gauss_wavefront_xy(nx=n, ny=n, ekev=EKEV, xMin=-half, xMax=half, yMin=-half,
                                                  yMax=half, sigX=sig, sigY=sig, d2waist=d2waist, xoff=xoff))


def _gauss_srw_radius(sig, d2waist, half, n=256):
    """Gaussian beam with the radius that SRW alone keeps, Rx = Ry = distance to the waist (WPG before
    2026.11, other sources): the field is curved as R(z), the parameter is wrong near the waist."""
    w = _gauss(sig, d2waist, half, n)
    w._srwl_wf.Rx = w._srwl_wf.Ry = d2waist
    return w


def _drift(w, L, pp):
    import wpg
    from wpg.optical_elements import Drift
    w = copy.deepcopy(w)
    bl = wpg.Beamline()
    bl.append(Drift(L), pp)
    bl.propagate(w)
    m = w.params.Mesh
    p = w.get_intensity()[..., 0].sum(0)
    x = np.linspace(m.xMin, m.xMax, m.nx)
    p = p / p.sum()
    return np.sqrt((p * (x - (p * x).sum()) ** 2).sum())


def _sigma(sig0, z):
    return sig0 * np.sqrt(1 + (z / (4 * np.pi * sig0 ** 2 / LAM)) ** 2)


def _after_mirror():
    import wpg
    from wpg.optical_elements import Drift, Mirror_elliptical, Use_PP
    w = _gauss(5e-6, 50., 6e-4)
    bl = wpg.Beamline()
    bl.append(Mirror_elliptical('x', 50., 5., 3e-3, 3e-3, 0.5), Use_PP(semi_analytical_treatment=1))
    bl.append(Drift(0.01), Use_PP(semi_analytical_treatment=1))
    bl.propagate(w)
    return w


# ---------------------------------------------------------------------------- report

def test_report_diverging_beam():
    """R = 50 m: the fitted radius is the Gaussian one z (1 + (zR/z)^2) = 50.1 m; the phase step 0.55 pi
    asks for the semi-analytical treatment, after removing the quadratic phase it is ~0."""
    from wpg.sampling import sampling_report
    from wpg.wpg_uti_wf import check_sampling
    rep = sampling_report(_gauss(5e-6, 50., 6e-4))
    zR = 4 * np.pi * 5e-6 ** 2 / LAM
    for ax in 'xy':
        assert rep[ax]['R_fit'] == pytest.approx(50. * (1 + (zR / 50.) ** 2), rel=2e-3)
        assert np.pi / 2 < rep[ax]['phase_step'] < 0.6 * np.pi
        assert rep[ax]['phase_step_residual'] < 0.05 * np.pi and rep[ax]['semi_ok']
        assert rep[ax]['edge'] < 1e-3
    text = check_sampling(_gauss(5e-6, 50., 6e-4))
    assert 'OK with semi_analytical_treatment=1' in text and 'disagrees' not in text


def test_report_converging_beam():
    """After an elliptical mirror (R = -5 m): the phase is wrapped (no fit), the SRW radius describes the
    field - the semi-analytical treatment works (the former check never passed for R < 0)."""
    from wpg.sampling import sampling_report
    rep = sampling_report(_after_mirror())
    assert rep['x']['R_params'] == pytest.approx(-5., rel=1e-2)
    assert np.isnan(rep['x']['R_fit']) and rep['x']['phase_step'] > 0.9 * np.pi
    assert rep['x']['semi_ok']


def test_report_wrong_srw_radius():
    """Waist 1 mm away: SRW keeps Rx = 1 mm, the field is nearly flat (R = 8.3e8 m): reported as a disagreement."""
    from wpg.sampling import sampling_report
    from wpg.wpg_uti_wf import check_sampling
    w = _gauss_srw_radius(1e-4, 1e-3, 6e-4)
    rep = sampling_report(w)
    zR = 4 * np.pi * 1e-4 ** 2 / LAM
    assert rep['x']['R_params'] == pytest.approx(1e-3)
    assert rep['x']['R_fit'] == pytest.approx(1e-3 * (1 + (zR / 1e-3) ** 2), rel=0.05)
    assert rep['x']['phase_step'] < 0.01 * np.pi and rep['x']['semi_ok'] is False
    assert 'disagrees with the curvature of the field' in check_sampling(w)


def test_calculate_fwhm_off_axis():
    """FWHM through the maximum: an off-axis beam has the same width as a centred one."""
    from wpg.wpg_uti_wf import calculate_fwhm
    c = calculate_fwhm(_gauss(5e-6, 50., 6e-4))
    off = calculate_fwhm(_gauss(5e-6, 50., 6e-4, xoff=2.5e-4))
    assert off['fwhm_x'] == pytest.approx(c['fwhm_x'], rel=0.02)
    assert off['fwhm_y'] == pytest.approx(c['fwhm_y'], rel=0.02)


# ---------------------------------------------------------------------------- radius and resizing

def test_set_radius_from_field_fixes_resizing():
    """With the SRW radius 0.01 m of a beam 1 cm from its waist, sampling 1.5 breaks a drift of 50 m even
    without the semi-analytical treatment; with the fitted radius the size is the analytic one."""
    from wpg.optical_elements import Use_PP
    from wpg.sampling import set_radius_from_field
    pp = Use_PP(semi_analytical_treatment=0, sampling=1.5)
    expect = _sigma(1e-4, 50.01)
    w = _gauss_srw_radius(1e-4, 0.01, 8e-4)
    assert _drift(w, 50., pp) > 3 * expect
    changed = set_radius_from_field(w)
    assert set(changed) == {'x', 'y'} and changed['x'][0] == pytest.approx(0.01)
    assert _drift(w, 50., pp) == pytest.approx(expect, rel=1e-3)


# ---------------------------------------------------------------------------- suggest

@pytest.mark.parametrize('sig0, d2waist, half, L, semi', [
    (5e-6, 50., 6e-4, 50., 0),   # Use_PP(): -1.0 %
    (5e-6, 1., 1.2e-4, 30., 0),  # Use_PP(): -11.5 %, the beam reaches the window border
    (5e-6, 50., 6e-4, 50., None),
])
def test_suggest_gaussian_drift(sig0, d2waist, half, L, semi):
    from wpg.optical_elements import Use_PP
    from wpg.sampling import suggest_propagation_parameters
    w = _gauss(sig0, d2waist, half)
    expect = _sigma(sig0, d2waist + L)
    pp, rep = suggest_propagation_parameters(w, L, semi_analytical=semi)
    assert rep['accepted']
    assert _drift(w, L, pp) == pytest.approx(expect, rel=2e-3)
    if semi == 0:
        assert abs(_drift(w, L, Use_PP(semi_analytical_treatment=0)) / expect - 1) > 5e-3
    else:
        assert rep['semi_analytical_treatment'] == 1


def test_suggest_wrong_srw_radius_warns():
    """A beam 1 cm from its waist: the function warns about the SRW radius and chooses no semi-analytical
    treatment; the wavefront of the user is not changed."""
    from wpg.sampling import suggest_propagation_parameters
    w = _gauss_srw_radius(1e-4, 0.01, 8e-4)
    with pytest.warns(UserWarning, match='set_radius_from_field'):
        pp, rep = suggest_propagation_parameters(w, 50.)
    assert rep['semi_analytical_treatment'] == 0 and set(rep['radius_fixed']) == {'x', 'y'}
    assert w._srwl_wf.Rx == pytest.approx(0.01)


def test_single_slice_of_pulse():
    """A pulse in the time domain: the trial uses one slice in the frequency domain, the one with most power."""
    import wpg
    from wpg.generators import build_gauss_wavefront
    from wpg.sampling import _single_slice
    w = wpg.Wavefront(build_gauss_wavefront(32, 32, 9, EKEV, -6e-4, 6e-4, -6e-4, 6e-4, 0.12e-15, 5e-6, 5e-6, 50.))
    s = _single_slice(w)
    assert s.params.Mesh.nSlices == 1 and s.params.wDomain == 'frequency'
    assert s.params.Mesh.sliceMin == pytest.approx(EKEV * 1e3, rel=1e-3)
    assert w.params.Mesh.nSlices == 9 and w.params.wDomain == 'time'
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        assert np.isfinite(s.get_intensity()).all()


# ---------------------------------------------------------------------------- radius of build_gauss_wavefront

@pytest.mark.parametrize('d2waist', [1e-3, 1., 50., -1.])
def test_gauss_radius_exact(d2waist):
    """build_gauss_wavefront_xy keeps the radius of the field R = z (1 + (zR/z)^2), not z (issue #133)."""
    from wpg.sampling import sampling_report
    w = _gauss(2e-5, d2waist, 2e-4)
    zR = 4 * np.pi * 2e-5 ** 2 / LAM
    for attr in ('Rx', 'Ry'):
        assert getattr(w._srwl_wf, attr) == pytest.approx(d2waist * (1 + (zR / d2waist) ** 2), rel=1e-9)
    rep = sampling_report(w)
    assert rep['x']['semi_ok'] and rep['y']['semi_ok']


def test_drift_near_waist():
    """9 keV, RMS 20 um (zR = 36.5 m), 1 m from the waist, drift 100 m with zoom 4 and sampling 0.5: with
    Rx = 1 m of SRW the RMS size was 788 um, the analytic one is 59 um."""
    from wpg.optical_elements import Use_PP
    pp = Use_PP(semi_analytical_treatment=1, zoom=4, sampling=0.5)
    expect = _sigma(2e-5, 101.)
    assert _drift(_gauss(2e-5, 1., 2e-4), 100., pp) == pytest.approx(expect, rel=1e-2)
    assert _drift(_gauss_srw_radius(2e-5, 1., 2e-4), 100., pp) > 5 * expect


@pytest.mark.parametrize('nz', [9, 10])
def test_pulse_in_time_domain(nz):
    """build_gauss_wavefront gives a pulse in the time domain: the wavelength of the report is the one of
    photonEnergy, not of the slice time (with nz odd the central slice at t = 0 gave R_fit = 0)."""
    import wpg
    from wpg.generators import build_gauss_wavefront
    from wpg.sampling import sampling_report, set_radius_from_field
    w = wpg.Wavefront(build_gauss_wavefront(64, 64, nz, EKEV, -2e-4, 2e-4, -2e-4, 2e-4, 0.12e-15, 2e-5, 2e-5, 1.))
    assert w.params.wDomain == 'time'
    exact = 1. * (1 + (4 * np.pi * 2e-5 ** 2 / LAM) ** 2)
    rep = sampling_report(w)
    assert rep['x']['R_params'] == pytest.approx(exact, rel=1e-9) and rep['x']['semi_ok']
    w._srwl_wf.Rx = w._srwl_wf.Ry = 1.
    changed = set_radius_from_field(w)
    assert set(changed) == {'x', 'y'}
    assert changed['x'][1] == pytest.approx(exact, rel=1e-2)
