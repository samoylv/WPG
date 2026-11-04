# -*- coding: utf-8 -*-
"""
Sampling of a wavefront and choice of propagation parameters (Use_PP zoom and sampling).

:func:`sampling_report` measures the sampling on the field itself: the phase step between
neighbouring pixels (Nyquist: below pi, here pi/2 is required for a margin), the radius of curvature
fitted to the phase of the field and the one SRW keeps for the semi-analytical treatment of the
quadratic phase (``params.Rx``, ``params.Ry``), the intensity at the window border.

:func:`suggest_propagation_parameters` chooses zoom and sampling for one optical element (usually a
drift) by trial propagations of one slice: a candidate is accepted when the result agrees with a run
on a finer and larger mesh. Sampling criteria of chirped fields in Fresnel propagation are discussed by
D. G. Voelz, M. C. Roggemann, Appl. Opt. 48, 6132 (2009), https://doi.org/10.1364/AO.48.006132.

.. module:: wpg.sampling
   :platform: Linux, Mac OSX, Windows
"""
import copy
import warnings
from array import array

import numpy as np
from scipy.integrate import trapezoid

#: h*c [eV*m]
HC_EV_M = 1.239841984e-6
#: intensity threshold (relative to the maximum) of the pixels taken into account
THRESHOLD = 1e-3
#: phase step [rad] from which the phase between neighbouring pixels is taken as wrapped (aliased)
ALIASED = 0.9 * np.pi


def _axes(wf):
    m = wf.params.Mesh
    return np.linspace(m.xMin, m.xMax, m.nx), np.linspace(m.yMin, m.yMax, m.ny)


def _field_2d(wf):
    """Complex field (ny, nx) of the polarization with more power, at the slice with the most power."""
    best = None
    for name in ('arrEhor', 'arrEver'):
        arr = getattr(wf.data, name)
        if arr is None:
            continue
        e = arr[..., 0].astype(float) + 1j * arr[..., 1].astype(float)  # (ny, nx, nz)
        power = (np.abs(e) ** 2).sum(axis=(0, 1))
        k = int(np.argmax(power))
        if best is None or power[k] > best[0]:
            best = (power[k], e[:, :, k], k)
    return best[1], best[2]


def _slice_energy_ev(wf, k):
    m = wf.params.Mesh
    if m.nSlices == 1:
        return wf.params.photonEnergy
    return np.linspace(m.sliceMin, m.sliceMax, m.nSlices)[k]


def _phase_steps(e, step, coords, axis, lam):
    """Phase steps between neighbours along an axis (0 - y, 1 - x) and a weighted fit of the
    phase gradient to a + b*coord: radius R = 2 pi / (lam * b) of exp(i pi coord^2 / (lam R))."""
    if axis == 1:
        p = e[:, 1:] * np.conj(e[:, :-1])
        mid = np.broadcast_to(0.5 * (coords[1:] + coords[:-1])[None, :], p.shape)
    else:
        p = e[1:, :] * np.conj(e[:-1, :])
        mid = np.broadcast_to(0.5 * (coords[1:] + coords[:-1])[:, None], p.shape)
    w = np.abs(p)
    mask = w > THRESHOLD * w.max()
    dphi = np.angle(p[mask])
    max_step = float(np.abs(dphi).max()) if dphi.size else 0.
    g, ww, xm = dphi / step, w[mask], mid[mask]
    if ww.sum() == 0:
        return max_step, np.inf
    xb = (ww * xm).sum() / ww.sum()
    var = (ww * (xm - xb) ** 2).sum()
    b = (ww * (xm - xb) * g).sum() / var if var > 0 else 0.
    radius = 2 * np.pi / (lam * b) if b != 0 else np.inf
    return max_step, radius


def sampling_report(wf):
    """
    Sampling of a wavefront measured on its field (the slice and polarization with the most power).

    :return: dict with keys 'x' and 'y', each a dict:
        'fwhm', 'rms' - beam size [m] (FWHM through the maximum, RMS of the projected intensity),
        'pixel', 'window' [m], 'edge' - maximum intensity at the window border relative to the maximum,
        'r_max' - largest distance from the centroid where the intensity exceeds THRESHOLD of the maximum [m],
        'phase_step' - largest phase step between neighbouring pixels where the intensity exceeds
        THRESHOLD [rad] (pi/2 or more: the field is undersampled; the value saturates at pi),
        'phase_step_residual' - the same after removing the quadratic phase with the radius SRW keeps
        for the semi-analytical treatment (params.Rx / Ry; None if it is 0),
        'R_params', 'R_fit' - radius of curvature [m]: SRW parameter and the fit to the phase of the field
        (nan if the phase step reaches ALIASED: the phase is wrapped, the fit is meaningless),
        'semi_ok' - the quadratic phase with the SRW radius describes the field (residual phase step
        below pi/2), i.e. ``Use_PP(semi_analytical_treatment=1)`` can be used; None if the SRW radius is 0,
        'pixels_per_fwhm'.
    """
    e, k = _field_2d(wf)
    lam = HC_EV_M / _slice_energy_ev(wf, k)
    x, y = _axes(wf)
    intensity = np.abs(e) ** 2
    iy, ix = np.unravel_index(np.argmax(intensity), intensity.shape)
    rep = {}
    for name, coords, axis, prof, cut, r_params, centre in (
            ('x', x, 1, intensity.sum(0), intensity[iy, :], wf._srwl_wf.Rx, wf._srwl_wf.xc),
            ('y', y, 0, intensity.sum(1), intensity[:, ix], wf._srwl_wf.Ry, wf._srwl_wf.yc)):
        step = (coords[-1] - coords[0]) / (coords.size - 1)
        p = prof / prof.sum()
        c = (p * coords).sum()
        above = coords[prof > THRESHOLD * prof.max()]
        phase_step, r_fit = _phase_steps(e, step, coords, axis, lam)
        if phase_step >= ALIASED:
            r_fit = np.nan  # wrapped phase: the fit is meaningless
        residual = None
        if r_params:
            q = np.exp(-1j * np.pi * (coords - centre) ** 2 / (lam * r_params))
            e_res = e * (q[None, :] if axis == 1 else q[:, None])
            residual = _phase_steps(e_res, step, coords, axis, lam)[0]
        fwhm = np.count_nonzero(cut > cut.max() / 2) * step
        rep[name] = dict(
            fwhm=fwhm, rms=float(np.sqrt((p * (coords - c) ** 2).sum())), pixel=step,
            window=coords[-1] - coords[0], r_max=float(np.abs(above - c).max()) if above.size else 0.,
            edge=float(max(prof[0], prof[-1]) / prof.max()),
            phase_step=phase_step, phase_step_residual=residual, R_params=r_params, R_fit=r_fit,
            semi_ok=None if residual is None else bool(residual < np.pi / 2),
            pixels_per_fwhm=fwhm / step)
    return rep


def check_sampling(wavefront):
    """
    Report on the sampling of a wavefront (text). See :func:`sampling_report` for the quantities.

    Checks: the window holds the beam (intensity at the border below 1e-3 of the maximum, window
    at least 3 FWHM); the phase step between neighbouring pixels is below pi/2 (Nyquist pi with a
    margin), for the field as it is and after removing the quadratic phase with the SRW radius
    (what the propagation with ``Use_PP(semi_analytical_treatment=1)`` works with). If the field is well
    sampled but not after removing that phase, the SRW radius disagrees with the field and the
    semi-analytical treatment fails: e.g. SRW gives Rx = 0.01 m for a Gaussian beam generated 1 cm from
    its waist (the true radius is ~1e8 m), and a drift of 50 m gave an RMS size of 0.5 m instead of 100 um.
    """
    rep = sampling_report(wavefront)
    fmt = '|{:9.3e}|{:9.3e}|{:9.3e}|{:9.3e}|{:10.3e}|{:10.3e}|{:9.2f}|'
    ret = 'WAVEFRONT SAMPLING REPORT\n'
    line = '+----------+---------+---------+---------+---------+----------+----------+---------+\n'
    ret += line
    ret += '|x/y       |FWHM [m] |px [m]   |ROI [m]  |edge I   |R SRW [m] |R fit [m] |dphi/pi  |\n'
    ret += line
    for label, a in (('Horizontal', rep['x']), ('Vertical  ', rep['y'])):
        ret += label + fmt.format(a['fwhm'], a['pixel'], a['window'], a['edge'], a['R_params'],
                                  a['R_fit'], a['phase_step'] / np.pi) + '\n'
    ret += line + '\n'
    for label, ax, a in (('Horizontal', 'x', rep['x']), ('Vertical', 'y', rep['y'])):
        width = 'width' if ax == 'x' else 'height'
        if a['window'] >= 3 * a['fwhm']:
            ret += '{} ROI > 3* FWHM({}) -> OK\n'.format(label, ax)
        else:
            ret += '{} ROI !> 3* FWHM({}) -> Increase ROI {} ({}).\n'.format(label, ax, width, ax)
        if a['edge'] <= THRESHOLD:
            ret += '{} intensity at the window border {:.1e} of max -> OK\n'.format(label, a['edge'])
        else:
            ret += '{} intensity at the window border {:.1e} of max -> the beam is cut, increase ROI ({}).\n'.format(
                label, a['edge'], ax)
        if a['phase_step'] < np.pi / 2:
            ret += '{} phase step {:.2f} pi < pi/2 -> OK for propagation without semi-analytical treatment\n'.format(
                label, a['phase_step'] / np.pi)
        else:
            ret += ('{} phase step {:.2f} pi >= pi/2 -> undersampled without semi-analytical treatment '
                    '(decrease pixel ({}) or use semi_analytical_treatment=1)\n').format(label, a['phase_step'] / np.pi, ax)
        if a['semi_ok'] is True:
            ret += ('{} phase step after removing the quadratic phase with the SRW radius {:.2f} pi -> OK with '
                    'semi_analytical_treatment=1\n').format(label, a['phase_step_residual'] / np.pi)
        elif a['semi_ok'] is False and a['phase_step'] < np.pi / 2:
            ret += ('{} SRW radius R = {:.3g} m disagrees with the curvature of the field R = {:.3g} m (phase step after '
                    'removing it {:.2f} pi) -> do not use semi_analytical_treatment=1, or set the radius of wf._srwl_wf '
                    'to the fitted one\n').format(label, a['R_params'], a['R_fit'], a['phase_step_residual'] / np.pi)
        elif a['semi_ok'] is False:
            ret += ('{} phase step after removing the quadratic phase with the SRW radius {:.2f} pi -> undersampled '
                    'also with semi_analytical_treatment=1, decrease pixel ({})\n').format(
                label, a['phase_step_residual'] / np.pi, ax)
        if a['pixels_per_fwhm'] > 100:
            ret += '{} {:.0f} pixels per FWHM: oversampled, sampling ({}) could be reduced\n'.format(
                label, a['pixels_per_fwhm'], ax)
    ret += '\nEND OF REPORT'
    return ret


def set_radius_from_field(wf):
    """
    Set the radius of curvature kept by SRW (wf._srwl_wf.Rx, Ry) to the one fitted to the phase of the
    field, along the axes where the SRW radius does not describe the field ('semi_ok' False) and the
    phase is sampled well enough for the fit (phase step below ALIASED).

    SRW uses this radius for the semi-analytical treatment and for resizing: with Rx = 0.01 m that SRW
    gives a Gaussian beam generated 1 cm from its waist (true radius ~8e7 m), a drift of 50 m with
    sampling 1.5 gave an RMS size of 502 um instead of 100 um even without the semi-analytical treatment.

    :return: dict {axis: (old radius, new radius)} of the changed axes
    """
    rep = sampling_report(wf)
    changed = {}
    for ax, attr in (('x', 'Rx'), ('y', 'Ry')):
        a = rep[ax]
        if a['semi_ok'] is False and np.isfinite(a['R_fit']):
            changed[ax] = (getattr(wf._srwl_wf, attr), float(a['R_fit']))
            setattr(wf._srwl_wf, attr, float(a['R_fit']))
    return changed


# ---------------------------------------------------------------------------- choice of Use_PP

def _single_slice(wf, k=None):
    """Copy of the wavefront with one slice in the frequency domain (the one with the most power)."""
    from srwpy.srwlib import srwl
    w = copy.deepcopy(wf)
    if w.params.wDomain == 'time':
        srwl.SetRepresElecField(w._srwl_wf, 'f')
    s = w._srwl_wf
    m = s.mesh
    if m.ne > 1:
        if k is None:
            _, k = _field_2d(w)
        energy = np.linspace(m.eStart, m.eFin, m.ne)[k]
        for name in ('arEx', 'arEy'):
            arr = np.asarray(getattr(s, name), dtype=np.float32).reshape(m.ny, m.nx, m.ne, 2)
            setattr(s, name, array('f', arr[:, :, k, :].ravel()))
        m.ne = 1
        m.eStart = m.eFin = energy
    return w


def _run(w, element, semi, zh, zv, sh, sv):
    import wpg
    from wpg.optical_elements import Use_PP
    w = copy.deepcopy(w)
    bl = wpg.Beamline()
    bl.append(element, Use_PP(semi_analytical_treatment=semi, zoom_h=zh, zoom_v=zv, sampling_h=sh, sampling_v=sv))
    bl.propagate(w)
    intensity = w.get_intensity()[..., 0]
    x, y = _axes(w)
    out = {}
    for name, coords, prof in (('x', x, intensity.sum(0)), ('y', y, intensity.sum(1))):
        step = (coords[-1] - coords[0]) / (coords.size - 1)
        density = prof / (prof.sum() * step)
        c = (density * coords).sum() * step
        out[name] = dict(coords=coords, density=density, rms=float(np.sqrt((density * (coords - c) ** 2).sum() * step)),
                         edge=float(max(prof[0], prof[-1]) / prof.max()), n=coords.size)
    return out


def _distance(a, b):
    """Total variation distance of two normalized profiles on different meshes (0 - same, 1 - disjoint)."""
    lo, hi = max(a['coords'][0], b['coords'][0]), min(a['coords'][-1], b['coords'][-1])
    grid = np.linspace(lo, hi, 4 * max(a['n'], b['n']))
    pa = np.interp(grid, a['coords'], a['density'])
    pb = np.interp(grid, b['coords'], b['density'])
    return 0.5 * trapezoid(np.abs(pa - pb), grid)


def suggest_propagation_parameters(wf, element, semi_analytical=None, tol=1e-2, factor=1.5, max_iter=8,
                                   max_points=8192):
    """
    Choose Use_PP zoom and sampling for one optical element by trial propagations of one slice.

    A candidate (zoom_h, zoom_v, sampling_h, sampling_v), starting from 1, is propagated and compared
    with a run on a mesh larger and finer by ``factor``. It is accepted when, along both axes, the RMS
    sizes agree within ``tol``, the projected profiles differ by less than ``tol`` (total variation
    distance) and the intensity at the border of the output window is below 1e-3 of the maximum.
    Otherwise zoom is increased along the axis where the beam reaches the border, sampling where it
    does not. Only zoom >= 1 is tried: zoom of SRW changes the window before the element, so a smaller
    window would cut the incident beam; crop the output with ``Empty`` and ``Use_PP(zoom=...)`` if needed.

    :param wf: wavefront (not changed); for a pulse the slice with the most power is used
    :param element: SRW or WPG optical element, a number - a drift of this length [m]
    :param semi_analytical: 0 or 1 (Use_PP semi_analytical_treatment); None - 1 if the field has a
        strong quadratic phase (phase step >= pi/4) and the SRW radius describes it (see
        :func:`sampling_report`, 'semi_ok'), else 0
    :param tol: relative accuracy of the RMS size and of the profile
    :param factor: step of zoom and sampling
    :param max_iter: maximum number of candidates
    :param max_points: maximum number of points along an axis after the element
    :return: (Use_PP, report dict: 'accepted', 'semi_analytical_treatment', 'history' of the candidates
        with the sizes and the reasons, 'input' - sampling_report of the slice, 'radius_fixed' - the SRW
        radius replaced by the fitted one on the trial copy, see :func:`set_radius_from_field`)
    """
    from wpg.optical_elements import Drift, Use_PP
    if np.isscalar(element):
        element = Drift(float(element))
    w1 = _single_slice(wf)
    radius_fixed = set_radius_from_field(w1)
    if radius_fixed:
        warnings.warn('The radius of curvature kept by SRW does not describe the field ({}): resizing and the '
                      'semi-analytical treatment fail with it. The parameters are chosen for the radius fitted to '
                      'the field: call wpg.sampling.set_radius_from_field(wf) before the propagation'.format(
                          ', '.join('{}: {:.3g} -> {:.3g} m'.format(k, *v) for k, v in radius_fixed.items())))
    rep0 = sampling_report(w1)
    semi_ok = all(rep0[a]['semi_ok'] is not False for a in 'xy')
    if semi_analytical is None:
        strong = any(rep0[a]['phase_step'] >= np.pi / 4 for a in 'xy')
        semi_analytical = int(strong and semi_ok)
    elif semi_analytical and not semi_ok:
        warnings.warn('The radius of curvature kept by SRW does not describe the field: the semi-analytical '
                      'treatment may fail (see check_sampling)')
    z = {'x': 1., 'y': 1.}
    s = {'x': 1., 'y': 1.}
    history = []
    accepted = False
    tested = (dict(z), dict(s))
    for _ in range(max_iter):
        tested = (dict(z), dict(s))
        a = _run(w1, element, semi_analytical, z['x'], z['y'], s['x'], s['y'])
        b = _run(w1, element, semi_analytical, z['x'] * factor, z['y'] * factor, s['x'] * factor, s['y'] * factor)
        step: dict = dict(zoom=dict(z), sampling=dict(s))
        ok_all = True
        for ax in 'xy':
            d_rms = abs(a[ax]['rms'] - b[ax]['rms']) / b[ax]['rms']
            d_prof = _distance(a[ax], b[ax])
            edge = a[ax]['edge']
            action = None
            if edge > THRESHOLD:
                z[ax] *= factor
                ok_all = False
                action = 'zoom'
            elif d_rms > tol or d_prof > tol:
                s[ax] *= factor
                ok_all = False
                action = 'sampling'
            step[ax] = dict(rms=a[ax]['rms'], rms_refined=b[ax]['rms'], d_rms=d_rms, d_profile=d_prof, edge=edge,
                            n=a[ax]['n'], action=action)
        history.append(step)
        if ok_all:
            accepted = True
            break
        if max(a['x']['n'], a['y']['n']) * factor ** 2 > max_points:
            warnings.warn('suggest_propagation_parameters: max_points reached before the accuracy {}'.format(tol))
            break
    if not accepted:
        warnings.warn('suggest_propagation_parameters: no candidate reached the accuracy {} in {} steps, '
                      'the last tested one is returned'.format(tol, len(history)))
    z, s = tested
    pp = Use_PP(semi_analytical_treatment=semi_analytical, zoom_h=z['x'], zoom_v=z['y'],
                sampling_h=s['x'], sampling_v=s['y'])
    return pp, dict(accepted=accepted, semi_analytical_treatment=semi_analytical, history=history,
                    input=rep0, radius_fixed=radius_fixed)
