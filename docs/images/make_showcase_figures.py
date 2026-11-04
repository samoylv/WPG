"""Figures for the README and the slides, computed with WPG.

Run from the repository root: python docs/images/make_showcase_figures.py [name ...]
(names: kb_focus, builds, sampling, reflectivity; all by default). Writes docs/images/*.png.

kb_focus: Gaussian beam 9 keV (RMS waist 5 um, 50 m upstream) focused by a Kirkpatrick-Baez pair of
elliptical mirrors (3 mrad, B4C coating), ideal and with random height errors of 1.5 nm RMS.
"""
import copy
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

import wpg  # noqa: E402
from wpg.generators import build_gauss_wavefront_xy  # noqa: E402
from wpg.optical_elements import Drift, Empty, Mirror_elliptical, Use_PP  # noqa: E402

warnings.simplefilter('ignore')
OUT = Path(__file__).parent
plt.rcParams.update({'font.size': 13, 'axes.titlesize': 14, 'figure.dpi': 110, 'savefig.bbox': 'tight'})
EKEV, TH = 9., 3e-3


def gauss(sig, d2waist, half, n=512):
    return wpg.Wavefront(build_gauss_wavefront_xy(nx=n, ny=n, ekev=EKEV, xMin=-half, xMax=half, yMin=-half,
                                                  yMax=half, sigX=sig, sigY=sig, d2waist=d2waist))


def intensity_map(w):
    m = w.params.Mesh
    return w.get_intensity()[..., 0], (m.xMin, m.xMax, m.yMin, m.yMax)


def kb_focus():
    from wpg.materials import Reflectivity
    from wpg.surface import Surface_errors, random_profile
    p_h, d_hv, f_v = 50., 0.6, 4.4          # source -> HKB, HKB -> VKB, VKB -> focus [m]
    angles = np.linspace(0.8 * TH, 1.2 * TH, 201)
    refl = Reflectivity.from_material('B4C', EKEV, angles, rho=2.52)

    def beamline(errors):
        bl = wpg.Beamline()
        semi = Use_PP(semi_analytical_treatment=1)
        bl.append(Mirror_elliptical('x', p_h, d_hv + f_v, TH, TH, 0.5, reflectivity=refl), semi)
        if errors:
            x, h = random_profile(0.5, 2001, 1.5e-9, -2.5, seed=1)
            bl.append(Surface_errors('x', TH, (x, h), width=1.2e-3), semi)
        bl.append(Drift(d_hv), semi)
        bl.append(Mirror_elliptical('y', p_h + d_hv, f_v, TH, TH, 0.5, reflectivity=refl), semi)
        if errors:
            x, h = random_profile(0.5, 2001, 1.5e-9, -2.5, seed=2)
            bl.append(Surface_errors('y', TH, (x, h), width=1.2e-3), semi)
        # propagation to the waist (SRW mode 4): the focal mesh comes from the Fourier transform
        bl.append(Drift(f_v), Use_PP(semi_analytical_treatment=4))
        bl.append(Empty(), Use_PP(zoom=0.4, sampling=1.))
        return bl

    # coarse pixel (4 per RMS size at the mirrors) and a wide window: the focal pixel is lambda*q/window = 0.05 um
    w0 = gauss(5e-6, p_h, 7e-3)
    results = {}
    for label, err in (('ideal KB mirrors', False), ('mirrors with 1.5 nm RMS height errors', True)):
        w = copy.deepcopy(w0)
        beamline(err).propagate(w)
        results[label] = intensity_map(w)
    I0, ext0 = intensity_map(gauss(5e-6, p_h, 5e-4, 256))  # the same beam on a fine mesh, for the picture only
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), layout='constrained')
    axes[0].imshow(I0 / I0.max(), extent=[1e3 * v for v in ext0], origin='lower', cmap='inferno')
    axes[0].set(title='(a) beam at the KB mirrors, 50 m (linear)', xlabel='x, mm', ylabel='y, mm')
    peak = max(v[0].max() for v in results.values())
    titles = ['(b) focus, ideal mirrors', '(c) focus, 1.5 nm RMS height errors']
    for ax, title, (I, ext) in zip(axes[1:], titles, results.values()):
        im = ax.imshow(np.log10(np.maximum(I / peak, 1e-5)), extent=[1e6 * v for v in ext], origin='lower',
                       cmap='inferno', vmin=-5, vmax=0)
        ax.set(title=title, xlabel='x, µm', ylabel='y, µm')
    cb = fig.colorbar(im, ax=axes[1:], shrink=0.9)
    cb.set_label('log₁₀ intensity, ideal focus peak = 0')
    fig.savefig(OUT / 'kb_focus.png')
    for label, (I, ext) in results.items():
        print(f'{label}: peak {I.max() / peak:.3f}, window {1e6 * (ext[1] - ext[0]):.2f} µm, n {I.shape}')


def builds():
    """A 5 m drift of a 5 keV pulse, 1024 x 1024 points, 128 slices, Use_PP() (no semi-analytical treatment),
    including the conversion to the frequency domain and back; Linux server with 32 cores and an RTX 3090,
    SRW master 3c38a0d, median of 3 runs (November 2026)."""
    labels = ['PyPI wheel\n(serial)', 'source build\n(serial)', 'OpenMP\n(32 cores)', 'CUDA\n(RTX 3090)']
    times = [27., 33., 9.6, 8.9]
    fig, ax = plt.subplots(figsize=(7.5, 4.2), layout='constrained')
    bars = ax.bar(labels, times, color=['#9e9e9e', '#9e9e9e', '#1f77b4', '#2ca02c'])
    for b, t in zip(bars, times):
        ax.text(b.get_x() + b.get_width() / 2, t + 0.6, f'{t:g} s', ha='center')
    ax.set(ylabel='time, s (median of 3 runs)', ylim=(0, 36),
           title='SRW builds: 5 m drift of a 1024 × 1024 × 128-slice pulse')
    ax.grid(True, axis='y', alpha=0.4)
    fig.savefig(OUT / 'srw_builds.png')


def reflectivity():
    """wpg.materials (xraylib) against the CXRO calculator (henke.lbl.gov, s-polarization, data/*.dat)."""
    from wpg import materials
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout='constrained')
    for name, formula, rho, ekev in (('B4C_9keV_s', 'B4C', 2.52, 9.), ('Pt_9keV_s', 'Pt', 21.45, 9.),
                                     ('Si_12keV_s', 'Si', 2.33, 12.)):
        deg, R = np.loadtxt(OUT / 'data' / f'{name}.dat', skiprows=2).T
        theta = np.deg2rad(deg)
        th = np.linspace(theta[0], theta[-1], 400)
        r_s, _ = materials.fresnel_reflectivity(materials.refractive_index(formula, ekev, rho), th)
        line, = ax.plot(1e3 * th, np.abs(r_s) ** 2, lw=2, label=f'{formula}, {ekev:g} keV, wpg.materials')
        ax.plot(1e3 * theta[::3], R[::3], 'o', mfc='none', color=line.get_color())
    ax.plot([], [], 'ko', mfc='none', label='CXRO calculator')
    ax.set(xlabel='grazing angle, mrad', ylabel='reflectivity', ylim=(0, 1.02),
           title='Mirror reflectivity: xraylib + Fresnel vs CXRO')
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=10, loc='center right')
    fig.savefig(OUT / 'mirror_reflectivity.png')


def sampling():
    """Gaussian beam 9 keV, waist 5 um 1 m upstream, window +-0.12 mm, drift 30 m without the
    semi-analytical treatment: Use_PP() against the parameters chosen by suggest_propagation_parameters."""
    from wpg.sampling import suggest_propagation_parameters
    lam = 1.239841984e-6 / (EKEV * 1e3)
    sig0, z0, L = 5e-6, 1., 30.
    w0 = gauss(sig0, z0, 1.2e-4, 256)
    pp, _ = suggest_propagation_parameters(w0, L, semi_analytical=0)
    zR = 4 * np.pi * sig0 ** 2 / lam
    sig = sig0 * np.sqrt(1 + ((z0 + L) / zR) ** 2)
    fig, ax = plt.subplots(figsize=(7.5, 4.6), layout='constrained')
    for label, use in (('Use_PP()', Use_PP(semi_analytical_treatment=0)),
                       (f'suggested: zoom {pp.zoom_h:.2f}', pp)):
        w = copy.deepcopy(w0)
        bl = wpg.Beamline()
        bl.append(Drift(L), use)
        bl.propagate(w)
        m = w.params.Mesh
        x = np.linspace(m.xMin, m.xMax, m.nx)
        p = w.get_intensity()[..., 0].sum(0)
        p = p / (p.sum() * (x[1] - x[0]))
        rms = np.sqrt((p * x ** 2).sum() * (x[1] - x[0]))
        ax.plot(1e6 * x, p * 1e-6, lw=2, label=f'{label}: RMS {1e6 * rms:.1f} µm')
    xs = np.linspace(-4 * sig, 4 * sig, 400)
    ax.plot(1e6 * xs, np.exp(-xs ** 2 / (2 * sig ** 2)) / (np.sqrt(2 * np.pi) * sig) * 1e-6, 'k--',
            label=f'analytic: RMS {1e6 * sig:.1f} µm')
    ax.set(xlabel='x, µm', ylabel='projected intensity, 1/µm',
           title='Drift of 30 m: the window chosen by suggest_propagation_parameters')
    ax.grid(True, alpha=0.4)
    ax.legend(fontsize=10, loc='upper left')
    fig.savefig(OUT / 'sampling_suggest.png')


if __name__ == '__main__':
    names = sys.argv[1:] or ['kb_focus', 'builds', 'reflectivity', 'sampling']
    for name in names:
        globals()[name]()
