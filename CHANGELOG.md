# Changelog

The format follows [Keep a Changelog 1.1.0](https://keepachangelog.com/en/1.1.0/). Versions are named by year and month.

## [Unreleased] - planned as 2026.11

Changes from the `develop` branch (commit 70baa30) to the `feature/wpg_update` branch.
SRW is no longer bundled with WPG: it comes from the `srwpy` package.
WPG now runs on current Python, NumPy, h5py and Matplotlib.

### Upgrade notes

Several functions gave wrong numbers and are fixed, so results of existing scripts can change:

- `Wavefront.divergence` now returns rad. For a 5 keV Gaussian beam with a 12 µm RMS waist it gives the expected 1.644 µrad. It gave −9.1·10³ (spatial frequency in 1/m with a wrong axis).
- `Wavefront.pulse_duration` gives 0.120 fs for a 0.120 fs RMS test pulse. It gave −0.039 fs.
- `get_pixel_size` is (max − min)/(n − 1): 1 % larger than before at n = 100.
- `build_gauss_wavefront` sets the field units to √(W/mm²), so values of a pulse in the frequency domain change (see Changed).
- A mirror pitch error now deflects the beam (see Changed). Mirrors without misalignment give the same result.
- `build_gauss_wavefront` and `build_gauss_wavefront_xy` keep the radius of curvature of the Gaussian beam instead of the distance to the waist (see Fixed). Propagation with the semi-analytical treatment or resizing changes for beams generated within a few Rayleigh lengths of the waist.

### Added

- Builds of SRW with OpenMP or CUDA on Linux: `make install SRW_MODE=omp|cuda` (SRW and WPG), `make srw SRW_MODE=cpu|omp|cuda|wheel [SRW_REF=...]` (only SRW, `wheel` returns the official wheel) or `build/build_srw.py`. `make install` without `SRW_MODE` keeps the installed SRW and installs the wheel only if `srwpy` is missing. The script clones `ochubar/SRW` (any branch, tag or fork) or uses a local checkout. It installs `srwpy` into the active environment and checks that the library is linked with OpenMP or cuFFT as requested.
  - A 5 m drift of a 5 keV pulse of 1024 × 1024 points and 128 slices (without the semi-analytical treatment, including the conversion to the frequency domain and back) on a Linux server with 32 cores and an RTX 3090, median of 3 runs: PyPI wheel 27 s, local serial build 33 s, OpenMP on 32 cores 9.6 s, CUDA on the GPU 8.9 s.
  - OpenMP and CUDA cannot be combined, because SRW uses FFTW2 for OpenMP and FFTW3 with cuFFT for CUDA.
  - Windows uses the PyPI wheel.
- Propagation on GPU with a CUDA build: `Beamline.propagate(wf, gpu=1)` and `propagate_wavefront(..., gpu=1)`, where `gpu` is the device number. See Known issues.
- Mirror reflectivity (`wpg.materials`, optional dependency `xraylib`: `pip install wpg[materials]`).
  - `Reflectivity.from_material('B4C', ekev, angles, rho=2.52)` (density in g/cm³) gives the Fresnel reflectivity of a thick coating with Névot–Croce roughness, for σ and π polarizations. It agrees with the CXRO calculator within 0.007 (absolute, in R) for B4C, Pt and Si at 9–12 keV and 1–6 mrad.
  - `reflectivity=` of the mirrors accepts a constant intensity reflectivity R or a `Reflectivity` table R(photon energy, grazing angle).
  - Helpers: `refractive_index`, `delta_beta`, `attenuation_length`, `density`.
- Mirror shapes `Mirror_hyperbolic`, `Mirror_parabolic` (focusing or collimating), `Mirror_spherical` and `Mirror_toroidal`, with the orientation arguments of `Mirror_elliptical`.
- Surface height errors (`wpg.surface`).
  - `Surface_errors(orient, theta, profile, width)` builds a transmission element from a 1D profile or a 2D map of heights, in the model of SRW (`srwl_opt_setup_surf_height_1d/2d`). It changes only the phase, so it is placed right after the mirror whose surface it describes.
  - `random_profile` and `random_surface` generate Gaussian random surfaces with a power-law spectrum and a given RMS height.
- Sampling tools (`wpg.sampling`).
  - `sampling_report` and `check_sampling` measure the phase step between neighbouring pixels (with and without the quadratic phase that the semi-analytical treatment removes) and the intensity at the window border.
  - `set_radius_from_field` sets the radius of curvature that SRW keeps from the phase of the field, for wavefronts whose radius does not describe the field (for example a Gaussian beam built by SRW directly, see Fixed).
  - `suggest_propagation_parameters(wf, element)` chooses `Use_PP` zoom and sampling by trial propagations of one slice. In two test drifts of Gaussian beams (9 keV, 50 m and 30 m) the default `Use_PP()` gave RMS sizes 1.0 % and 11.5 % below the analytic ones. With the chosen parameters the difference is below 0.01 %.
- `Transmission(amplitude, opd, rx, ry)` builds a transmission element from two arrays, the amplitude transmission and the optical path difference (issue #83). The field after it is multiplied by amplitude·exp(i·k·opd).
- `wpg_uti_wf.plot_slices(wf, threshold)` plots the intensity maps of the slices whose power is above the threshold, in mm for a wavefront in R-space and in µrad in Q-space, with the time or the photon energy of each slice in the title. It comes from the slice plots of L. Samoylova and computes the intensity once without copying the wavefront: 198 maps of a SACLA pulse of 201 × 201 points and 2662 slices took 40 s and 3.6 GB of memory.
- `Mirror_plane` accepts `x0` (the profile point placed on the beam axis) and any whitespace delimiter in the profile file.
- Diagnostics for bug reports: `python -m wpg_diagnose` (or the `wpg-diagnose` command) prints the Python environment, the version and location of WPG, the origin of SRW with the command that repeats its build, OpenMP and CUDA in SRW, the GPU and the versions of the dependencies, as Markdown to paste into an issue (`--json` for JSON). Among the problems it lists an activated environment that is not the running Python, WPG or SRW installed in the user site-packages (`~/.local`) and an OpenMP runtime loaded from another installation. It runs also when `import wpg` fails. `--test` adds a self-test of about 10 s: a drift of a Gaussian beam against the analytic size, OpenMP with one thread against all threads, GPU against CPU.
  - `build/build_srw.py` saves the build description with the SHA-256 checksum of `srwlpy` next to the library, so the diagnostics recognize a description left from another build.
  - `import wpg` without SRW says how to install `srwpy` and how to report the problem.
- About 160 unit tests (`pytest`), run by GitHub Actions on Linux and Windows, and the execution of 7 sample notebooks (`pytest -m notebooks`, run by hand: the notebooks take several minutes and are not run by GitHub Actions).

### Changed

- SRW comes from the `srwpy` package: the PyPI wheel 4.2.1 or a local build (see Added). The SRW copy of 2019 in `wpg/srw` is removed.
- Python 3.10 or newer is required. Tested with Python 3.10 (Linux) and 3.12 (Windows, NumPy 2.5, h5py 3.16, Matplotlib 3.11). GitHub Actions are configured for Python 3.10 and 3.13.
- Packaging is in `pyproject.toml` (setuptools 77 or newer, the license as the SPDX expression `LicenseRef-WPG` with the file `LICENSE`). Install with `pip install -e ".[test]" --config-settings editable_mode=compat`: the compat mode of the setuptools editable install lets IDEs and static analyzers find the package.
- Field names of the wavefront parameters are declared in one module, `wpg/glossary.py`, with their types and units. Autocompletion and type checkers recognize `wf.params.Mesh.nx` and the like, and a misspelled field name raises `AttributeError` instead of creating a new attribute. The shortcuts `wf.nx`, `wf.Rx` and so on are read-only properties: assignment like `wf.nx = 100` raises an error, set `wf.params.Mesh.nx` instead. The HDF5 format is unchanged, files of the old version load identically.
- `build_gauss_wavefront`
  - The field is in √(W/mm²) (SRW field-units flag `unitElFld = 2`), so the spectrum of a pulse integrates to its energy. Previously SRW took the field as photons per second per 0.1 % bandwidth, and for a 1 mJ pulse the spectrum integrated to 6.2·10¹² instead of 10⁻³ J.
  - With `nz = 1` it builds a monochromatic wavefront and issues a warning, because SRW left the final photon energy `eFin` undefined (NaN) for a single time slice.
- `Use_PP` raises `ValueError` for a wrong switch value, an unknown `semi_analytical_treatment` mode or a malformed `srw_pp`, instead of keeping them silently.
- `Beamline`
  - `append` raises `TypeError` for an element that is neither an SRW optical element nor `Empty` or `Screen`. Previously such an element was dropped together with its propagation parameters.
  - `Screen` in a beamline writes the wavefront to its file after its resizing. Previously only the resizing was applied and no file was written.
  - `propagate` warns for a wavefront in the time domain. SRW propagates every slice as if its time were a photon energy, and the result is wrong: a 5 keV Gaussian pulse of 20 slices had an RMS size of 875 µm after a 5 m drift instead of the expected 452 µm. Convert to the frequency domain first with `srwl.SetRepresElecField(wf._srwl_wf, 'f')`.
- A mirror pitch error deflects the beam. `Mirror_elliptical`, `Mirror_hyperbolic` and `Mirror_parabolic` take the design grazing angle `thetaE` and the actual angle `theta0`. With `fixed_axis=True` (the default) the output optical axis is the reflected direction for `thetaE`, so a misalignment `theta0 != thetaE` deflects the beam by 2·(`theta0` − `thetaE`) off this axis, as in a real beamline. SRW alone builds the output axis from the actual mirror normal: the beam stayed on the axis and only the aberrations changed. `fixed_axis=False` restores that behaviour.
- `check_sampling` reports the sampling measured on the field (see Fixed), so the text of the report is different.
- `import wpg` no longer silences all warnings of the process, which also hid warnings of user code. Only the informational warnings of the glossary are silenced.
- Docker: one `Docker/Dockerfile` replaces the three images based on miniconda of 2019. Build arguments: `SRW_MODE` (wheel, cpu, omp or cuda) and `BASE` (the base image, e.g. a CUDA image for `SRW_MODE=cuda`).
- Documentation moved from Read the Docs (not rebuilt since 2019) to GitHub Pages, https://samoylv.github.io/WPG/. GitHub Actions build it on every push and pull request and publish it from `develop`.

### Deprecated

- `wpg.srwlib`: kept as an alias of `srwpy.srwlib` with a `DeprecationWarning`. Import `srwpy.srwlib` in new code.

### Removed

- `wpg/srw`: the bundled SRW of 2019 with binaries for Python 3.6.
- `setup.py`, `pytest.ini` and `requirements_win.txt`, replaced by `pyproject.toml`.
- `build/Makefile`: the build of the SRW fork of S. Yakubov (branch `openmp_memoryfix`) with a manually built FFTW2. The fork fixed a memory error of the OpenMP version. The fix has been in SRW since 2022 (commit bdc833c), and on the current SRW master the error is not reproduced (checked with AddressSanitizer).
- `Docker/base`, `Docker/develop` and `Docker/openmp`.
- `s2e/prop/utest/BeamlinesTest.py` and `Test.py`: they imported a module that no longer exists, the s2e tests are in `tests/test_helpers_fixes.py`. `ReferenceTest.py` (comparison with reference files on the XFEL server) is kept.

### Fixed

- Plots of `wpg_uti_wf` are faster for long pulses. `plot_wf` computes the intensity once instead of four times. `look_at_q_space` builds the angular map by FFT of the field on all cores (`angular_intensity`) instead of copying the wavefront and converting it with SRW, which converts the slices one by one in one thread. For a SACLA pulse of 201 × 201 points and 2662 slices on a 32-core server with SRW built with OpenMP (8 threads): `plot_wf` 0.70 → 0.27 s, `look_at_q_space` 6.97 → 1.71 s. For an odd number of points SRW took the angular step 1/((n − 1)·dx) instead of 1/(n·dx) and gave NaN at some points, so the angular FWHM printed for this pulse was 0.26 µrad in x and 4.4 µrad in y, now 4.4 µrad in both.
- Memory of `Wavefront` (issue #63). The parameter objects `wf.params`, `wf.params.Mesh` and `wf.data` kept a reference to the wavefront, so a wavefront was freed only by the cyclic garbage collector, late for large fields. In a loop of 40 propagations of a 10 MB pulse the memory of the process grew from 159 to 440 MB (9.5 MB per step), also with `del`. Now `wf.params` and `wf.data` are views created at each access, a wavefront is freed as soon as it is not used, and the memory stays at 100–110 MB.
- Beam parameters of `Wavefront`
  - `divergence` and the profiles along qx and qy took the intensity before the transform to angular space, swapped the axes and returned spatial frequency instead of angle (see Upgrade notes).
  - `get_axis`, `get_pixel_size` and `get_profile` converted the user's wavefront to another representation and did not restore the time or frequency domain. With an odd `nx` they filled the field with NaN.
  - `pulse_duration` and `get_beam_width` fitted in seconds starting from zero and could return a wrong or negative width.
  - `get_efield('vertical')` failed.
- `build/build_srw.py` failed with a short commit hash in `SRW_REF` (GitHub serves a single commit only by its full hash), now it fetches the history and resolves the hash.
- HDF5: `load_hdf5` kept the custom fields of the previously loaded file. Attributes holding lists of strings or bytes are now written and read correctly.
- `Beamline`
  - Operations by name (`index`, `edit_*`, `remove_element`, `crop_beamline`) ignored elements after `Empty`, and `edit_propagation_parameters` and `crop_beamline` failed with `AttributeError`.
  - `Empty` as the first or the only element crashed SRW.
  - Drifts inside a nested `SRWLOptC` did not shift `zCoord`.
- Optical elements and generators
  - `Mirror_plane_2d` laid the height map out so that a tilt along the mirror became a transverse tilt.
  - The normal of a mirror with `roll` was not a unit vector (its length was 1.00005 at 10 mrad).
  - `propagate_wavefront` changed its input wavefront.
  - `CRL` with the refractive index decrement and the attenuation length vs photon energy given as numpy arrays failed with `TypeError`, because SRW recognizes only lists (issue #92).
  - `build_gauss_wavefront` ignored the vertical beam quality factor `_my`, and `_mx` alone crashed SRW.
  - `build_gauss_wavefront` and `build_gauss_wavefront_xy` kept the radius of curvature that SRW sets, the distance to the waist z, instead of z·(1 + (z_R/z)²) with the Rayleigh length z_R (issue #133). SRW computes the field with the right curvature, but uses that radius for the semi-analytical treatment and resizing. A 9 keV beam with a 20 µm RMS waist (z_R = 36.5 m), generated 1 m from its waist on 256 × 256 points over 400 µm, had an RMS size of 148 µm after a 100 m drift with `Use_PP(semi_analytical_treatment=1)` and 1171 µm with zoom 4 and sampling 0.5. The analytic size is 58.86 µm. With the exact radius both give 58.8–58.9 µm.
  - `source.sase_pulse` mixed up the axes and the time and frequency domains and produced wrong vertical and circular polarizations (the remaining problems are in Known issues).
- GENESIS converters
  - `genesis.py` and `genesis_v3.py` lost the last slice.
  - `genesis_v2.py` swapped x and y and scaled once more the field that ocelot already gives in √W. After the fix the power agrees with the ocelot reader within a relative 1·10⁻⁵.
- Helpers
  - `useful_code.backpropagation.fit_gaussian` returned the x and y widths and centres swapped, because the first axis of `np.indices` is y. The function returned by `gaussian` now takes (y, x) instead of (x, y). `gaussian_parameters` stored in older files have x and y swapped.
  - `wfrutils.calculate_fwhm_x` took the centre from the x axis instead of y and returned FWHM = 0, and `plot_wfront` multiplied the power by the photon energy regardless of the field units.
  - `wpg_uti_wf.calculate_fwhm` measured an off-axis beam through the centre of the array instead of its maximum.
- `check_sampling` gave false alarms. Its Fresnel-zone test never succeeded for a converging beam (negative radius), used the SRW radius of curvature, which is wrong near the waist (see `build_gauss_wavefront` above), and reported a beam sampled above the Nyquist limit as undersampled.
- `AuxTransmAddSurfHeightProfileScaled` (used by `Mirror_plane` and the tutorials) is vectorized: the same result, 1.8 s → 0.03 s for a 1500 × 100 mesh. A profile with decreasing positions gave zeros.
- s2e: the history link to the parent file used a path relative to the current folder, the batch mode failed with `PicklingError` and hung on Linux with the OpenMP build of SRW (its worker processes were forked after SRW had started OpenMP threads, now they are started with `spawn`), and `diagnostics.py` failed on missing keys and on `plt.hold`.
- Incompatibilities with NumPy 2, h5py 3, Matplotlib 3 and Python 3.12: `copy=False`, `Dataset.value`, `tostring`, `pyplot.hold`, `imp`, Python 2 prints, `cPickle` in two sample notebooks (issue #128).

### Known issues

- SRW on GPU ignores the drift length for several slices without the semi-analytical treatment. A single slice, and several slices with `Use_PP(semi_analytical_treatment=1)`, agree with CPU within 10⁻⁶ of the maximum.
- SRW alone (`srwl.CalcElecFieldGaussian` called directly) sets the radius of curvature of a Gaussian beam to the distance to the waist, which breaks the semi-analytical treatment near the waist. The WPG generators correct it (see Fixed). For such wavefronts `check_sampling` reports the disagreement, and `wpg.sampling.set_radius_from_field(wf)` corrects the radius, also for a pulse in the time domain.
- Elements with the same name in a `Beamline`: operations by name find the last one.
- `source.sase_pulse`: the field normalization does not match its √(W/mm²) label, and the spectral grid is not conjugate to the time grid.
- The absolute normalization of `converters/genesis.py` and `genesis_v3.py` is not checked against a GENESIS file of known power.

[Unreleased]: https://github.com/samoylv/WPG/compare/develop...feature/wpg_update
