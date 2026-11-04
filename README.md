# WPG 

WPG, WavePropaGator, is an interactive simulation framework for X-ray wavefront propagation. 

WPG provides intuitive interface to the [SRW library](https://github.com/ochubar/SRW). The [application examples](http://wpg.readthedocs.org/en/latest/tutorials.html) are  mainly oriented on [European XFEL](http://www.xfel.eu) design parameters. To learn more details, see [online documentation pages](http://wpg.readthedocs.org/en/latest/index.html).

# Installation

WPG needs Python 3.10 or newer. SRW is not bundled with WPG any more: it is installed as the
[`srwpy`](https://pypi.org/project/srwpy/) package (official SRW builds for Linux, Mac OS X and Windows).

```bash
git clone https://github.com/samoylv/WPG.git
cd WPG
python -m pip install -e ".[test]" --config-settings editable_mode=compat
python -m pytest
```

`editable_mode=compat` lets IDEs and static analyzers (Pylance, pyright, mypy) find the editable WPG package.
Wavefront fields are declared in `wpg/glossary.py`, so `wf.params.Mesh.nx`, `wf.data.arrEhor` etc.
are known to autocompletion together with their types and units.

With conda:

```bash
conda env create -n wpg -f environment.yml
conda activate wpg
python -m pip install --no-deps -e .
```

Code written for older WPG versions with `from wpg.srwlib import ...` still works, but prints a
deprecation warning: import from `srwpy.srwlib` instead.

## SRW with OpenMP or CUDA (Linux)

The official `srwpy` wheels are serial. To use the parallel SRW versions, build `srwpy` from the
[SRW sources](https://github.com/ochubar/SRW) into the active Python environment:

```bash
make install SRW_MODE=omp     # SRW with OpenMP (FFTW2) + WPG
make install SRW_MODE=cuda    # SRW with CUDA (FFTW3 and cuFFT) + WPG
make srw SRW_MODE=cpu         # only SRW, serial build from the latest sources
```

* `SRW_REF` - SRW branch, tag or commit (`master` by default), `SRW_SRC` - existing SRW checkout instead of cloning;
* `SRW_CUDA_ARCHS=all` compiles CUDA code for all supported GPUs (80 86 89 90 120) instead of the GPU of this machine;
* `PYTHON` - Python of the target environment (`python3` by default).

The build script is [`build/build_srw.py`](build/build_srw.py) (`python build/build_srw.py --help`).
It needs gcc and git, CMake is installed into the environment by the script,
the CUDA build also needs the CUDA toolkit (`CUDA_PATH`, `/usr/local/cuda` by default).
The script checks that the installed `srwlpy` is linked with OpenMP or CUDA as requested and saves the build
description to `build/srw_build_info.json`.
OpenMP and CUDA cannot be combined in one build: SRW uses FFTW2 for the OpenMP version and FFTW3 with
cuFFT for the CUDA version.

With the CUDA build, propagation runs on GPU only when requested: `Beamline.propagate(wfr, gpu=1)`
(`gpu` is the number of the GPU device). This is experimental: SRW 4.2.1 on GPU agrees with CPU
for a single slice and, for several slices, with `Use_PP(semi_analytical_treatment=1)`, but for several
slices without semi-analytical treatment it ignores the drift length.

## Docker

See [Docker/readme.md](Docker/readme.md).

## Run samples

```bash
cd samples
jupyter notebook
```

Try the samples from `samples/Tutorials`.

# References

If you use the WPG for your research, we would appreciate it if you would refer to the following paper:

* Samoylova, L., Buzmakov, A., Chubar, O. & Sinn, H. WavePropaGator: Interactive framework for X-ray FEL optics design and simulations. // Journal of Applied Crystallography 08/2016; 49(4) pp. 1347-1355. DOI:10.1107/S160057671600995X http://journals.iucr.org/j/issues/2016/04/00/zd5006/index.html


# Primary works that use WPG 

* Yoon, Chun Hong, et al. "A comprehensive simulation framework for imaging single particles and biomolecules at the European X-ray Free-Electron Laser." Scientific reports 6 (2016): 24791.
* Fortmann-Grote, Carsten, et al. "Start-to-end simulation of single-particle imaging using ultra-short pulses at the European X-ray Free-Electron Laser." IUCrJ 4.5 (2017): 560-568.
* Manetti, Maurizio, et al. "XFEL Photon Pulses Database (XPD) for modeling XFEL experiments." (2016).
* Fortmann-Grote, C., et al. "Simex: Simulation of experiments at advanced light sources." arXiv preprint arXiv:1610.05980 (2016).
* Faenov, A. Ya, et al. "Advanced high resolution x-ray diagnostic for HEDP experiments." Scientific reports 8.1 (2018): 16407.
* Sinn, H., et al. "The SASE1 X-ray beam transport system." Journal of synchrotron radiation 26.3 (2019).
* Ruiz-Lopez, Mabel, et al. "Wavefront-propagation simulations supporting the design of a time-delay compensating monochromator beamline at FLASH2." Journal of synchrotron radiation 26.3 (2019).
* Roling, S., et al. "Time-dependent wave front propagation simulation of a hard x-ray split-and-delay unit: Towards a measurement of the temporal coherence properties of x-ray free electron lasers." Physical Review Special Topics-Accelerators and Beams 17.11 (2014): 110705.
* V. Kärcher, S. Roling, L. Samoylova, A. Buzmakov, U. Zastrau, K. Appel, M.V. Yurkov, E. Schneidmiller, F. Siewert, and H. Zacharias "Simulating wavefront propagation for a beam line with split-and-delay unit and compund refractive lenses at the European X-ray Free Electron Laser" // In press
* Guest, Trey. W., R. Bean, R. Kammering, G. van Riessen, A. P. Mancuso, and B. Abbey. “A Phenomenological Model of the X-Ray Pulse Statistics of a High-Repetition-Rate X-Ray Free-Electron Laser.” IUCrJ 10, no. 6 (November 1, 2023). https://doi.org/10.1107/S2052252523008242.


