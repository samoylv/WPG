# WPG 

WPG, WavePropaGator, is an interactive simulation framework for X-ray wavefront propagation. 

WPG provides intuitive interface to the [SRW library](https://github.com/ochubar/SRW). The [application examples](http://wpg.readthedocs.org/en/latest/tutorials.html) are  mainly oriented on [European XFEL](http://www.xfel.eu) design parameters. To learn more details, see [online documentation pages](http://wpg.readthedocs.org/en/latest/index.html).

# Fast installation
## Python requiements

From relese 2019.12 we drop support of python2. Please use python 3.6 or 3.7.

## Get sources
Download sources https://github.com/samoylv/WPG/archive/develop.zip or clone git repository 
```bash
git clone https://github.com/samoylv/WPG.git
```

## Install python environment

Install conda from https://docs.conda.io/en/latest/miniconda.html or https://www.anaconda.com/distribution/ or from command line on Linux
```bash
wget -c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
```

Install dependencies
```bash
conda create -n wpg36 --file requirements.txt -c conda-forge -y
conda activate wpg36
```

## For Maxwell server users

For DESY Maxwell server you usually should install conda env to your GPSF partition, for example
```bash
conda create --prefix /gpfs/exfel/data/user/buzmakov/conda_env/wpg36 --file requirements.txt -c conda-forge -y
```

## Build Linux version
```bash
cd WPG
make
```

if you need SRW with OpenMP support (currently crashed on Maxwell in some cases)


```bash
cd WPG
OPENMP_MODE=omp make
```

## Get Windows binaries

WPG contain original SRW binaries for python 3.6 x64 with OpenMP support.

```bash
conda create -n wpg36 --file requirements_win.txt -c conda-forge -y
conda activate wpg36
```

Other windows binaries of SRW may be copied from original SRW repository (https://github.com/ochubar/SRW/tree/master/env/work/srw_python/lib) and placed in WPG/wpg/srw folder

## Run smaples

Run
````bash
cd WPG
jupyter notebook
````
Try run samples from WPG/samples/Tutorials





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


