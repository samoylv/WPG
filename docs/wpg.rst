Introduction
============

About wpg
---------

The simulations based on wave optics have become indispensable for beamline design for highly coherent novel X-ray sources such as X-ray Free Electron Lasers (XFEL). We present a new interactive framework package for coherent and partially coherent X-ray wavefront propagation simulations – “WaveProperGator” (WPG). The package has been developed at European XFEL to facilitate for the end users (beamline scientists and XFEL users) the designing, optimizing and improving X-ray optics to meet their experimental requirements. Our package uses SRW C/C++ library and its Python binding for wavefront propagation simulations. The tool allows for changing source and optics parameters and visualizing the results interactively. The framework is cross-platform: it runs reliably on Linux, MS Windows 7, and Mac OS X. Using IPython as a web-frontend make the users possible run the code at a remote server as well as at their local personal computer. One can use popular Python libraries (such as scipy, numpy, matplotlib) for pre- and post-processing as well as for visualization of the simulation results. The wavefronts are saved in hdf5 format for the eventual further processing and start-to-end simulations of experiments. The HDF5 format allows for keeping the calculation history within a single file, thus facilitating communication between various scientific groups, as well as cross-checking with other simulation results. The WPG source code together with guidelines for installation and application examples are available on the web. Several application examples, specific for XFEL, will be presented.

Installation
------------

Useful links
------------

:mod:`wpg` module
----------------------------

.. automodule:: wpg