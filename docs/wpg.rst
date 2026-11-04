About wpg
=========

The simulations based on wave optics have become indispensable for beamline design for highly coherent novel X-ray sources such as X-ray Free Electron Lasers (XFEL). 

We present a new interactive framework package for coherent and partially coherent X-ray wavefront propagation simulations – “WaveProperGator” (WPG). 

The package has been developed at European XFEL to facilitate for the end users (beamline scientists and XFEL users) the design, 
and optimization of X-ray optics to meet their experimental requirements. Our package uses the SRW C/C++ library and its Python binding for wavefront propagation simulations. The tool allows for changing source and optics parameters and visualizing the results interactively. 

The framework is cross-platform: it runs reliably on Linux, MS Windows 7, and Mac OS X. Using IPython as a web-frontend enables users to run the code on a remote server as well as on their local personal computer. One can use popular Python libraries (such as scipy, numpy, matplotlib) for pre- and post-processing as well as for visualization of the simulation results. 

The wavefronts are saved in hdf5 format for the eventual further processing and start-to-end simulations of experiments. The HDF5 format allows for keeping the calculation history within a single file, thus facilitating communication between various scientific groups, as well as cross-checking with other simulation results. The WPG source code together with guidelines for installation and application examples are available on the web. 

Several application examples, specific for XFEL, will be presented.

Source code aviable at GitHub https://github.com/samoylv/WPG

References
==========

If you use the WPG for your research, we would appreciate it if you would refer to the following papers:

1. Samoylova, L., Buzmakov, A., Chubar, O. & Sinn, H. WavePropaGator: Interactive framework for X-ray FEL optics design and simulations. // Journal of Applied Crystallography 08/2016; 49(4) pp. 1347-1355. DOI:10.1107/S160057671600995X 
http://journals.iucr.org/j/issues/2016/04/00/zd5006/index.html

See also 
========

1. WPG paper http://journals.iucr.org/j/issues/2016/04/00/zd5006/zd5006.pdf
2. XPD: XFEL Photon pulses Database
  https://in.xfel.eu/xpd/ and XPD API https://in.xfel.eu/xpd/doc/Rest_API.html
3. Simex project https://www.eucall.eu/organisation/work_packages/wp4___simex/

Getting started
===============

.. image:: _static/FELsource_out.png
    :scale: 100 %

Installation
------------

WPG is a pure Python package. The SRW library is not bundled any more,
it is installed as the ``srwpy`` package (official SRW builds from PyPI).
Python 3.10 or newer is required.

Standard installation (Linux, Mac OS X, MS Windows)
+++++++++++++++++++++++++++++++++++++++++++++++++++

Get the sources and install WPG together with its dependencies (``srwpy``, numpy, scipy, h5py, matplotlib):

	.. code:: sh

		git clone https://github.com/samoylv/WPG.git
		cd WPG
		python -m pip install -e .

The same with conda:

	.. code:: sh

		conda env create -n wpg -f environment.yml
		conda activate wpg
		python -m pip install --no-deps -e .

SRW built from sources with OpenMP or CUDA (Linux)
++++++++++++++++++++++++++++++++++++++++++++++++++

The official ``srwpy`` wheels are built without OpenMP and CUDA.
To use the parallel SRW versions, build ``srwpy`` from the
`SRW sources <https://github.com/ochubar/SRW>`_ into the active Python environment:

	.. code:: sh

		make srw SRW_MODE=omp      # OpenMP (FFTW2)
		make srw SRW_MODE=cuda     # CUDA (FFTW3 and cuFFT)
		make srw SRW_MODE=cpu      # serial build from the latest sources
		make install               # SRW (official wheel by default) and WPG

``SRW_REF`` selects the SRW branch, tag or commit (``master`` by default),
``SRW_SRC`` points to an existing SRW checkout instead of cloning.
The build script is ``build/build_srw.py``, run it with ``--help`` for all options.
OpenMP and CUDA cannot be combined in one build: SRW uses FFTW2 for the OpenMP version
and FFTW3 with cuFFT for the CUDA version.
Building requires gcc, git and CMake (installed automatically into the environment by the script),
the CUDA build also requires the CUDA toolkit (``CUDA_PATH``, ``/usr/local/cuda`` by default).

Run the samples
+++++++++++++++

   .. code:: sh

       cd samples
       jupyter notebook


If you have encounter errors while running the notebook containing messages such as

	.. code::

		ValueError: unknown locale: UTF-8

please try setting the following environment variables:

	.. code:: sh

	    export LC_ALL=en_US.UTF-8
	    export LANG=en_US.UTF-8
	    
Add these lines to your ~/.profile to make them permanent.	    

On xfel server
++++++++++++++

You can run the notebook server on a remote machine (e.g. a powerful cluster with heaps of memory) and connect to the server from your local machine via ssh tunnel. Do achieve this, do the following steps.

1.) Log in to the remote machine and issue the command

.. code:: sh

	jupyter-notebook --no-browser --port YOUR_UNIQUE_PORT_NUMBER --notebook-dir WPG_WORKING_DIRECTORY  

Here, YOUR\_UNIQUE\_PORT\_NUMBER is the port number that the server will use to communicate to the client later. Pick a number larger than 1024.

To avoid different users use the same port number, please visit `Simulation <https://docs.xfel.eu/share/page/site/simulation/dashboard>`_ web-site to register your port number in the `table  <https://docs.xfel.eu/share/page/site/simulation/wiki-page?title=List_of_port_numbers_in_use>`_.

	
2.) Setup ssh tunnel to the server. 
Leave the terminal window open from which you connected to the server and open a new terminal window. Issue the command

   .. code:: sh

       #If you have Linux or Mac OS X local machine, use the following command

       ssh -l <your_user_name> -f -N <server_name>  -LYOUR_UNIQUE_PORT_NUMBER:localhost:YOUR_UNIQUE_PORT_NUMBER

On Windows you can use putty and setup ssh tunnel as described here: `putty.pdf <https://github.com/samoylv/WPG/wiki/putty.pdf>`__

Open your local browser with the following web address: http://localhost:YOUR_UNIQUE_PORT_NUMBER


Useful links
------------

For visualization and browsing HDF5 files.

Please download HDFView tool from
[[http://www.hdfgroup.org/hdf-java-html/hdfview/]]

:mod:`wpg` module
=================

.. automodule:: wpg
