About wpg
=========

The simulations based on wave optics have become indispensable for beamline design for highly coherent novel X-ray sources such as X-ray Free Electron Lasers (XFEL). 

We present a new interactive framework package for coherent and partially coherent X-ray wavefront propagation simulations – “WaveProperGator” (WPG). 

The package has been developed at European XFEL to facilitate for the end users (beamline scientists and XFEL users) the designing, optimizing and improving X-ray optics to meet their experimental requirements. Our package uses SRW C/C++ library and its Python binding for wavefront propagation simulations. The tool allows for changing source and optics parameters and visualizing the results interactively. 

The framework is cross-platform: it runs reliably on Linux, MS Windows 7, and Mac OS X. Using IPython as a web-frontend make the users possible run the code at a remote server as well as at their local personal computer. One can use popular Python libraries (such as scipy, numpy, matplotlib) for pre- and post-processing as well as for visualization of the simulation results. 

The wavefronts are saved in hdf5 format for the eventual further processing and start-to-end simulations of experiments. The HDF5 format allows for keeping the calculation history within a single file, thus facilitating communication between various scientific groups, as well as cross-checking with other simulation results. The WPG source code together with guidelines for installation and application examples are available on the web. 

Several application examples, specific for XFEL, will be presented.

Getting started
===============

.. image:: _static/FELsource_out.png
    :scale: 100 %

Installation
------------

On Ubuntu Desktop
+++++++++++++++++

Install dependencies

	.. code:: sh

	   sudo add-apt-repository -y ppa:jtaylor/ipython
	   sudo apt-get update
	   sudo apt-get install -y build-essential python-dev unzip python-numpy python-matplotlib 
	   sudo apt-get install -y python-pip python-scipy python-h5py ipython-notebook

Select the directory, in which WPG will be located

   .. code:: sh

       cd <your_working_directory>

Download and build library 

 .. code:: sh

		wget http://github.com/samoylv/WPG/archive/develop.zip

Extract package:

	.. code:: sh

		unzip develop.zip

Change the directory:

	.. code:: sh

		cd WPG-develop

Build library. This will download and build FFTW2 and SRW:

	.. code:: sh

		make all

Run web interface.

   .. code:: sh

       cd samples
       ipython notebook

If web page not pop up automatically, open your browser in http://localhost:8888

Mac OS X
++++++++

Install dependencies. You should have XCode and `MacPorts <http://www.macports.org/install.php>`__ installed.
   
   .. code:: sh
   
	   sudo port install py27-numpy py27-h5py py27-matplotlib

For ipython notebook:

	.. code:: sh

		sudo port install py27-tornado py27-zmq py27-nose py27-jinja2
		py27-sphinx py27-pygments py27-readline py27-ipython py27-scipy

For wget automatic downloading:

	.. code:: sh

		sudo port install wget

Select the directory, in which WPG will be located

	.. code:: sh

		sh cd

Download and build library

	.. code:: sh

		wget --no-check-certificate http://github.com/samoylv/WPG/archive/master.zip

Extract package

	.. code:: sh

		unzip master.zip

Change the directory

	.. code:: sh

		cd WPG-develop

Build library. This will download and build FFTW2 and SRW

	.. code:: sh

		make all

Run web interface.

   .. code:: sh

       cd samples
       ipython notebook

-  If web page not pop up automatically, open your browser in http://localhost:8888

If ypu have some errors runing IPython notebook, which ends with

	.. code::

		ValueError: unknown locale: UTF-8

, please add to file `~/.profile` the following commands and try restart your session or reboot:

	.. code:: sh

	    export LC_ALL=en_US.UTF-8
	    export LANG=en_US.UTF-8

On xfel server
++++++++++++++

You can directly run web interface using YOUR\_UNIQUE\_PORT\_NUMBER, and specifiy the path to the WPG installation on the XFEL server.

YOUR\_UNIQUE\_PORT\_NUMBER should be greater than 1024. 

Please visit `Simulation <https://docs.xfel.eu/share/page/site/simulation/dashboard>`_ web-site to register your port number in the `table  <https://docs.xfel.eu/share/page/site/simulation/wiki-page?title=List_of_port_numbers_in_use>`_ and get further information.

	.. code:: sh

	   ipython notebook --no-browser --port=YOUR_UNIQUE_PORT_NUMBER --notebook-dir=<your_working_directory> &

Setup ssh tunnel to the server. **Please use another LOCAL terminal window!**

   .. code:: sh

       #If you have Linux or Mac OS X local machine, use the following command

       ssh -l <your_account> -f -N <server_name>  -LYOUR_UNIQUE_PORT_NUMBER:localhost:YOUR_UNIQUE_PORT_NUMBER

On Windows you can use putty and setup ssh tunnel in the following way: `putty.pdf <https://github.com/samoylv/WPG/wiki/putty.pdf>`__

Open your local browser with the following web address: http://localhost:YOUR_UNIQUE_PORT_NUMBER


On MS Windows
+++++++++++++

You should have installed python2.7 with modules numpy, matplotlib, h5py and ipython. If these modules have not been installed yet, you can download a free python bundle with preinstalled packages `here <http://continuum.io/downloads>`_ 

Download `WPG package <https://github.com/samoylv/WPG/archive/master.zip>`_ and unpack it.

Download `SRW library <https://github.com/ochubar/SRW/archive/master.zip>`_ and unpack it in any folder.

Copy the following files from the SRW folder to WPG folder:	

- SRW-master/env/work/srw_python/lib/srwlpy2_x64.pyd to WGP/wpg

Rename srwlpy2_x64.pyd to srwlpy.pyd

Run ```ipython notebook ``` in ```WGP/samples```

If you have SRW already installed
+++++++++++++++++++++++++++++++++

Just copy ``srwlpy.so`` in 'wpg'
folder

Useful links
------------

For visualization and browsing HDF5 files.

Please download HDFVeiw tool from
[[http://www.hdfgroup.org/hdf-java-html/hdfview/]]

:mod:`wpg` module
=================

.. automodule:: wpg
