About wpg
=========

The simulations based on wave optics have become indispensable for beamline design for highly coherent novel X-ray sources such as X-ray Free Electron Lasers (XFEL). 

We present a new interactive framework package for coherent and partially coherent X-ray wavefront propagation simulations – “WaveProperGator” (WPG). 

The package has been developed at European XFEL to facilitate for the end users (beamline scientists and XFEL users) the design, 
and optimization of X-ray optics to meet their experimental requirements. Our package uses the SRW C/C++ library and its Python binding for wavefront propagation simulations. The tool allows for changing source and optics parameters and visualizing the results interactively. 

The framework is cross-platform: it runs reliably on Linux, MS Windows 7, and Mac OS X. Using IPython as a web-frontend enables users to run the code on a remote server as well as on their local personal computer. One can use popular Python libraries (such as scipy, numpy, matplotlib) for pre- and post-processing as well as for visualization of the simulation results. 

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

	   sudo apt-get update
	   sudo apt-get install -y build-essential python-dev unzip python-numpy python-matplotlib 
	   sudo apt-get install -y python-pip python-scipy python-h5py
	   sudo pip install jupyter 


Download sources:

 .. code:: sh

		wget http://github.com/samoylv/WPG/archive/develop.zip

Extract package (this will create a new directory "WPG-develop"):

	.. code:: sh

		unzip develop.zip

Change into the new directory:

	.. code:: sh

		cd WPG-develop

Build library. This will download and build FFTW2 and SRW:

	.. code:: sh

		make all

Start the jupyter notebook server. This should launch a new jupyter session in your web browser.

   .. code:: sh

       cd samples
       jupyter-notebook

If the notebook does not appear in a fresh web browser tab (window), direct your browoser to _`<http://localhost:8888 <http://localhost:8888>`_.

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

Download sources:

	.. code:: sh

		wget --no-check-certificate http://github.com/samoylv/WPG/archive/develop.zip

Extract package

	.. code:: sh

		unzip develop.zip

Change the directory

	.. code:: sh

		cd WPG-develop

Build library. This will download and build FFTW2 and SRW

	.. code:: sh

		make all

Run jupyter notebook server:

   .. code:: sh

       cd samples
       jupyter-notebook

- If the notebook does not appear in a fresh web browser tab (window), direct you browoser to _`<http://localhost:8888 <http://localhost:8888>`_.


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


On MS Windows
+++++++++++++

You should have installed python2.7 with modules numpy, matplotlib, h5py and ipython. If these modules have not been installed yet, you can download a free python bundle with preinstalled packages `here <http://continuum.io/downloads>`_ 

Download `WPG package <https://github.com/samoylv/WPG/archive/develop.zip>`_ and unpack it.

Download `Our small brunch of SRW library <https://github.com/buzmakov/SRW/archive/feature/srw_lite.zip>`_ and unpack it in any folder.

Copy the following files from the SRW folder to WPG folder:	

- SRW-feature-srw_lite/env/work/srw_python/lib/srwlpy2_x64.pyd to WGP/wpg

Rename srwlpy2_x64.pyd to srwlpy.pyd

Run ```ipython notebook ``` in ```WGP/samples```

If you have SRW already installed
+++++++++++++++++++++++++++++++++

Either add the directory SRW_Dev/env/work/srw_python to your $PYTHONPATH or copy the file srwlpy.so under that directory to WPG/wpg/.


Useful links
------------

For visualization and browsing HDF5 files.

Please download HDFView tool from
[[http://www.hdfgroup.org/hdf-java-html/hdfview/]]

:mod:`wpg` module
=================

.. automodule:: wpg
