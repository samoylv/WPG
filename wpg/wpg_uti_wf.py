# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__author__ = 'A. Buzmakov, L. Samoylova'

import time
import numpy
import copy
import pylab

# Import standart libraries and addnig "../wavefront" directory to python
# search path
#import os
#import sys
#sys.path.insert(0, os.path.join('..','..'))

from wpg.wavefront import Wavefront
from wpg.beamline import Beamline
# from srwlib import *
import wpg.srwlib
from wpg.srwlib import srwl


def print_mesh(wfr):
    """
    Print out wfr wavefront mesh.
    """    
 
    wf_mesh = wfr.params.Mesh;
    w_space=wfr.params.wSpace
    print(w_space)
    if (w_space=='R-space'):
        print( 'nx {:5d}  range_x [{:.1e}, {:.1e}] mm'.format(wf_mesh.nx,wf_mesh.xMin*1e3,wf_mesh.xMax*1e3))
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mm'.format(wf_mesh.ny,wf_mesh.yMin*1e3,wf_mesh.yMax*1e3))   
    if (w_space=='Q-space'):
        print('nx {:5d}  range_x [{:.1e}, {:.1e}] mrad'.format(wf_mesh.nx,wf_mesh.qxMin*1e3,wf_mesh.qxMax*1e3))
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mrad'.format(wf_mesh.ny,wf_mesh.qyMin*1e3,wf_mesh.qyMax*1e3))   
    return

def averaged_intensity(wf,bPlot=True):
    import pylab as plt
    """
    calculate the slice-to-slice integral intensity averaged over a meaningful range, mainly needed for processing spiky FEL source 

    :params: wf: wavefront structure 
    :params: bPlot: if True plot temporary pulse structure in the meaningful range 
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold
    """
    J2eV = 6.24150934e18
    #total0=wf.get_intensity().sum();
    dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)*1.e15/(wf.params.Mesh.nSlices-1)
    dx = (wf.params.Mesh.xMax - wf.params.Mesh.xMin)/(wf.params.Mesh.nx - 1)
    dy = (wf.params.Mesh.yMax - wf.params.Mesh.yMin)/(wf.params.Mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0); # I(t/slice_num)
    int0 = int0*(J2eV/wf.params.photonEnergy*dt*1.e-15*dx*dy*1.e4)
    #print(int0.shape)#,int0
    #print(J2eV/wf.params.photonEnergy*dt*1e-15*1e-4)
    int0max = max(int0)
    threshold = int0max * 0.01
    aw = numpy.argwhere(int0 > threshold)
    #print( aw.shape)
    int0_mean = int0[min(aw):max(aw)] # meaningful range of pulse
    #total0 = total0*J2eV/wf.params.photonEnergy*dx*dy*dt*1e-15*1e4/(dx*dy*1e12)# units: [ph/um^2], intrinsic wf: cm^2
    if bPlot: 
        Nph2W=wf.params.photonEnergy/(J2eV*dt*1e-15) # transfer Nphotons per px to Watts
        plt.figure();
        plt.plot(int0*Nph2W);
        plt.plot(numpy.arange(min(aw), max(aw)),int0_mean*Nph2W,'ro');plt.show()
    averaged = int0_mean.sum()/len(int0_mean)
    print('number of meaningful slices:',len(int0_mean))
    return averaged

def plot_t_wf(wf):
    import pylab as plt
    """
    plot wavefront in time domain (obligatory?) and R-space

    :params: wf: wavefront structure 
    """
    wf_intensity = wf.get_intensity().sum(axis=-1)
    average = averaged_intensity(wf,bPlot = True)
    nslices = wf.params.Mesh.nSlices
    dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)/(nslices-1)
    t0 = dt*nslices/2 + wf.params.Mesh.sliceMin
    plt.figure(figsize=(5,5),dpi=200)
    xmin,xmax,ymax,ymin = wf.get_limits()
    plt.imshow(wf_intensity, extent=[xmin*1e3,xmax*1e3,ymax*1e3,ymin*1e3])
    plt.xlabel('$mm$',fontsize=20); plt.ylabel('$mm$',fontsize=20);
    plt.axis('tight')
    plt.title('relative intensity={:03.3g}, t0={:03.2f} fs'.format(wf_intensity.sum()/average, t0*1.e15))
    plt.show()

def plot_t_wf_a(wf):
    import pylab
    """
    plot wavefront in Q-space

    :params: wf: wavefront structure 
    """
    wf_intensity = wf.get_intensity().sum(axis=-1)
    nslices = wf.params.Mesh.nSlices
    pylab.figure(figsize=(5,5),dpi=200)
    xmin,xmax,ymax,ymin = wf.get_limits()
    pylab.imshow(wf_intensity, extent=[xmin*1e6,xmax*1e6,ymax*1e6,ymin*1e6])
    pylab.xlabel('$\mu rad$',fontsize=20); pylab.ylabel('$\mu rad$',fontsize=20);
    pylab.axis('tight')
    #pylab.title('intensity={:03.3g} photons t={:03.2f} fs'.format(wf_intensity.sum(), dt*nslices/2))
    pylab.show()
    return

def look_at_q_space(wf, output_file = None):
    """
    change wavefront representation from R- to Q-space and store it in output file.
    
    :param wf: Wavefront object in R-space representation
    :param output_file: if parameter present - store propagaed wavefront to file
    :return: propagated wavefront object:
    """
    if True: #isinstance(wf, Wavefront):
        #print wf
        wfr = Wavefront(srwl_wavefront=wf._srwl_wf)
    else:
        print( 'no wavefront specified...')
        return

    if not wf.params.wSpace=='R-space':
        print( 'space should be in R-space, but not '+ wf.params.wSpace)
        return
    srwl_wf = wfr._srwl_wf
    srwl_wf_a = copy.deepcopy(srwl_wf)
    srwl.SetRepresElecField(srwl_wf_a, 'a')
    wf_a = Wavefront(srwl_wf_a)
    if not output_file is None:
        print('store wavefront to HDF5 file: '+ output_file+'...')
        wf_a.store_hdf5(output_file); print('done')

    print(  calculate_fwhm(wf_a))
    plot_t_wf_a(wf_a)
    return


def propagate_wavefront(wavefront, beamline, output_file = None):
    """
    Propagate wavefront and store it in output file.
    
    :param wavefront: Wavefront object or path to HDF5 file
    :param beamline: SRWLOptC container of beamline
    :param output_file: if parameter present - store propagaed wavefront to file
    :return: propagated wavefront object:
    """
    
    if not isinstance(beamline, Beamline):
        bl = Beamline(beamline)
    else:
        bl = beamline
    print(bl)

    if isinstance(wavefront, Wavefront):
        wfr = Wavefront(srwl_wavefront=wavefront._srwl_wf)
    else:
        print('*****reading wavefront from h5 file...')
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)

    print_mesh(wfr)
    print('*****propagating wavefront (with resizing)...')
    bl.propagate(wfr)

    if output_file is not None:
        print( 'save hdf5:', output_file)
        wfr.store_hdf5(output_file)
    print( 'done')
    return wfr


def calculate_fwhm(wfr):
    """
    Calculate FWHM of the beam calculating number of point bigger then max/2 throuhgt center of the image

    :param wfr:  wavefront
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m]
    """
#    intens = wfr.get_intensity(polarization='total')
    intens = wfr.get_intensity(polarization='total').sum(axis=-1);


    mesh = wfr.params.Mesh
    if (wfr.params.wSpace=='R-space'):
        dx = (mesh.xMax-mesh.xMin)/mesh.nx
        dy = (mesh.yMax-mesh.yMin)/mesh.ny
    elif (wfr.params.wSpace=='Q-space'):
        dx = (mesh.qxMax-mesh.qxMin)/mesh.nx
        dy = (mesh.qyMax-mesh.qyMin)/mesh.ny
    else:
        return

    x_center = intens[intens.shape[0]//2,:]
    fwhm_x = len(x_center[x_center>x_center.max()/2])*dx

    y_center = intens[:,intens.shape[1]//2]
    fwhm_y = len(y_center[y_center>y_center.max()/2])*dy
    if (wfr.params.wSpace=='Q-space'):
        print( wfr.params.wSpace)
        wl = 12.39*1e-10/(wfr.params.photonEnergy*1e-3)      #WaveLength
        fwhm_x = fwhm_x*wl
        fwhm_y = fwhm_y*wl

    return {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y}

def get_intensity_on_axis(wfr):
    """
    Calculate intensity (e.g. spectrum in frequency domain) along (x=y=0)
    :param wfr:  wavefront
    :return: [z,s0] in [a.u.] 
    """

    wf_intensity = wfr.get_intensity(polarization='horizontal')
    mesh = wfr.params.Mesh;
    dim = numpy.shape(wf_intensity) #array dimensions # <-to avoid wrong dimension assignment
    sz = numpy.zeros(shape=(mesh.nSlices, 2), dtype='float64')
    sz[:,0] = numpy.linspace(mesh.sliceMin, mesh.sliceMax, mesh.nSlices);
    sz[:,1] = wf_intensity[dim[0]/2,dim[1]/2, :] / wf_intensity.max() # <-to avoid wrong dimension assignment

    return sz
