# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__author__ = 'A. Buzmakov, L. Samoylova'

import time
import numpy
#import pylab

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

    wf_mesh = wfr.params.Mesh
    print('nx {:5d}  range_x [{:.1e}, {:.1e}]'.format(
        wf_mesh.nx, wf_mesh.xMin, wf_mesh.xMax))
    print('ny {:5d}  range_y [{:.1e}, {:.1e}]'.format(
        wf_mesh.ny, wf_mesh.yMin, wf_mesh.yMax))


def propagate_wavefront(wavefront, beamline, output_file=None):
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
        wfr = Wavefront(srwl_wavefront=wavefront._srw_wf)
    else:
        print('*****reading wavefront from h5 file...')
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)

    print_mesh(wfr)
    print('*****propagating wavefront (with resizing)...')
    bl.propagate(wfr)

    # if output_file is not None:
    #     print('save hdf5:', output_file)
    #     wfr.store_hdf5(output_file)
    # print('done')
    return wfr


def calculate_fwhm(wfr):
    """
    Calculate FWHM of the beam calculating number of point bigger then max/2 throuhgt center of the image

    :param wfr:  wavefront
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m]
    """
#    intens = wfr.get_intensity(polarization='total')
    intens = wfr.get_intensity(polarization='total').sum(axis=-1)

    mesh = wfr.params.Mesh
    dx = (mesh.xMax - mesh.xMin)/mesh.nx
    dy = (mesh.yMax - mesh.yMin)/mesh.ny

    x_center = intens[intens.shape[0]//2, :]
    fwhm_x = len(x_center[x_center > x_center.max()/2])*dx

    y_center = intens[:, intens.shape[1]//2]
    fwhm_y = len(y_center[y_center > y_center.max()/2])*dy
    return {'fwhm_x': fwhm_x, 'fwhm_y': fwhm_y}


def get_intensity_on_axis(wfr):
    """
    Calculate intensity (spectrum in frequency domain) along z-axis (x=y=0)
    :param wfr:  wavefront
    :return: [z,s0] in [a.u.] if frequency domain
    """

    wf_intensity = wfr.get_intensity(polarization='horizontal')
    mesh = wfr.params.Mesh
    zmin = mesh.sliceMin
    zmax = mesh.sliceMax
    sz = numpy.zeros((mesh.nSlices, 2), dtype='float64')
    sz[:, 0] = numpy.linspace(zmin, zmax, mesh.nSlices)
    sz[:, 1] = wf_intensity[mesh.nx/2, mesh.ny/2, :] / wf_intensity.max()

    return sz
