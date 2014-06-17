# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov, L. Samoylova'

import time
import numpy
import pylab

# Import standart libraries and addnig "../wavefront" directory to python
# search path
import os
import sys
sys.path.insert(0, os.path.join('..','..'))

from wpg import Wavefront, Beamline
# from srwlib import *
import wpg.srwlib
from wpg.srwlib import srwl\

J2EV = 6.24150934e18


def print_beamline(bl):
    mbl = Beamline(bl)
    print mbl


def create_numpy_array_from_rows(rows, slices=None):
    # slice size (Re, Im)
    N = len(rows[0]) / 2
    if slices is None:
        slices = range(len(rows) / N)
    slice_count = len(slices)
    # 3d array
    y = numpy.zeros(shape=(N, 2 * N, slice_count), dtype='float32')

    for si, s in enumerate(slices):
        for ii in range(N):
            y[ii, :, si] = rows[s * N + ii]
    return y


def plot_1d(profile, title_fig, title_x, title_y):
    # pylab.figure()
    pylab.plot(profile[0], profile[1])
    pylab.xlabel(title_x)
    pylab.ylabel(title_y)
    pylab.title(title_fig)
    pylab.grid(True)


def plot_2d(amap, xmin, xmax, ymin, ymax, title_fig, title_x, title_y):
    # pylab.figure()
    pylab.imshow(amap, extent=(ymin, ymax, xmin, xmax))
    pylab.colorbar()
    # pylab.axis('tight')
    pylab.xlabel(title_x)
    pylab.ylabel(title_y)
    pylab.title(title_fig)


def calculate_peak_pos(mwf):
    # irradiance
    irr = mwf.get_intensity(slice_number=0, polarization='vertical')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    nc = numpy.where(irr == irr_max)
    irr_x = irr[ny / 2, :]
    irr_y = irr[:, nx / 2]
    x0 = numpy.max(x_axis[numpy.where(irr_x == max(irr_x))])
    y0 = numpy.max(y_axis[numpy.where(irr_y == max(irr_y))])
    return [x0, y0]


def get_mesh(mwf):
    wf_mesh = mwf.params.Mesh
    nx = wf_mesh.nx
    ny = wf_mesh.ny
    [xmin, xmax, ymin, ymax] = [wf_mesh.xMin,
                                wf_mesh.xMax, wf_mesh.yMin, wf_mesh.yMax]
    return [nx, ny, xmin, xmax, ymin, ymax]


def plot_wfront(mwf, title_fig, isHlog, isVlog, i_x_min, i_y_min, orient, onePlot, bPlotPha=None,saveDir=None):
    """
        Plot 2D wavefront (a slice).
        
        :param mwf: 2D wavefront structure 
        :param title_fig: Figure title
        :param isHlog: if True, plot the horizontal cut in logarithmic scale
        :param isVlog: if True, plot the vertical cut in logarithmic scale
        :param i_x_min: Intensity threshold for horizontral cut, i.e. x-axis  limits are [min(where(i_x<i_x_min):max(where(i_x<i_x_min)]
        :param i_y_min: Intensity threshold for vertical cut,
        :param orient: 'x' for returning horizontal cut, 'y' for vertical cut
        :param onePlot: if True, put intensity map and plot of cuts on one plot, as  subplots
        :param bPlotPha: if True, plot the cuts of WF phase
        :return: 2-column array containing horizontal or vertical cut data in dependence of 'orient' parameter
        """
    if isHlog:
        print 'FWHMx[um]:', calculate_fwhm_x(mwf) * 1e6
    else:
        print 'FWHMx [mm]:', calculate_fwhm_x(mwf) * 1e3
    if isVlog:
        print 'FWHMy [um]:', calculate_fwhm_y(mwf) * 1e6
    else:
        print 'FWHMy [mm]:', calculate_fwhm_y(mwf) * 1e3
    [xc, yc] = calculate_peak_pos(mwf)
    print 'Coordinates of center, [mm]:', xc * 1e3, yc * 1e3
    ii = mwf.get_intensity(slice_number=0, polarization='vertical')
    # [LS14-06-02] 
    # for 2D Gaussian the intrincic SRW GsnBeam wave field units Nph/mm^2/0.1%BW 
    # to get fluence W/mm^2 
    # (Note: coherence time for Gaussian beam duration should be specified):  
    ii = ii*mwf.params.photonEnergy/J2EV#*1e3
    imax = numpy.max(ii)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    ph = mwf.get_phase(slice_number=0, polarization='vertical')
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    print 'stepX, stepY [um]:', dx * 1e6, dy * 1e6, '\n'
    xa = numpy.linspace(xmin, xmax, nx); 
    ya = numpy.linspace(ymin, ymax, ny); 

    print 'Total power (integrated over full range): %g [GW]' %(ii.sum(axis=0).sum(axis=0)*dx*dy*1e6*1e-9) 
    print 'Peak power calculated using FWHM:         %g [GW]' %(imax*1e-9*1e6*2*numpy.pi*(calculate_fwhm_x(mwf)/2.35)*(calculate_fwhm_y(mwf)/2.35))
    print 'Max irradiance: %g [GW/mm^2]'    %(imax*1e-9) 
    
    pylab.figure(figsize=(21,6))
    if onePlot:
        pylab.subplot(131)
    [x1, x2, y1, y2] = mwf.get_limits()
    pylab.imshow(ii, extent=[x1 * 1e3, x2 * 1e3, y1 * 1e3, y2 * 1e3])
    pylab.set_cmap('bone')
    #pylab.set_cmap('hot')
    pylab.axis('tight')
    #pylab.colorbar(orientation='horizontal')
    pylab.xlabel('x (mm)')
    pylab.ylabel('y (mm)')
    pylab.title(title_fig)

    irr_y = ii[:, numpy.max(numpy.where(xa == xc))]
    irr_x = ii[numpy.max(numpy.where(ya == yc)), :]
    pha_y = ph[:, numpy.max(numpy.where(xa == xc))]
    pha_x = ph[numpy.max(numpy.where(ya == yc)), :]

    if onePlot:
        pylab.subplot(132)
    else:
        pylab.figure()
    if isVlog and max(irr_y) > 0:
        #ya = ya*1e6
        pylab.semilogy(ya * 1e6, irr_y, '-vk')
        pylab.xlabel('(um)')
        pylab.xlim(min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e6, max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e6)
    else:
        #ya = ya*1e3
        pylab.plot(ya * 1e3, irr_y)
        pylab.xlabel('y (mm)')
        pylab.xlim(min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e3, max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3)
    pylab.ylabel('Irradiance (W/$mm^2$)')
    pylab.title('Vertical cut,  xc = ' + str(int(xc * 1e6)) + ' um')
    pylab.grid(True)
    if onePlot:
        pylab.subplot(133)
    else:
        pylab.figure()
    if isHlog and max(irr_x) > 0:
        #xa = xa*1e6
        pylab.semilogy(xa * 1e6, irr_x, '-vr')
        pylab.xlabel('x, (um)')
        pylab.xlim(min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e6, max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e6)
    else:
        #xa = xa*1e3
        pylab.plot(xa * 1e3, irr_x)
        pylab.xlabel('x (mm)')
        pylab.xlim(min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e3, max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3)
    pylab.ylabel('Irradiance (W/$mm^2$)')
    pylab.title('Horizontal cut, yc = ' + str(int(yc * 1e6)) + ' um')
    pylab.grid(True)
    
    if saveDir is not None: 
        epsname="%s/%s.eps" % (saveDir,title_fig.split("at ")[1].split(" m")[0])
        pylab.savefig(epsname)
        #pylab.close(epsfig)
    
    if bPlotPha:
        pylab.figure()
        pylab.plot(ya * 1e3, pha_y, '-ok')
        pylab.xlim(min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e3, max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3)
        pylab.ylim(-numpy.pi, numpy.pi)
        pylab.xlabel('y (mm)')
        pylab.title('phase, vertical cut, x=0')
        pylab.grid(True)

        pylab.figure()
        pylab.plot(xa * 1e3, pha_x, '-or')
        pylab.xlim(min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e3, max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3)
        pylab.ylim(-numpy.pi, numpy.pi)
        pylab.xlabel('x (mm)')
        pylab.title('phase, horizontal cut, y=0')
        pylab.grid(True)

    if orient == 'x':
        dd = numpy.zeros(shape=(nx, 2), dtype=float)
        dd[:, 0] = xa
        #for idx in range(nx): dd[idx, 1] = sum(ii[:, idx])
        dd[:, 1] = irr_x
    if orient == 'y':
        dd = numpy.zeros(shape=(ny, 2), dtype=float)
        dd[:, 0] = ya
        #for idx in range(ny): dd[idx, 1] = sum(ii[idx, :])
        dd[:, 1] = irr_y
    return dd


def calculate_fwhm(dd):
    irr_x = dd[:, 1]
    irr_max = numpy.max(irr_x)
    x_axis = dd[:, 0]
    fwhm = max(x_axis[numpy.where(irr_x >= irr_max / 2)]) - \
        min(x_axis[numpy.where(irr_x >= irr_max / 2)])
    return fwhm


def calculate_mediane(dd):
    irr_x = dd[:, 1]
    irr_max = numpy.max(irr_x)
    x_axis = dd[:, 0]
    mediane = (max(x_axis[numpy.where(irr_x >= irr_max / 2)])
               + min(x_axis[numpy.where(irr_x >= irr_max / 2)])) / 2
    return mediane


def calculate_fwhm_x(mwf):
    # irradiance
    irr = mwf.get_intensity(slice_number=0, polarization='vertical')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    [xc, yc] = calculate_peak_pos(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    irr_x = irr[numpy.max(numpy.where(y_axis == yc)), :]
    fwhm = 0.
    idx = numpy.where(irr_x >= irr_max / 2)
    if numpy.size(idx) > 0:
        fwhm = max(x_axis[numpy.where(irr_x >= irr_max / 2)]) - min(
            x_axis[numpy.where(irr_x >= irr_max / 2)])
    return fwhm


def calculate_fwhm_y(mwf):
    irr = mwf.get_intensity(slice_number=0, polarization='vertical')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    [xc, yc] = calculate_peak_pos(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    irr_y = irr[:, numpy.max(numpy.where(x_axis == xc))]
    fwhm = 0.
    idx = numpy.where(irr_y >= irr_max / 2)
    if numpy.size(idx) > 0:
        fwhm = max(y_axis[numpy.where(irr_y >= irr_max / 2)]) - min(
            y_axis[numpy.where(irr_y >= irr_max / 2)])
    return fwhm


def propagate_run(ifname, ofname, optBL, bSaved=False):
    """
        Propagate wavefront through a beamline and save the result (optionally).
        
        :param ifname: input hdf5 file name with wavefront to be propagated 
        :param ofname: output hdf5 file name
        :param optBL: beamline
        :param bSaved: if True, save propagated wavefront in h5 file
        :return: propagated wavefront
        """
    print_beamline(optBL)
    startTime = time.time()
    print '*****reading wavefront from h5 file...'
    w2 = Wavefront()
    w2.load_hdf5(ifname + '.h5')
    wfr = w2._srwl_wf
    print '*****propagating wavefront (with resizing)...'
    srwl.PropagElecField(wfr, optBL)
    mwf = Wavefront(wfr)
    print '[nx, ny, xmin, xmax, ymin, ymax]', get_mesh(mwf)
    if bSaved:
        print 'save hdf5:', ofname + '.h5'
        mwf.store_hdf5(ofname + '.h5')
    print 'done'
    print 'propagation lasted:', round((time.time() - startTime) / 6.) / 10., 'min'
    return wfr
