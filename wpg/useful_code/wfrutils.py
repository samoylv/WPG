# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

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
from wpg.wpg_uti_wf import propagate_wavefront
# from srwlib import *
import wpg.srwlib
from wpg.srwlib import srwl

J2EV = 6.24150934e18


def print_beamline(bl):
    if isinstance(bl, Beamline):
        print(bl)
    elif isinstance(bl, wpg.srwlib.SRWLOptC):
        mbl = Beamline(bl)
        print(mbl)
    else:
        raise ValueError(
            'Input type must be wpg.srwlib.SRWLOptC or wpg.Beamline, given: {}'.format(
                type(bl))
            )


def create_numpy_array_from_rows(rows, slices=None):
    # slice size (Re, Im)
    N = len(rows[0]) / 2
    if slices is None:
        slices = list(range(len(rows) / N))
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
    irr = mwf.get_intensity(slice_number=0, polarization='horizontal')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    nc = numpy.where(irr == irr_max)
    irr_x = irr[ny // 2, :]
    irr_y = irr[:, nx // 2]
    x0 = numpy.max(x_axis[numpy.where(irr_x == numpy.max(irr_x))])
    y0 = numpy.max(y_axis[numpy.where(irr_y == numpy.max(irr_y))])
    return [x0, y0]


def get_mesh(mwf):
    wf_mesh = mwf.params.Mesh
    nx = wf_mesh.nx
    ny = wf_mesh.ny
    [xmin, xmax, ymin, ymax] = [wf_mesh.xMin,
                                wf_mesh.xMax, wf_mesh.yMin, wf_mesh.yMax]
    return [nx, ny, xmin, xmax, ymin, ymax]


def show_slices_hsv(wf, slice_numbers=None, pretitle=''):
    """
        Show slices: intensity, phase, gaussian approximation parameters and cuts.
        @TBD:All gaussian parameters in pixels now. Should be fixed.
        @TBD: Add normalization to averaged slice intensity

        :params wf: wpg.Wavefront
        :params slice_numbers: slices to be shown, may by list, int, or None (for all slices)
        :params pretitle: string to be add in the beginning of the title line
        """

    from matplotlib.colors import hsv_to_rgb
    from wpg.useful_code.backpropagation import fit_gaussian, gaussian

    J2eV = 6.24150934e18
    wf_intensity = wf.get_intensity(polarization='horizontal')
    wf_phase = wf.get_phase(polarization='horizontal')
    if wf.params.wSpace == 'R-space':
        pulse_energy = wf.get_intensity().sum(axis=0).sum(axis=0).sum(axis=0)
        energyJ = calc_pulse_energy(wf)
        dx = (wf.params.Mesh.xMax - wf.params.Mesh.xMin) / \
            (wf.params.Mesh.nx - 1)
        dy = (wf.params.Mesh.yMax - wf.params.Mesh.yMin) / \
            (wf.params.Mesh.ny - 1)
    elif wf.params.wSpace == 'Q-space':
        dx = (wf.params.Mesh.qxMax - wf.params.Mesh.qxMin) / \
            (wf.params.Mesh.nx - 1)
        dy = (wf.params.Mesh.qyMax - wf.params.Mesh.qyMin) / \
            (wf.params.Mesh.ny - 1)
    else:
        raise TypeError('wSpace should be "R-space" or "Q-space"')

    dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / \
        (wf.params.Mesh.nSlices - 1)
    print('dt', dt)

    if slice_numbers is None:
        slice_numbers = range(wf_intensity.shape[-1])

    if isinstance(slice_numbers, int):
        slice_numbers = [slice_numbers, ]

    intense = wf_intensity.sum(0).sum(0)
    intense = numpy.squeeze(intense)
    intense = intense*dx*dy*1e6*1e-9  # [GW],  dx,dy [mm]
    print('Z coord: {0:.4f} m.'.format(wf.params.Mesh.zCoord))

    pylab.figure()
    if wf.params.wDomain == 'time':
        pylab.plot(numpy.linspace(wf.params.Mesh.sliceMin, wf.params.Mesh.sliceMax,
                                  wf.params.Mesh.nSlices)*1e15, intense)
        pylab.plot(slice_numbers, intense[slice_numbers], color='g', linestyle='None',
                   markersize=5, marker='o', markerfacecolor='w', markeredgecolor='g')
        pylab.title(pretitle+' Instanteneous power')
        pylab.xlim(wf.params.Mesh.sliceMin*1e15, wf.params.Mesh.sliceMax*1e15)
        pylab.xlabel('fs')
        pylab.ylabel('[GW]')
    else:  # if wDomain=='frequency'
        pylab.plot(numpy.linspace(-wf.params.Mesh.nSlices*dt/2, wf.params.Mesh.nSlices*dt/2,
                                  wf.params.Mesh.nSlices)/wf.params.photonEnergy*1e3, intense)
        pylab.plot((slice_numbers*dt-wf.params.Mesh.nSlices*dt/2)/wf.params.photonEnergy*1e3, intense[slice_numbers], color='g', linestyle='None',
                   markersize=5, marker='o', markerfacecolor='w', markeredgecolor='g')
        pylab.title(pretitle+' Spectrum')
        pylab.xlabel('$\Delta \omega / \omega _0 10^{3}$')
        pylab.ylabel('[a.u.]')
    pylab.show()

    total_intensity = wf_intensity.sum(axis=-1)
    data = total_intensity*dt

    fit_result = fit_gaussian(data)
    fit_result = fit_gaussian(data)
    (height, center_x, center_y, width_x, width_y) = fit_result['params']
    rsquared = fit_result['rsquared']
    fit = gaussian(height, center_x, center_y, width_x, width_y)
    fit_data = fit(*numpy.indices(data.shape))

    if wf.params.wDomain == 'time' and wf.params.wSpace == 'R-space':
        print('Total pulse intinsity {:.2f} [mJ]'.format(
            energyJ*1e3))
    print( '''Gaussian approximation parameters:
        center_x : {0:.2f}um.\t center_y : {1:.2f}um.
        width_x  : {2:.2f}um\t width_y : {3:.2f}um.
        rsquared : {4:0.4f}.'''.format((center_x-numpy.floor(wf.params.Mesh.nx/2))*dx*1e6,
                                       (center_y -
                                        numpy.floor(wf.params.Mesh.ny/2))*dy*1e6,
                                       width_x*dx*1e6, width_y*dy*1e6, rsquared))

    if wf.params.wSpace == 'R-space':
        x_axis = numpy.linspace(
            wf.params.Mesh.xMin, wf.params.Mesh.xMax, wf.params.Mesh.nx)
    elif wf.params.wSpace == 'Q-space':
        x_axis = numpy.linspace(
            wf.params.Mesh.qxMin, wf.params.Mesh.qxMax, wf.params.Mesh.nx)
    else:
        raise TypeError('wSpace should be "R-space" or "Q-space"')
    y_axis = x_axis

    pylab.figure(figsize=(15, 7))
    pylab.subplot(121)
    pylab.imshow(
        data*dx*dy*1e6*J2eV/wf.params.photonEnergy, extent=wf.get_limits())
    pylab.colorbar(orientation='horizontal')
    if wf.params.wSpace == 'R-space':
        pylab.title('Nphotons per ' + str(numpy.floor(dx*1e6)) +
                    'x'+str(numpy.floor(dx*1e6))+' $\mu m ^2$ pixel')

    pylab.subplot(122)
    pylab.plot(y_axis*1e6,     data[:, int(center_x)]*1e3, 'b', label='Y-cut')
    pylab.hold(True)
    pylab.plot(
        y_axis*1e6, fit_data[:, int(center_x)]*1e3, 'b:', label='Gaussian fit')
    pylab.hold(True)
    pylab.plot(x_axis*1e6,     data[int(center_y), :]*1e3,  'g', label='X-cut')
    pylab.hold(True)
    pylab.plot(
        x_axis*1e6, fit_data[int(center_y), :]*1e3,  'g--', label='Gaussian fit')
    pylab.xlabel('[$\mu$m]')
    pylab.ylabel('mJ/mm$^2$')
    pylab.grid(True)
    pylab.legend()

    pylab.show()

    for sn in slice_numbers:
        data = wf_intensity[:, :, sn]
        data = data*dt
        phase = wf_phase[:, :, sn]
        fit_result = fit_gaussian(data)
        fit_result = fit_gaussian(data)
        (height, center_x, center_y, width_x, width_y) = fit_result['params']
        rsquared = fit_result['rsquared']
        fit = gaussian(height, center_x, center_y, width_x, width_y)
        fit_data = fit(*numpy.indices(data.shape))
        #$center_x = int(wf.params.Mesh.nSlices/2); center_y = center_x

        print('Slice number: {}'.format(sn))
        print( '''Gaussian approximation parameters:
            center_x : {0:.2f}um.\t center_y : {1:.2f}um.
            width_x  : {2:.2f}um\t width_y : {3:.2f}um.
            rsquared : {4:0.4f}.'''.format((center_x-numpy.floor(wf.params.Mesh.nx/2))*dx*1e6,
                                           (center_y -
                                            numpy.floor(wf.params.Mesh.ny/2))*dy*1e6,
                                           width_x*dx*1e6, width_y*dy*1e6, rsquared))

        pylab.figure(figsize=(15, 7))

        pylab.subplot(121)
        # number of photons in a slice of thickness dt
        intensity = data*dx*dy*1e6*J2eV/wf.params.photonEnergy*1e-6
        phase = wf_phase[:, :, sn]

        H = intensity
        V = phase
        S = numpy.ones_like(V)
        # V=(V-V.min())/(V.max()-V.min())

        h = (H-H.min())/(H.max()-H.min())
        v = V/(2*numpy.pi)+0.5

        HSV = numpy.dstack((v, S, h))
        RGB = hsv_to_rgb(HSV)
        pylab.imshow(RGB)
        pylab.title('Nphotons x10$^6$ per ' + str(numpy.floor(dx*1e6)) +
                    'x'+str(numpy.floor(dx*1e6))+' $\mu m ^2$ pixel')
#         pylab.contour(fit_data, cmap=pylab.cm.copper)

#         pylab.subplot(142)
#         pylab.imshow(wf_phase[:, :, sn])

#         pylab.colorbar(orientation='horizontal')

        pylab.subplot(143)
        pylab.plot(
            y_axis*1e6,     data[:, int(center_x)]*1e3,  'b', label='Y-cut')
        pylab.hold(True)
        pylab.plot(
            y_axis*1e6, fit_data[:, int(center_x)]*1e3,  'b:', label='Gaussian fit')
        pylab.hold(True)
        pylab.plot(
            x_axis*1e6,     data[int(center_y), :]*1e3,  'g', label='X-cut')
        pylab.hold(True)
        pylab.plot(
            x_axis*1e6, fit_data[int(center_y), :]*1e3,  'g--', label='Gaussian fit')
        pylab.xlabel('[$\mu$m]')
        pylab.ylabel('mJ/mm$^2$')
        pylab.grid(True)
        pylab.legend()

        pylab.subplot(144)
        pylab.plot(
            y_axis*1e6, phase[:, int(center_x)], label='Y-cut', marker='d', markersize=4)
        pylab.plot(
            x_axis*1e6, phase[int(center_y), :], label='X-cut', marker='o', markersize=4)
        pylab.xlabel('[$\mu$m]')
        pylab.legend()

        pylab.show()

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
        print('FWHMx[um]:', calculate_fwhm_x(mwf) * 1e6)
    else:
        print('FWHMx [mm]:', calculate_fwhm_x(mwf) * 1e3)
    if isVlog:
        print('FWHMy [um]:', calculate_fwhm_y(mwf) * 1e6)
    else:
        print('FWHMy [mm]:', calculate_fwhm_y(mwf) * 1e3)
    [xc, yc] = calculate_peak_pos(mwf)
    print('Coordinates of center, [mm]:', xc * 1e3, yc * 1e3)
    ii = mwf.get_intensity(slice_number=0, polarization='horizontal')
    # [LS14-06-02] 
    # for 2D Gaussian the intrincic SRW GsnBeam wave field units Nph/mm^2/0.1%BW 
    # to get fluence W/mm^2 
    # (Note: coherence time for Gaussian beam duration should be specified):  
    ii = ii*mwf.params.photonEnergy/J2EV#*1e3
    imax = numpy.max(ii)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    ph = mwf.get_phase(slice_number=0, polarization='horizontal')
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    print('stepX, stepY [um]:', dx * 1e6, dy * 1e6, '\n')
    xa = numpy.linspace(xmin, xmax, nx); 
    ya = numpy.linspace(ymin, ymax, ny); 

    if mwf.params.wEFieldUnit != 'arbitrary':
        print('Total power (integrated over full range): %g [GW]' %(ii.sum(axis=0).sum(axis=0)*dx*dy*1e6*1e-9)) 
        print('Peak power calculated using FWHM:         %g [GW]' %(imax*1e-9*1e6*2*numpy.pi*(calculate_fwhm_x(mwf)/2.35)*(calculate_fwhm_y(mwf)/2.35)))
        print('Max irradiance: %g [GW/mm^2]'    %(imax*1e-9)) 
        label4irradiance = 'Irradiance (W/$mm^2$)'
    else:
        ii = ii / imax
        label4irradiance = 'Irradiance (a.u.)'
    
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
    if isVlog and numpy.max(irr_y) > 0:
        #ya = ya*1e6
        pylab.semilogy(ya * 1e6, irr_y, '-vk')
        pylab.xlabel('(um)')
        pylab.xlim(numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e6, numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e6)
    else:
        #ya = ya*1e3
        pylab.plot(ya * 1e3, irr_y)
        pylab.xlabel('y (mm)')
        pylab.xlim(numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e3, numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3)
    pylab.ylim(0,numpy.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Vertical cut,  xc = ' + str(int(xc * 1e6)) + ' um')
    pylab.grid(True)
    if onePlot:
        pylab.subplot(133)
    else:
        pylab.figure()
    if isHlog and numpy.max(irr_x) > 0:
        #xa = xa*1e6
        pylab.semilogy(xa * 1e6, irr_x, '-vr')
        pylab.xlabel('x, (um)')
        pylab.xlim(numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e6, numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e6)
    else:
        #xa = xa*1e3
        pylab.plot(xa * 1e3, irr_x)
        pylab.xlabel('x (mm)')
        pylab.xlim(numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e3, numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3)
    pylab.ylim(0,numpy.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Horizontal cut, yc = ' + str(int(yc * 1e6)) + ' um')
    pylab.grid(True)
    
    if saveDir is not None: 
        epsname="%s/%s.eps" % (saveDir,title_fig.split("at ")[1].split(" m")[0])
        pylab.savefig(epsname)
        #pylab.close(epsfig)
    
    if bPlotPha:
        pylab.figure()
        pylab.plot(ya * 1e3, pha_y, '-ok')
        pylab.xlim(numpy.min(ya[numpy.where(irr_y >= imax * i_y_min)])
                   * 1e3, numpy.max(ya[numpy.where(irr_y >= imax * i_y_min)]) * 1e3)
        pylab.ylim(-numpy.pi, numpy.pi)
        pylab.xlabel('y (mm)')
        pylab.title('phase, vertical cut, x=0')
        pylab.grid(True)

        pylab.figure()
        pylab.plot(xa * 1e3, pha_x, '-or')
        pylab.xlim(numpy.min(xa[numpy.where(irr_x >= imax * i_x_min)])
                   * 1e3, numpy.max(xa[numpy.where(irr_x >= imax * i_x_min)]) * 1e3)
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


def stat1(x, y):
    """
    Calculate statistic moments of y(x) data.

    :param x: variable
    :param y: distribution y(x)
    :return: 'expected value' and 'variance'
    """
    mu = numpy.average(x, weights=y)  # expected value
    var = numpy.sqrt(numpy.average((x - mu)**2, weights=y))  # variance
    return mu, var


def calculate_fwhm(dd):
    irr_x = dd[:, 1]
    irr_max = numpy.max(irr_x)
    x_axis = dd[:, 0]
    fwhm = numpy.max(x_axis[numpy.where(irr_x >= irr_max / 2)]) - \
        numpy.min(x_axis[numpy.where(irr_x >= irr_max / 2)])
    return fwhm


def calculate_mediane(dd):
    irr_x = dd[:, 1]
    irr_max = numpy.max(irr_x)
    x_axis = dd[:, 0]
    mediane = (numpy.max(x_axis[numpy.where(irr_x >= irr_max / 2)])
               + numpy.min(x_axis[numpy.where(irr_x >= irr_max / 2)])) / 2
    return mediane


def calculate_fwhm_x(mwf):
    # irradiance
    irr = mwf.get_intensity(slice_number=0, polarization='horizontal')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    [xc, yc] = calculate_peak_pos(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    irr_x = irr[numpy.max(numpy.where(y_axis == yc)), :]
    fwhm = 0.
    idx = numpy.where(irr_x >= irr_max / 2)
    if numpy.size(idx) > 0:
        fwhm = numpy.max(x_axis[numpy.where(irr_x >= irr_max / 2)]) - numpy.min(
            x_axis[numpy.where(irr_x >= irr_max / 2)])
    return fwhm


def calculate_fwhm_y(mwf):
    irr = mwf.get_intensity(slice_number=0, polarization='horizontal')
    irr_max = numpy.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    [xc, yc] = calculate_peak_pos(mwf)
    x_axis = numpy.linspace(xmin, xmax, nx)
    y_axis = numpy.linspace(ymin, ymax, ny)
    irr_y = irr[:, numpy.max(numpy.where(x_axis == xc))]
    fwhm = 0.
    idx = numpy.where(irr_y >= irr_max / 2)
    if numpy.size(idx) > 0:
        fwhm = numpy.max(y_axis[numpy.where(irr_y >= irr_max / 2)]) - numpy.min(
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
    print('*****reading wavefront from h5 file...')
    w2 = Wavefront()
    w2.load_hdf5(ifname + '.h5')
    wfr = w2._srwl_wf
    print('*****propagating wavefront (with resizing)...')
    srwl.PropagElecField(wfr, optBL)
    mwf = Wavefront(wfr)
    print('[nx, ny, xmin, xmax, ymin, ymax]', get_mesh(mwf))
    if bSaved:
        print('save hdf5:', ofname + '.h5')
        mwf.store_hdf5(ofname + '.h5')
    print('done')
    print('propagation lasted:', round((time.time() - startTime) / 6.) / 10., 'min')
    return wfr

