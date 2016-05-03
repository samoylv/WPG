# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import copy
import numpy
import pylab

from wpg.beamline import Beamline
from wpg.srwlib import srwl
from wpg.wavefront import Wavefront

__author__ = 'A. Buzmakov, L. Samoylova, C. Fortmann-Grote'


def print_mesh(wf):
    """
    Print out wfr wavefront mesh.
    """

    wf_mesh = wf.params.Mesh
    w_space = wf.params.wSpace
    print(w_space)
    if (w_space == 'R-space'):
        print('nx {:5d}  range_x [{:.1e}, {:.1e}] mm'.format(
            wf_mesh.nx, wf_mesh.xMin*1e3, wf_mesh.xMax*1e3)
        )
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mm'.format(
            wf_mesh.ny, wf_mesh.yMin*1e3, wf_mesh.yMax*1e3)
        )
    if (w_space == 'Q-space'):
        print('nx {:5d}  range_x [{:.1e}, {:.1e}] mrad'.format(
            wf_mesh.nx, wf_mesh.qxMin*1e3, wf_mesh.qxMax*1e3)
        )
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mrad'.format(
            wf_mesh.ny, wf_mesh.qyMin*1e3, wf_mesh.qyMax*1e3)
        )
    return


def calc_pulse_energy(wf):
    """
    calculate energy of  in time domain
    params: wf: wavefront structure
    return: pulse energy value in [J]
    """
    J2eV = 6.24150934e18
    if wf.params.wDomain != 'time':
        print('Pulse energy cannot be calculated for {:s} domain'.format(
            wf.params.wDomain))
        return None
    else:
        dx = (wf.params.Mesh.xMax - wf.params.Mesh.xMin) / \
            (wf.params.Mesh.nx - 1)
        dy = (wf.params.Mesh.yMax - wf.params.Mesh.yMin) / \
            (wf.params.Mesh.ny - 1)
        dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / \
            (wf.params.Mesh.nSlices - 1)
        pulse_energy = wf.get_intensity().sum(axis=0).sum(axis=0).sum(axis=0)
        pulse_energy_J = pulse_energy*dx*dx*1e6*dt
        print('Number of photons per pulse: {:e}'.format(
            pulse_energy_J*J2eV/wf.params.photonEnergy))
        return pulse_energy_J


def averaged_intensity(wf, bPlot=True):
    """
    calculate the slice-to-slice integral intensity averaged over a meaningful range, mainly needed for processing spiky FEL source

    :params: wf: wavefront structure
    :params: bPlot: if True plot temporary pulse structure in the meaningful range
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold
    """
    J2eV = 6.24150934e18
    # total0=wf.get_intensity().sum();
    dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin) * \
        1.e15/(wf.params.Mesh.nSlices-1)
    dx = (wf.params.Mesh.xMax - wf.params.Mesh.xMin)/(wf.params.Mesh.nx - 1)
    dy = (wf.params.Mesh.yMax - wf.params.Mesh.yMin)/(wf.params.Mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)  # I(t/slice_num)
    int0 = int0*(J2eV/wf.params.photonEnergy*dt*1.e-15*dx*dy*1.e4)
    # print(int0.shape)#,int0
    # print(J2eV/wf.params.photonEnergy*dt*1e-15*1e-4)
    int0max = max(int0)
    threshold = int0max * 0.01
    aw = numpy.argwhere(int0 > threshold)
    #print( aw.shape)
    int0_mean = int0[min(aw):max(aw)]  # meaningful range of pulse
    # total0 =
    # total0*J2eV/wf.params.photonEnergy*dx*dy*dt*1e-15*1e4/(dx*dy*1e12)#
    # units: [ph/um^2], intrinsic wf: cm^2
    if bPlot:
        # transfer Nphotons per px to Watts
        Nph2W = wf.params.photonEnergy/(J2eV*dt*1e-15)
        pylab.figure()
        pylab.plot(int0*Nph2W)
        pylab.plot(numpy.arange(min(aw), max(aw)), int0_mean*Nph2W, 'ro')
        pylab.show()
    averaged = int0_mean.sum()/len(int0_mean)
    print('number of meaningful slices:', len(int0_mean))
    return averaged


def plot_t_wf(wf, save='', range_x=None, range_y=None):
    """
    Plot wavefront in  R-space.

    :params: wf: wavefront structure

    :params: save: Whether to save the figure on disk
    :type:  string for filename. Empty string '' means don't save.
    :default: '', do not save the figure.

    :params: range_x: x-axis range.
    :type: float
    :default: None, take entire x range.

    :params: range_y: y-ayis range.
    :type: float
    :default: None, take entire y range.
    """
    import matplotlib.pyplot as plt
    # Get the wavefront and integrate over time.
    wf_intensity = wf.get_intensity().sum(axis=-1)

    # Get average and time slicing.
    average = averaged_intensity(wf, bPlot=True)
    nslices = wf.params.Mesh.nSlices
    dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)/(nslices-1)
    t0 = dt*nslices/2 + wf.params.Mesh.sliceMin

    # Setup a figure.
    figure = plt.figure(figsize=(10, 10), dpi=100)
    plt.axis('tight')
    # Set colormap. Inferno is not available in all matplotlib versions, fall
    # back to gnuplot.
    try:
        plt.set_cmap('viridis')
    except:
        plt.set_cmap('gnuplot')

    # Profile plot.
    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    # Get limits.
    xmin, xmax, ymax, ymin = wf.get_limits()

    # Plot profile as 2D colorcoded map.
    profile.imshow(
        wf_intensity, extent=[xmin*1.e3, xmax*1.e3, ymax*1.e3, ymin*1.e3])
    profile.set_aspect('equal', 'datalim')

    # Get x and y ranges.
    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    #x = numpy.linspace(xmin*1.e3,xmax*1.e3,wf_intensity.shape[0])
    #y = numpy.linspace(ymin*1.e3,ymax*1.e3,wf_intensity.shape[1])
    x = numpy.linspace(xmin*1.e3, xmax*1.e3, wf_intensity.shape[1])
    y = numpy.linspace(ymin*1.e3, ymax*1.e3, wf_intensity.shape[0])

    # Labels.
    profile.set_xlabel('$mm$', fontsize=12)
    profile.set_ylabel('$mm$', fontsize=12)

    # x-projection plots above main plot.
    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    print(x.shape, wf_intensity.sum(axis=0).shape)

    x_projection.plot(x, wf_intensity.sum(axis=0), label='x projection')

    # Set range according to input.
    if range_x is None:
        profile.set_xlim([xmin*1.e3, xmax*1.e3])
    else:
        profile.set_xlim([-range_x/2., range_x/2.])

    # Set title.
    x_projection.set_title('relative intensity={:03.3g}, t0={:03.2f} fs'.format(
        wf_intensity.sum()/average, t0*1.e15))

    # y-projection plot right of main plot.
    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label='y projection')

    # Hide minor tick labels, they disturb here.
    plt.minorticks_off()

    # Set range according to input.
    if range_y is None:
        profile.set_ylim([ymin*1.e3, ymax*1.e3])
    else:
        profile.set_ylim([-range_y/2., range_y/2.])

    # If requested, save to disk, otherwise show in interactive window.
    if save != '':
        # Add parameters.
        plt.savefig(save)
    else:
        plt.show()


def plot_t_wf_a(wf, save='', range_x=None, range_y=None):
    """
    Plot wavefront in Q-space.

    :params: wf: wavefront structure

    :params: save: Whether to save the figure on disk
    :type:  string for filename. Empty string '' means don't save.
    :default: '', do not save the figure.

    :params: range_x: x-axis range.
    :type: float
    :default: None, take entire x range.

    :params: range_y: y-ayis range.
    :type: float
    :default: None, take entire y range.

    """
    import matplotlib.pyplot as plt
    wf_intensity = wf.get_intensity().sum(axis=-1)
    plt.figure(figsize=(10, 10), dpi=100)
    plt.axis('tight')
    try:
        plt.set_cmap('viridis')
    except:
        plt.set_cmap('gnuplot')

    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    xmin, xmax, ymax, ymin = wf.get_limits()

    profile.imshow(
        wf_intensity, extent=[xmin*1.e6, xmax*1.e6, ymax*1.e6, ymin*1.e6])
    profile.set_aspect('equal', 'datalim')

    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    #x = numpy.linspace(xmin*1.e6,xmax*1.e6,wf_intensity.shape[0])
    #y = numpy.linspace(ymin*1.e6,ymax*1.e6,wf_intensity.shape[1])
    x = numpy.linspace(xmin*1.e6, xmax*1.e6, wf_intensity.shape[1])
    y = numpy.linspace(ymin*1.e6, ymax*1.e6, wf_intensity.shape[0])
    profile.set_xlabel(r'$\mu$rad', fontsize=12)
    profile.set_ylabel(r'$\mu$rad', fontsize=12)

    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    print(x.shape, wf_intensity.sum(axis=0).shape)
    x_projection.plot(x, wf_intensity.sum(axis=0), label='x projection')
    if range_x is None:
        profile.set_xlim([xmin*1.e6, xmax*1.e6])
    else:
        profile.set_xlim([-range_x/2., range_x/2.])

    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label='y projection')

    # Hide minor tick labels.
    plt.minorticks_off()

    if range_y is None:
        profile.set_ylim([ymin*1.e6, ymax*1.e6])
    else:
        profile.set_ylim([-range_y/2., range_y/2.])

    if save != '':
        plt.savefig(save)
    else:
        plt.show()


def look_at_q_space(wf, output_file=None, save='', range_x=None, range_y=None):
    """
    change wavefront representation from R- to Q-space and store it in output file.

    :params wf: Wavefront object in R-space representation
    :params output_file: if parameter present - store propagaed wavefront to file

    :params save: Whether to save the figure on disk
    :type:  string for filename. Empty string '' means don't save.
    :default: '', do not save the figure.

    :params range_x: x-axis range.
    :type: float
    :default: None, take entire x range.

    :params range_y: y-ayis range.
    :type: float
    :default: None, take entire y range.

    :return: propagated wavefront object:
    """
    wfr = Wavefront(srwl_wavefront=wf._srwl_wf)

    if not wf.params.wSpace == 'R-space':
        print('space should be in R-space, but not ' + wf.params.wSpace)
        return
    srwl_wf = wfr._srwl_wf
    srwl_wf_a = copy.deepcopy(srwl_wf)
    srwl.SetRepresElecField(srwl_wf_a, 'a')
    wf_a = Wavefront(srwl_wf_a)
    if output_file is not None:
        print('store wavefront to HDF5 file: ' + output_file+'...')
        wf_a.store_hdf5(output_file)
        print('done')

    print(calculate_fwhm(wf_a))
    plot_t_wf_a(wf_a, save=save, range_x=range_x, range_y=range_y)
    return


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
        wfr = Wavefront(srwl_wavefront=wavefront._srwl_wf)
    else:
        print('*****reading wavefront from h5 file...')
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)

    print_mesh(wfr)
    print('*****propagating wavefront (with resizing)...')
    bl.propagate(wfr)

    if output_file is not None:
        print('save hdf5:', output_file)
        wfr.store_hdf5(output_file)
    print('done')
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
    if (wfr.params.wSpace == 'R-space'):
        dx = (mesh.xMax-mesh.xMin)/mesh.nx
        dy = (mesh.yMax-mesh.yMin)/mesh.ny
    elif (wfr.params.wSpace == 'Q-space'):
        dx = (mesh.qxMax-mesh.qxMin)/mesh.nx
        dy = (mesh.qyMax-mesh.qyMin)/mesh.ny
    else:
        return

    x_center = intens[intens.shape[0]//2, :]
    fwhm_x = len(x_center[x_center > x_center.max()/2])*dx

    y_center = intens[:, intens.shape[1]//2]
    fwhm_y = len(y_center[y_center > y_center.max()/2])*dy
    if (wfr.params.wSpace == 'Q-space'):
        print(wfr.params.wSpace)
        wl = 12.39*1e-10/(wfr.params.photonEnergy*1e-3)  # WaveLength
        fwhm_x = fwhm_x*wl
        fwhm_y = fwhm_y*wl

    return {'fwhm_x': fwhm_x, 'fwhm_y': fwhm_y}


def get_intensity_on_axis(wfr):
    """
    Calculate intensity (e.g. spectrum in frequency domain) along (x=y=0)
    :param wfr:  wavefront
    :return: [z,s0] in [a.u.]
    """

    wf_intensity = wfr.get_intensity(polarization='horizontal')
    mesh = wfr.params.Mesh
    # array dimensions # <-to avoid wrong dimension assignment
    dim = numpy.shape(wf_intensity)
    sz = numpy.zeros(shape=(mesh.nSlices, 2), dtype='float64')
    sz[:, 0] = numpy.linspace(mesh.sliceMin, mesh.sliceMax, mesh.nSlices)
    # <-to avoid wrong dimension assignment
    sz[:, 1] = wf_intensity[dim[0]/2, dim[1]/2, :] / wf_intensity.max()

    return sz
