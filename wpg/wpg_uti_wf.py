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

    :param wf: wpg.Wavefront structure
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

    :param wf: wpg.Wavefront structure
    :return: pulse energy value, J
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


def averaged_intensity(wf, bPlot=False):
    """
    wrapper for integral_intensity() for backward compatibility

    """
    integral_intensity(wf, bPlot)

def integral_intensity(wf, threshold=0.01, bPlot=True):
    """
    plot the slice-to-slice integral intensity averaged over a meaningful range

    :param wf: wavefront structure
    :param threshold: defined the threshold for slices, integrated_slice_intensity_max*threshold
    :param bPlot: if True plot temporary structure or spectrum in the meaningful range
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold, mainly needed for processing spiky FEL source

    """
    J2eV = 6.24150934e18
    # total0=wf.get_intensity().sum();
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin)/(mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin)/(mesh.ny - 1)
    int0 = wf.get_intensity().sum(axis=0).sum(axis=0)  # I(slice_num)
    int0 = int0*(dx*dy*1.e6) # wf amplitude units sqrt(W/mm^2)

    # Get center pixel numbers.
    center_nx = int(mesh.nx/2)
    center_ny = int(mesh.ny/2)
    int0_00 = wf.get_intensity()[center_ny, center_nx, :]
    int0max = max(int0)

    # Get meaningful slices.
    aw = [a[0] for a in numpy.argwhere(int0 > int0max*threshold)]
    int0_mean = int0[min(aw):max(aw)]  # meaningful range of pulse
    if bPlot:
        dSlice = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices - 1)
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices)*dSlice+ mesh.sliceMin,int0)
        pylab.plot(numpy.arange(min(aw), max(aw))*dSlice + mesh.sliceMin, int0_mean, 'ro')
        if(wf.params.wDomain=='time'):
            pylab.title('Power');pylab.xlabel('s');pylab.ylabel('W')
        else: #frequency domain
            pylab.title('Spectral Energy');pylab.xlabel('eV');pylab.ylabel('J/eV')
        pylab.show()
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices)*dSlice+ mesh.sliceMin,int0_00)
        pylab.plot(numpy.arange(min(aw), max(aw))*dSlice + mesh.sliceMin, int0_00[min(aw):max(aw)], 'ro')
        if(wf.params.wDomain=='time'):
            pylab.title('On-Axis Power Density');pylab.xlabel('s');pylab.ylabel('W/mm^2')
        else: #frequency domain
            pylab.title('On-Axis Spectral Fluence');pylab.xlabel('eV');pylab.ylabel('J/eV/mm^2')
        pylab.show()
    averaged = int0_mean.sum()/len(int0_mean)
    print('number of meaningful slices:', len(int0_mean))
    if(wf.params.wDomain=='time'):
        dt = (mesh.sliceMax - mesh.sliceMin)/(mesh.nSlices - 1)
        print('Pulse energy {:1.2g} J'.format(int0_mean.sum()*dt))
    return averaged

def plot_t_wf(wf, save='', range_x=None, range_y=None,im_aspect='equal'):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map() for backward compatibility

    """
    integral_intensity(wf, bPlot=True)
    plot_intensity_map(wf, save, range_x, range_y,im_aspect)

def plot_wf(wf, save='', range_x=None, range_y=None,im_aspect='equal'):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map() for backward compatibility

    """
    integral_intensity(wf, bPlot=True)
    plot_intensity_map(wf, save, range_x, range_y,im_aspect)

def plot_intensity_map(wf, save='', range_x=None, range_y=None,im_aspect='equal'):
    """
    Plot wavefront in  R-space.

    :param wf: wavefront structure
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
    :param im_aspect: aspect for 2D image, string or float number, see matplotlib set_aspect().
    """
    import matplotlib.pyplot as plt
    # Get the wavefront and integrate over time.
    wf_intensity = wf.get_intensity().sum(axis=-1)

    # Get average and time slicing.
    #average = averaged_intensity(wf, bPlot=True)
    nslices = wf.params.Mesh.nSlices
    if (nslices>1):
        dt = (wf.params.Mesh.sliceMax-wf.params.Mesh.sliceMin)/(nslices-1)
        t0 = dt*nslices/2 + wf.params.Mesh.sliceMin
    else:
        t0 = (wf.params.Mesh.sliceMax+wf.params.Mesh.sliceMin)/2

    # Setup a figure.
    figure = plt.figure(figsize=(10, 10), dpi=100)
    plt.axis('tight')
    # Profile plot.
    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    # Get limits.
    xmin, xmax, ymax, ymin = wf.get_limits()

    # Plot profile as 2D colorcoded map.
    profile.imshow(
        wf_intensity, extent=[xmin*1.e3, xmax*1.e3, ymax*1.e3, ymin*1.e3], cmap="YlGnBu_r")
    profile.set_aspect(im_aspect, 'datalim')


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
    if(wf.params.wDomain=='time'):
        x_projection.set_title('t0={:03.1g} s '.format(t0))
    else: #frequency domain
        x_projection.set_title('E0={:05.2g} eV'.format(t0))

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


#def plot_t_wf_a(wf, save='', range_x=None, range_y=None):
    """
    Plot wavefront in Q-space.

    :param wf: wavefront structure
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

def look_at_q_space(wf, output_file=None, save='', range_x=None, range_y=None):
    """
    a wrapper for backward compatibility

    """
    plot_intensity_qmap(wf, output_file, save, range_x, range_y)

def plot_intensity_qmap(wf, output_file=None, save='', range_x=None, range_y=None,im_aspect='equal'):
    """
    change wavefront representation from R- to Q-space and plot it the resulting wavefront.

    :param wf: Wavefront object in R-space representation
    :param output_file: if parameter present - store wavefront in Q-space to the file
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
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

    #plot_t_wf_a(wf_a, save=save, range_x=range_x, range_y=range_y)
    import matplotlib.pyplot as plt
    wf_intensity = wf_a.get_intensity().sum(axis=-1)
    plt.figure(figsize=(10, 10), dpi=100)
    plt.axis('tight')

    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
    xmin, xmax, ymax, ymin = wf_a.get_limits()

    profile.imshow(
        wf_intensity, extent=[xmin*1.e6, xmax*1.e6, ymax*1.e6, ymin*1.e6], cmap="YlGnBu_r")
    profile.set_aspect(im_aspect, 'datalim')

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

    return


def propagate_wavefront(wavefront, beamline, output_file=None):
    """
    Propagate wavefront and store it in output file.

    :param wavefront: wpg.Wavefront object or path to HDF5 file
    :param beamline: SRWLOptC container of beamline
    :param output_file: if parameter present - store propagaed wavefront to file
    :return: propagated wavefront object
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
        wl = 12.398*1e-10/(wfr.params.photonEnergy*1e-3)  # WaveLength
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
