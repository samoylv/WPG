# -*- coding: utf-8 -*-

import copy
import numpy
import matplotlib.pyplot as pylab

from wpg.beamline import Beamline
from srwpy.srwlib import srwl
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
            wf_mesh.nx, wf_mesh.xMin * 1e3, wf_mesh.xMax * 1e3)
        )
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mm'.format(
            wf_mesh.ny, wf_mesh.yMin * 1e3, wf_mesh.yMax * 1e3)
        )
    if (w_space == 'Q-space'):
        print('nx {:5d}  range_x [{:.1e}, {:.1e}] mrad'.format(
            wf_mesh.nx, wf_mesh.qxMin * 1e3, wf_mesh.qxMax * 1e3)
        )
        print('ny {:5d}  range_y [{:.1e}, {:.1e}] mrad'.format(
            wf_mesh.ny, wf_mesh.qyMin * 1e3, wf_mesh.qyMax * 1e3)
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
        pulse_energy_J = pulse_energy * dx * dy * 1e6 * dt
        print('Number of photons per pulse: {:e}'.format(
            pulse_energy_J * J2eV / wf.params.photonEnergy))
        return pulse_energy_J


def averaged_intensity(wf, bPlot=False):
    """
    wrapper for integral_intensity() for backward compatibility

    """
    return integral_intensity(wf, bPlot=bPlot)


def integral_intensity(wf, threshold=0.01, bPlot=True, intensity=None):
    """
    plot the slice-to-slice integral intensity averaged over a meaningful range

    :param wf: wavefront structure
    :param threshold: defined the threshold for slices, integrated_slice_intensity_max*threshold
    :param bPlot: if True plot temporary structure or spectrum in the meaningful range
    :param intensity: wf.get_intensity() if already computed (computing it takes seconds for a long pulse)
    :return: intensity averaged over 'meaningful' slices, i.e. above 1% threshold, mainly needed for processing spiky FEL source

    """
    J2eV = 6.24150934e18
    # total0=wf.get_intensity().sum();
    mesh = wf.params.Mesh
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    if intensity is None:
        intensity = wf.get_intensity()
    int0 = intensity.sum(axis=0).sum(axis=0)  # I(slice_num)
    int0 = int0 * (dx * dy * 1.e6)  # wf amplitude units sqrt(W/mm^2)

    # Get center pixel numbers.
    center_nx = int(mesh.nx / 2)
    center_ny = int(mesh.ny / 2)
    int0_00 = intensity[center_ny, center_nx, :]
    int0max = max(int0)

    # Get meaningful slices.
    aw = [a[0] for a in numpy.argwhere(int0 > int0max * threshold)]
    int0_mean = int0[min(aw):max(aw) + 1]  # meaningful range of pulse
    if bPlot:
        if mesh.nSlices > 1:
            dSlice = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        else:
            dSlice = 0
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices) * dSlice + mesh.sliceMin, int0)
        pylab.plot(numpy.arange(min(aw), max(aw) + 1) *
                   dSlice + mesh.sliceMin, int0_mean, 'ro')
        if(wf.params.wDomain == 'time'):
            pylab.title('Power')
            pylab.xlabel('s')
            pylab.ylabel('W')
        else:  # frequency domain
            pylab.title('Spectral Energy')
            pylab.xlabel('eV')
            pylab.ylabel('J/eV')
        pylab.show()
        pylab.figure()
        pylab.plot(numpy.arange(mesh.nSlices) *
                   dSlice + mesh.sliceMin, int0_00)
        pylab.plot(numpy.arange(min(aw), max(aw) + 1) * dSlice +
                   mesh.sliceMin, int0_00[min(aw):max(aw) + 1], 'ro')
        if(wf.params.wDomain == 'time'):
            pylab.title('On-Axis Power Density')
            pylab.xlabel('s')
            pylab.ylabel('W/mm^2')
        else:  # frequency domain
            pylab.title('On-Axis Spectral Fluence')
            pylab.xlabel('eV')
            pylab.ylabel('J/eV/mm^2')
        pylab.show()
    averaged = int0_mean.sum() / len(int0_mean)
    print('number of meaningful slices:', len(int0_mean))
    if(wf.params.wDomain == 'time') and mesh.nSlices > 1:
        dt = (mesh.sliceMax - mesh.sliceMin) / (mesh.nSlices - 1)
        print('Pulse energy {:1.2g} J'.format(int0_mean.sum() * dt))
    return averaged


def plot_t_wf(wf, save='', range_x=None, range_y=None, im_aspect='equal'):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map()
    for backward compatibility

    """
    plot_wf(wf, save, range_x, range_y, im_aspect)


def plot_wf(wf, save='', range_x=None, range_y=None, im_aspect='equal'):
    """
    a wrapper, calls integral_intensity() and plot_intensity_map()
    for backward compatibility; the intensity is computed once for both

    """
    intensity = wf.get_intensity()
    integral_intensity(wf, bPlot=True, intensity=intensity)
    plot_intensity_map(wf, save, range_x, range_y, im_aspect, intensity=intensity)


def plot_intensity_map(wf, save='', range_x=None, range_y=None,
                       im_aspect='equal', intensity=None):
    """
    Plot wavefront in  R-space.

    :param wf: wavefront structure
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.
    :param im_aspect: aspect for 2D image, string or float number, see matplotlib set_aspect().
    :param intensity: wf.get_intensity() if already computed
    """
    import matplotlib.pyplot as plt
    if intensity is None:
        intensity = wf.get_intensity()
    # Get the wavefront and integrate over time.
    wf_intensity = intensity.sum(axis=-1)

    # Get average and time slicing.
    # average = averaged_intensity(wf, bPlot=True)
    nslices = wf.params.Mesh.nSlices
    if (nslices > 1):
        dt = (wf.params.Mesh.sliceMax - wf.params.Mesh.sliceMin) / (nslices - 1)
        t0 = (wf.params.Mesh.sliceMax + wf.params.Mesh.sliceMin) / 2
    else:
        t0 = (wf.params.Mesh.sliceMax + wf.params.Mesh.sliceMin) / 2

    # Setup a figure.
    plt.figure(figsize=(10, 10), dpi=100)
    # Profile plot.
    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    # Get limits.
    xmin, xmax, ymax, ymin = wf.get_limits()

    # Plot profile as 2D colorcoded map.
    profile.imshow(
        wf_intensity, extent=[xmin * 1.e3, xmax *
                              1.e3, ymax * 1.e3, ymin * 1.e3],
        cmap="YlGnBu_r")
    profile.set_aspect(im_aspect)

    # Get x and y ranges.
    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    #x = numpy.linspace(xmin*1.e3,xmax*1.e3,wf_intensity.shape[0])
    #y = numpy.linspace(ymin*1.e3,ymax*1.e3,wf_intensity.shape[1])
    x = numpy.linspace(xmin * 1.e3, xmax * 1.e3, wf_intensity.shape[1])
    y = numpy.linspace(ymin * 1.e3, ymax * 1.e3, wf_intensity.shape[0])

    # Labels.
    profile.set_xlabel('$mm$', fontsize=12)
    profile.set_ylabel('$mm$', fontsize=12)

    # x-projection plots above main plot.
    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    print(x.shape, wf_intensity.sum(axis=0).shape)
    fwhm = calculate_fwhm(wf, intensity=intensity)

    print("FWHM in x = %4.3e m." % (fwhm['fwhm_x']))
    print("FWHM in y = %4.3e m." % (fwhm['fwhm_y']))

    x_projection.plot(x, wf_intensity.sum(axis=0), label='x projection')

    # Set range according to input.
    if range_x is None:
        profile.set_xlim([xmin * 1.e3, xmax * 1.e3])
    else:
        profile.set_xlim([-range_x / 2., range_x / 2.])

    # Set title.
    if(wf.params.wDomain == 'time'):
        x_projection.set_title('t = {:.4g} fs'.format(t0 * 1e15))
    else:  # frequency domain
        x_projection.set_title('E = {:.3f} eV'.format(t0))

    # y-projection plot right of main plot.
    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label='y projection')

    # Hide minor tick labels, they disturb here.
    plt.minorticks_off()

    # Set range according to input.
    if range_y is None:
        profile.set_ylim([ymin * 1.e3, ymax * 1.e3])
    else:
        profile.set_ylim([-range_y / 2., range_y / 2.])

    # If requested, save to disk, otherwise show in interactive window.
    if save != '':
        # Add parameters.
        plt.savefig(save)
    else:
        plt.show()


# def plot_t_wf_a(wf, save='', range_x=None, range_y=None):
#     """
#     Plot wavefront in Q-space.

#     :param wf: wavefront structure
#     :params: save: Whether to save the figure on disk
#     :type:  string for filename. Empty string '' means don't save.
#     :default: '', do not save the figure.

#     :params: range_x: x-axis range.
#     :type: float
#     :default: None, take entire x range.

#     :params: range_y: y-ayis range.
#     :type: float
#     :default: None, take entire y range.

#     """

def _fft_workers():
    """Number of CPU cores available to this process (for scipy.fft)."""
    import os
    try:
        return len(os.sched_getaffinity(0))
    except AttributeError:  # Windows, macOS
        return os.cpu_count() or 1


def angular_intensity(wf, chunk=64):
    """
    Intensity in angles summed over the slices, by FFT of the field of every slice (scipy.fft, all cores).

    The same as the intensity after ``srwl.SetRepresElecField(wf, 'a')`` summed over the slices, in the same
    units, but without a copy of the wavefront and much faster: SRW converts the slices one by one in one
    thread (6.5 s for 201 x 201 x 2662 slices with 1 and with 8 OpenMP threads). Checked against SRW for even
    nx, ny (the same mesh, normalized maps agree within 1e-6). For odd nx or ny SRW takes the angular step
    1/((n - 1) dx) instead of 1/(n dx) of the discrete Fourier transform and gives NaN at some points, here
    the step is 1/(n dx).

    :param wf: Wavefront in R-space
    :param chunk: number of slices transformed at once (limits the memory)
    :return: (intensity (ny, nx), limits (xmin, xmax, ymax, ymin) [rad], as Wavefront.get_limits in Q-space)
    """
    import scipy.fft
    m = wf.params.Mesh
    nx, ny, ns = m.nx, m.ny, m.nSlices
    dx, dy = (m.xMax - m.xMin) / (nx - 1), (m.yMax - m.yMin) / (ny - 1)
    workers = _fft_workers()
    total = numpy.zeros((ny, nx))
    for name in ('arrEhor', 'arrEver'):
        arr = getattr(wf.data, name)
        if arr is None:
            continue
        for k0 in range(0, ns, chunk):
            part = arr[:, :, k0:k0 + chunk, :]
            e = numpy.empty(part.shape[:3], dtype=numpy.complex64)
            e.real, e.imag = part[..., 0], part[..., 1]
            e = scipy.fft.fft2(e, axes=(0, 1), workers=workers, overwrite_x=True)
            total += (e.real ** 2 + e.imag ** 2).sum(axis=-1)
    # the continuous Fourier transform: the field times dx dy, as SRW
    total = numpy.fft.fftshift(total) * (dx * dy) ** 2
    wl = wf.wavelength
    fx = numpy.fft.fftshift(numpy.fft.fftfreq(nx, dx))
    fy = numpy.fft.fftshift(numpy.fft.fftfreq(ny, dy))
    return total, (fx[0] * wl, fx[-1] * wl, fy[-1] * wl, fy[0] * wl)


def look_at_q_space(wf, output_file=None, save='',
                    range_x=None, range_y=None):
    """
    a wrapper for backward compatibility

    """
    plot_intensity_qmap(wf, output_file, save, range_x, range_y)


def plot_intensity_qmap(wf, output_file=None, save='',
                        range_x=None, range_y=None, im_aspect='equal'):
    """
    change wavefront representation from R- to Q-space and plot it the resulting wavefront.

    :param wf: Wavefront object in R-space representation
    :param output_file: if parameter present - store wavefront in Q-space to the file
    :param save: string for filename. Empty string '' means don't save.
    :param range_x: x-axis range, _float_. If None, take entire x range.
    :param range_y: y-ayis range, float. If None, take entire y range.

    The map is computed by :func:`angular_intensity` (FFT on all cores, no copy of the wavefront). Only with
    ``output_file`` the wavefront is copied and converted by SRW, to store it in Q-space.
    """
    if not wf.params.wSpace == 'R-space':
        print('space should be in R-space, but not ' + wf.params.wSpace)
        return
    if output_file is not None:
        srwl_wf_a = copy.deepcopy(wf._srwl_wf)
        srwl.SetRepresElecField(srwl_wf_a, 'a')
        print('store wavefront to HDF5 file: ' + output_file+'...')
        Wavefront(srwl_wf_a).store_hdf5(output_file)
        print('done')
        del srwl_wf_a

    wf_intensity, (xmin, xmax, ymax, ymin) = angular_intensity(wf)
    # FWHM through the maximum, as calculate_fwhm in Q-space
    iy, ix = numpy.unravel_index(numpy.argmax(wf_intensity), wf_intensity.shape)
    row, col = wf_intensity[iy, :], wf_intensity[:, ix]
    print({'fwhm_x': float(numpy.count_nonzero(row > row.max() / 2) * (xmax - xmin) / (wf_intensity.shape[1] - 1)),
           'fwhm_y': float(numpy.count_nonzero(col > col.max() / 2) * (ymax - ymin) / (wf_intensity.shape[0] - 1))})

    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 10), dpi=100)

    profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)

    profile.imshow(
        wf_intensity, extent=[xmin * 1.e6, xmax *
                              1.e6, ymax * 1.e6, ymin * 1.e6],
        cmap="YlGnBu_r")
    profile.set_aspect(im_aspect)

    # [LS:2016-03-17]
    # change shape dimension, otherwise, in case nx!=ny ,
    # 'x, y should have the same dimension' error from py plot
    # x = numpy.linspace(xmin*1.e6,xmax*1.e6,wf_intensity.shape[0])
    # y = numpy.linspace(ymin*1.e6,ymax*1.e6,wf_intensity.shape[1])
    x = numpy.linspace(xmin * 1.e6, xmax * 1.e6, wf_intensity.shape[1])
    y = numpy.linspace(ymin * 1.e6, ymax * 1.e6, wf_intensity.shape[0])
    profile.set_xlabel(r'$\mu$rad', fontsize=12)
    profile.set_ylabel(r'$\mu$rad', fontsize=12)

    x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
    print(x.shape, wf_intensity.sum(axis=0).shape)
    x_projection.plot(x, wf_intensity.sum(axis=0), label='x projection')
    if range_x is None:
        profile.set_xlim([xmin * 1.e6, xmax * 1.e6])
    else:
        profile.set_xlim([-range_x / 2., range_x / 2.])

    y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
    y_projection.plot(wf_intensity.sum(axis=1), y, label='y projection')

    # Hide minor tick labels.
    plt.minorticks_off()

    if range_y is None:
        profile.set_ylim([ymin * 1.e6, ymax * 1.e6])
    else:
        profile.set_ylim([-range_y / 2., range_y / 2.])

    if save != '':
        plt.savefig(save)
    else:
        plt.show()

    return


def plot_slices(wf, threshold=0.2, polarization='total', save='', range_x=None, range_y=None,
                im_aspect='equal'):
    """
    Plot the intensity maps of the slices whose power is above ``threshold`` of the most powerful slice.

    The axes follow the representation of the wavefront: coordinates in mm (R-space) or angles in urad
    (Q-space, e.g. a copy after ``srwl.SetRepresElecField(wf._srwl_wf, 'a')``; use an even number of points:
    with an odd nx or ny SRW returns a wrong field, https://github.com/ochubar/SRW/issues/70). The title of a
    map is the time of the slice [fs] or its photon energy [eV]. The intensity is computed once for all slices
    and the wavefront is not copied. From the slice plots of L. Samoylova.

    :param wf: Wavefront
    :param threshold: power of a slice (intensity summed over x, y) relative to the most powerful one
    :param polarization: 'total', 'horizontal' or 'vertical'
    :param save: '' - show the figures, else the prefix of the files: <save><slice number>.png
    :param range_x: x-axis range [mm in R-space, urad in Q-space], None - the whole mesh
    :param range_y: y-axis range, as range_x
    :param im_aspect: aspect of the map, see matplotlib set_aspect()
    :return: dict: 'selected_slices' (numbers), 'slice_power' (relative to the maximum), 'slice_values'
        (time [s] or photon energy [eV] of the selected slices)
    """
    import warnings
    import matplotlib.pyplot as plt
    m = wf.params.Mesh
    qspace = wf.params.wSpace == 'Q-space'
    if qspace and (m.nx % 2 or m.ny % 2):
        warnings.warn('plot_slices: {} x {} points in Q-space: SRW converts a mesh with an odd number of points '
                      'to angles wrongly (https://github.com/ochubar/SRW/issues/70)'.format(m.nx, m.ny))
    intensity = wf.get_intensity(polarization=polarization)
    power = intensity.sum(axis=(0, 1))
    power = power / power.max()
    selected = numpy.flatnonzero(power >= threshold)
    values = numpy.linspace(m.sliceMin, m.sliceMax, m.nSlices)
    xmin, xmax, ymax, ymin = wf.get_limits()  # [m] in R-space, [rad] in Q-space
    scale, unit = (1e6, r'$\mu$rad') if qspace else (1e3, 'mm')
    x = numpy.linspace(xmin, xmax, m.nx) * scale
    y = numpy.linspace(ymin, ymax, m.ny) * scale
    for k in selected:
        cut = intensity[:, :, k]
        plt.figure(figsize=(10, 10), dpi=100)
        profile = plt.subplot2grid((3, 3), (1, 0), colspan=2, rowspan=2)
        profile.imshow(cut, extent=[x[0], x[-1], y[-1], y[0]], cmap="YlGnBu_r")
        profile.set_aspect(im_aspect)
        profile.set_xlabel(unit, fontsize=12)
        profile.set_ylabel(unit, fontsize=12)
        x_projection = plt.subplot2grid((3, 3), (0, 0), sharex=profile, colspan=2)
        x_projection.plot(x, cut.sum(axis=0))
        if wf.params.wDomain == 'time':
            x_projection.set_title('slice {}: t = {:.4g} fs'.format(k, values[k] * 1e15))
        else:
            x_projection.set_title('slice {}: E = {:.3f} eV'.format(k, values[k]))
        y_projection = plt.subplot2grid((3, 3), (1, 2), rowspan=2, sharey=profile)
        y_projection.plot(cut.sum(axis=1), y)
        plt.minorticks_off()
        profile.set_xlim([x[0], x[-1]] if range_x is None else [-range_x / 2., range_x / 2.])
        profile.set_ylim([y[0], y[-1]] if range_y is None else [-range_y / 2., range_y / 2.])
        # FWHM through the maximum of the slice, as calculate_fwhm
        iy, ix = numpy.unravel_index(numpy.argmax(cut), cut.shape)
        fwhm = [numpy.count_nonzero(line > line.max() / 2) * step
                for line, step in ((cut[iy, :], x[1] - x[0]), (cut[:, ix], y[1] - y[0]))]
        print('slice {}: power {:.3f} of the maximum, FWHM {:.4g} x {:.4g} {}'.format(
            k, power[k], fwhm[0], fwhm[1], 'urad' if qspace else 'mm'))
        if save != '':
            plt.savefig('{}{}.png'.format(save, k))
            plt.close()
        else:
            plt.show()
    return {'selected_slices': selected, 'slice_power': power[selected], 'slice_values': values[selected]}


def propagate_wavefront(wavefront, beamline, output_file=None, gpu=None):
    """
    Propagate wavefront and store it in output file.

    :param wavefront: wpg.Wavefront object or path to HDF5 file
    :param beamline: SRWLOptC container of beamline
    :param output_file: if parameter present - store propagaed wavefront to file
    :param gpu: number of the GPU device (1 - first GPU), None - CPU, see Beamline.propagate
    :return: propagated wavefront object
    """

    if not isinstance(beamline, Beamline):
        bl = Beamline(beamline)
    else:
        bl = beamline
    print(bl)

    if isinstance(wavefront, Wavefront):
        # copy: the caller's wavefront must not be propagated in place
        wfr = Wavefront(srwl_wavefront=copy.deepcopy(wavefront._srwl_wf))
    else:
        print('*****reading wavefront from h5 file...')
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)

    print_mesh(wfr)
    print('*****propagating wavefront (with resizing)...')
    bl.propagate(wfr, gpu=gpu)

    if output_file is not None:
        print('save hdf5:', output_file)
        wfr.store_hdf5(output_file)
    print('done')
    return wfr


def calculate_fwhm(wfr, intensity=None):
    """
    Calculate FWHM of the beam calculating number of point bigger then max / 2 throuhgt center of the image

    :param wfr:  wavefront
    :param intensity: wfr.get_intensity() (total polarization) if already computed
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m] in R-space, in [rad] in Q-space
    """
    if intensity is None:
        intensity = wfr.get_intensity(polarization='total')
    intens = intensity.sum(axis=-1)

    mesh = wfr.params.Mesh
    if (wfr.params.wSpace == 'R-space'):
        dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)
        dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)
    elif (wfr.params.wSpace == 'Q-space'):
        dx = (mesh.qxMax - mesh.qxMin) / (mesh.nx - 1)
        dy = (mesh.qyMax - mesh.qyMin) / (mesh.ny - 1)
    else:
        return

    # row and column through the maximum (the centre of the array missed an off-axis beam)
    iy, ix = numpy.unravel_index(numpy.argmax(intens), intens.shape)
    x_center = intens[iy, :]
    fwhm_x = len(x_center[x_center > x_center.max() / 2]) * dx

    y_center = intens[:, ix]
    fwhm_y = len(y_center[y_center > y_center.max() / 2]) * dy
    if (wfr.params.wSpace == 'Q-space'):
        wl = 12.398 * 1e-10 / (wfr.params.photonEnergy * 1e-3)  # WaveLength
        fwhm_x = fwhm_x * wl
        fwhm_y = fwhm_y * wl

    return {'fwhm_x': fwhm_x, 'fwhm_y': fwhm_y}


def get_intensity_on_axis(wfr):
    """
    Calculate intensity (e.g. spectrum in frequency domain) along (x=y=0)

    :param wfr:  wavefront
    :return: [z,s0] in [a.u.]
    """

    wf_intensity = wfr.get_intensity(polarization='total')
    mesh = wfr.params.Mesh
    # array dimensions # <-to avoid wrong dimension assignment
    dim = numpy.shape(wf_intensity)
    sz = numpy.zeros(shape=(mesh.nSlices, 2), dtype='float64')
    sz[:, 0] = numpy.linspace(mesh.sliceMin, mesh.sliceMax, mesh.nSlices)
    # <-to avoid wrong dimension assignment
    on_axis = wf_intensity[dim[0] // 2, dim[1] // 2, :]
    sz[:, 1] = on_axis / on_axis.max()

    return sz

def check_sampling(wavefront):
    """
    Report on the sampling of a wavefront (text), see :func:`wpg.sampling.check_sampling`.

    The former criterion (Fresnel zone lambda*R/(2*FWHM) within 7-10 pixels) gave false alarms: a
    negative R of a converging beam never passed, R of SRW near the waist is the distance to it, and
    oversampling was reported as an error. The phase step between neighbouring pixels is measured now.
    """
    from wpg.sampling import check_sampling as _check
    return _check(wavefront)

    