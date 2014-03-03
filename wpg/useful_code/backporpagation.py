import os
import errno
import gc
import pylab as plt
import numpy as np

from wpg.srwlib import SRWLOptC
from wpg import Wavefront
from wpg import Beamline # new class for fixing bugs of srw classes

from wpg import optical_elements 

def mkdir_p(path):
    """
    Create directory tree, if not exists (mkdir -p)
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def gaussian(height, center_x, center_y, width_x, width_y):
    """Returns a gaussian function with the given parameters"""
    width_x = float(width_x)
    width_y = float(width_y)
    return lambda x, y: height * np.exp(
        -(((center_x - x) / width_x) ** 2 + ((center_y - y) / width_y) ** 2) / 2)


def moments(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution by calculating its
    moments """
    total = data.sum()
    X, Y = np.indices(data.shape)
    x = (X * data).sum() / total
    y = (Y * data).sum() / total
    col = data[:, int(y)]
    width_x = np.sqrt(
        abs((np.arange(col.size) - y) ** 2 * col).sum() / col.sum())
    row = data[int(x), :]
    width_y = np.sqrt(
        abs((np.arange(row.size) - x) ** 2 * row).sum() / row.sum())
    height = data.max()
    return height, x, y, width_x, width_y


def fit_gaussian(data):
    """Returns (height, x, y, width_x, width_y)
    the gaussian parameters of a 2D distribution found by a fit"""
    params = moments(data)
    errorfunction = lambda p: np.ravel(
        gaussian(*p)(*np.indices(data.shape)) - data)
    x, cov_x, infodict, mesg, ier = scipy.optimize.leastsq(
        errorfunction, params, full_output=True)

    res = gaussian(*x)(*np.indices(data.shape))
    error = infodict['fvec']
    ss_err = ((error) ** 2).sum()
    ss_tot = ((res - res.mean()) ** 2).sum()
    rsquared = 1 - (ss_err / ss_tot)
    return {'params': x, 'rsquared': rsquared}


def fit_gaussian_pulse(wf, polarization='vertical'):
    """
    Calculate gaussian parameters for all slices in wavefront and add it to custom_fields['misc']['gaussain_parameters']
    """
    # TODO: add other polarizations
    wf_intensity = wf.get_intensity(None, polarization)
    slices_number = wf.params.Mesh.nSlices

    height = np.zeros(slices_number)
    center_x = np.zeros(slices_number)
    center_y = np.zeros(slices_number)
    width_x = np.zeros(slices_number)
    width_y = np.zeros(slices_number)
    rsquared = np.zeros(slices_number)
    slices_intensity = wf_intensity.sum(0).sum(0)
    pulse_intensity = wf.get_intensity().sum(axis=-1)

    # gaussinan fitting of slices
    for sn in range(slices_number):
        data = np.squeeze(wf_intensity[:, :, sn])
        fit_result = fit_gaussian(data)
        (height[sn], center_x[sn], center_y[sn],
         width_x[sn], width_y[sn]) = fit_result['params']
        rsquared[sn] = fit_result['rsquared']
        # TODO: normalize to real data
    del data
    del fit_result

    # gaussian fitting of total intensity
    data_pulse = pulse_intensity
    fit_result_pulse = fit_gaussian(data_pulse)
    (height_pulse, center_x_pulse, center_y_pulse,
     width_x_pulse, width_y_pulse) = fit_result_pulse['params']
    rsquared_pulse = fit_result_pulse['rsquared']

    fit_pulse = gaussian(height_pulse, center_x_pulse,
                         center_y_pulse, width_x_pulse, width_y_pulse)
    fit_data_pulse = fit_pulse(*np.indices(data_pulse.shape))

    # intensity inside gaussian fitting
    pulse_total_intensity_0_25 = np.sum(
        pulse_intensity * (fit_data_pulse > height_pulse * 0.25))
    pulse_total_intensity_0_50 = np.sum(
        pulse_intensity * (fit_data_pulse > height_pulse * 0.50))
    pulse_total_intensity_0_75 = np.sum(
        pulse_intensity * (fit_data_pulse > height_pulse * 0.75))

    if not 'misc' in wf.custom_fields:
        wf.custom_fields['misc'] = {}

    # TODO: check wf.id
    # TODO: custom fields to defaultdict? def tree(): return defaultdict(tree)
    if not 'gaussian_parameters' in wf.custom_fields['misc']:
        wf.custom_fields['misc']['gaussian_parameters'] = {}

    wf.custom_fields['misc']['gaussian_parameters']['height'] = height
    wf.custom_fields['misc']['gaussian_parameters']['center_x'] = center_x
    wf.custom_fields['misc']['gaussian_parameters']['center_y'] = center_y
    wf.custom_fields['misc']['gaussian_parameters']['width_x'] = width_x
    wf.custom_fields['misc']['gaussian_parameters']['width_y'] = width_y
    wf.custom_fields['misc']['gaussian_parameters']['rsquared'] = rsquared
    wf.custom_fields['misc']['gaussian_parameters'][
        'slices_intensity'] = slices_intensity

    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_height'] = height_pulse
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_center_x'] = center_x_pulse
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_center_y'] = center_y_pulse
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_width_x'] = width_x_pulse
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_width_y'] = width_y_pulse
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_rsquared'] = rsquared_pulse

    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_intensity'] = pulse_intensity
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_total_intensity_0_25'] = pulse_total_intensity_0_25
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_total_intensity_0_50'] = pulse_total_intensity_0_50
    wf.custom_fields['misc']['gaussian_parameters'][
        'pulse_total_intensity_0_75'] = pulse_total_intensity_0_75


def show_slices(wfr, slice_numbers=None):
    """
    Show slices: intensity, phase, gaussian approximation parameters and cuts.
    All gaussina parameters in pixels now. Should be fixed.
    
    :params wfr: wpg.Wavefront
    :params slice_numbers: slices to be shown, may by list, int, or None (for all slices)
    """
    wf_intensity = wfr.get_intensity(polarization='vertical')
    wf_phase = wfr.get_phase(polarization='vertical')

#     print wf_intensity.shape

    if slice_numbers is None:
        slice_numbers = range(wf_intensity.shape[-1])

    if isinstance(slice_numbers, int):
        slice_numbers = [slice_numbers, ]

    intense = wf_intensity.sum(0).sum(0)
    intense = np.squeeze(intense)
#     print intense.shape
    print 'Z coord: {0:.4f} m.'.format(wfr.params.Mesh.zCoord)
    plt.figure()
    plt.plot(slice_numbers, intense[slice_numbers])
    plt.title('Slices intensity')
    plt.show()

    total_intensity = wf_intensity.sum(axis=-1)
    data = total_intensity
    fit_result = fit_gaussian(data)
    fit_result = fit_gaussian(data)
    (height, center_x, center_y, width_x, width_y) = fit_result['params']
    rsquared = fit_result['rsquared']
    fit = gaussian(height, center_x, center_y, width_x, width_y)
    fit_data = fit(*np.indices(data.shape))

    print 'Total pulse intinsity'
    print '''Gaussian approximation parameters:
            center_x : {0:.1f}px.\t center_y : {1:.1f}px.
            width_x : {2:.1f}px\twidth_y : {3:.1f}px.
            rsquared: {4:0.4f}.'''.format(center_x, center_y, width_x, width_y, rsquared)

    plt.figure(figsize=(20, 7))
    plt.subplot(131)
    plt.imshow(data)
    plt.colorbar(orientation='horizontal')

    plt.contour(fit_data, cmap=plt.cm.copper)

    plt.subplot(132)
    plt.plot(data[:, int(center_x)], label='Data')
    plt.hold(True)
    plt.plot(fit_data[:, int(center_x)], label='Gaussian fit')
    plt.legend()

    plt.subplot(133)
    plt.plot(data[int(center_y), :], label='Data')
    plt.hold(True)
    plt.plot(fit_data[int(center_y), :], label='Gaussian fit')
    plt.legend()

    plt.show()

    for sn in slice_numbers:
        data = wf_intensity[:, :, sn]
        fit_result = fit_gaussian(data)
        fit_result = fit_gaussian(data)
        (height, center_x, center_y, width_x, width_y) = fit_result['params']
        rsquared = fit_result['rsquared']
        fit = gaussian(height, center_x, center_y, width_x, width_y)
        fit_data = fit(*np.indices(data.shape))

        print 'Slice number: {}'.format(sn)
        print '''Gaussian approximation parameters:
                center_x : {0:.1f}px.\t center_y : {1:.1f}px.
                width_x : {2:.1f}px\twidth_y : {3:.1f}px.
                rsquared: {4:0.4f}.'''.format(center_x, center_y, width_x, width_y, rsquared)

        plt.figure(figsize=(20, 7))
        plt.subplot(141)
        plt.imshow(data)
        plt.colorbar(orientation='horizontal')

        plt.contour(fit_data, cmap=plt.cm.copper)

        plt.subplot(142)
        plt.imshow(wf_phase[:, :, sn])

        plt.colorbar(orientation='horizontal')

        plt.subplot(143)
        plt.plot(data[:, int(center_x)], label='Data')
        plt.hold(True)
        plt.plot(fit_data[:, int(center_x)], label='Gaussian fit')
        plt.legend()

        plt.subplot(144)
        plt.plot(data[int(center_y), :], label='Data')
        plt.hold(True)
        plt.plot(fit_data[int(center_y), :], label='Gaussian fit')
        plt.legend()

        plt.show()


def forward_propagate(root_dir, distance, propagation_parameters):
    """
    Forward_propagate_wavefront
    the result will saved in root_dir\distance\distance.h5 file
    
    :param root_dir: directory, where '0.h' file located
    :param distance: distance to forward propagate initial wvefront
    :param propagation_parameters: SRW propagation parameters
    """

    out_dir = os.path.join(root_dir, '{:0.4f}'.format(distance))
    mkdir_p(out_dir)

    out_file_name = '{:0.4f}.h5'.format(distance)
    out_path = os.path.join(out_dir, out_file_name)

    if os.path.exists(out_path):
        print 'File exists: {}. Skiping.'.format(out_path)
        return out_path

    ppDrift0 = propagation_parameters

    drift0 = optical_elements.Drift(distance)
    srwl_bl0 = SRWLOptC([drift0, ], [ppDrift0, ])
    bl0 = Beamline(srwl_bl0)

    # forward propagate to L0 meters
    wf_L0 = Wavefront()
    wf_L0.load_hdf5(os.path.join(root_dir, '0.h5'))

    wpg.srwlib.srwl.SetRepresElecField(wf_L0._srwl_wf, 'f')
    bl0.propagate(wf_L0)
    wpg.srwlib.srwl.SetRepresElecField(wf_L0._srwl_wf, 't')
    fit_gaussian_pulse(wf_L0)
    wf_L0.store_hdf5(out_path)

    print 'Save file : {}'.format(out_path)

    del wf_L0
    return out_path


def back_propagate(params):
    '''
    Propagate pulse from file params[0] at the distance params[1] and save result to HDF5 file.
    If output files exists - skip calculations.    
    '''
    (input_path, distance, propagation_parameters) = params
    input_dir, input_file_name = os.path.split(input_path)

    out_file_name = '{}_{:0.4f}.h5'.format(
        '.'.join(input_file_name.split('.')[:-1]), distance)
    out_path = os.path.join(input_dir, out_file_name)

    if os.path.exists(out_path):
        return

    wf_L1 = Wavefront()
    wf_L1.load_hdf5(input_path)

    drift1 = optical_elements.Drift(distance)

    srwl_bl1 = SRWLOptC([drift1, ], [propagation_parameters, ])
    bl1 = Beamline(srwl_bl1)

    wpg.srwlib.srwl.SetRepresElecField(wf_L1._srwl_wf, 'f')
    bl1.propagate(wf_L1)
    wpg.srwlib.srwl.SetRepresElecField(wf_L1._srwl_wf, 't')

    fit_gaussian_pulse(wf_L1)
    wf_L1.store_hdf5(out_path)
    del wf_L1

    gc.collect()

    return out_path
