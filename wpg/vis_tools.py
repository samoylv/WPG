# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov'

import pylab as plt
import numpy as np
import wavefront


def get_limits2d(wf, axis):
    m = wf.srwl_wfr.mesh
    if axis == 'z':
        return [m.xStart, m.xFin, m.yStart, m.yFin]
    elif axis == 'x':
        return [m.eStart, m.eFin, m.yStart, m.yFin]
    elif axis == 'y':
        return [m.eStart, m.eFin, m.xStart, m.xFin]


def get_slice(data, slice_number, axis=None):
    """
    Return 2D slice from 3D array along given axis

    :param data:
    :param slice_numb:
    :param axis:
    :return: :raise:
    """
    if not axis in ['x', 'y', 'z']:
        raise ValueError("Axis parameter must be 'x','y' or 'z', but '{0}' given".format(axis))
    if not type(data) is np.ndarray:
        raise TypeError('Type of data must be numpy.ndarray, given type is {0}'.format(type(data)))
    if not len(data.shape) == 3:
        raise TypeError('Dimension of data must be 3, given dimension is {0}'.format(len(data.shape)))

    if slice_number == 'center':
        if axis == 'z':
            res = data[:, :, int(data.shape[2] / 2)]
        elif axis == 'x':
            res = data[int(data.shape[0] / 2), :, :]
        elif axis == 'y':
            res = data[:, int(data.shape[1] / 2), :]
    else:
        if axis == 'z':
            res = data[:, :, slice_number]
        elif axis == 'x':
            res = data[slice_number, :, :]
        elif axis == 'y':
            res = data[:, slice_number, :]

    return np.squeeze(res)


def show_intensity(data, slice_number=None, polarization=None, axis=None, extent=None, cmap=None, vmin=None, vmax=None):
    return show_2d_data(func=data.get_intensity, data=data, slice_number=slice_number,
                        polarization=polarization, axis=axis, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)


def show_phase(data, slice_number=None, polarization=None, axis=None, extent=None, cmap=None, vmin=None, vmax=None):
    return show_2d_data(func=data.get_phase, data=data, slice_number=slice_number,
                        polarization=polarization, axis=axis, extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)


def show_2d_data(func, data, slice_number=None, axis=None, extent=None, cmap=None, vmin=None, vmax=None,
                 polarization=None):
    if axis is None:
        axis = 'z'

    if type(data) is np.ndarray:
        if not len(data.shape) in [2, 3]:
            raise TypeError('Dimension of data must be 2 or 3, given dimension is {0}'.format(len(data.shape)))
        if len(data.shape) == 3:
            slice_data = get_slice(data, slice_number=slice_number, axis=axis)
    elif type(data) is wavefront.Wavefront:
        func_data = func(polarization=polarization)
        slice_data = get_slice(func_data, slice_number=slice_number, axis=axis)
        if extent is None:
            extent = data.get_limits(axis=axis)
    else:
        raise TypeError(
            'Type of data must be numpy.ndarray or wavefront.Wavefront, given type is {0}'.format(type(data))
        )

    plt.imshow(slice_data, cmap=cmap, vmin=vmin, vmax=vmax, extent=extent)

    if axis == 'z':
        plt.xlabel('m')
        plt.ylabel('m')
    elif axis == 'y':
        plt.xlabel('s')
        plt.ylabel('m')
    elif axis == 'y':
        plt.xlabel('m')
        plt.ylabel('s')
    else:
        pass

    if not axis == 'z':
        plt.axes().set_aspect('auto')