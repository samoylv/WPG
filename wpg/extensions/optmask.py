import math
import os
import scipy


import os
import sys

sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations")
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")

sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations")
sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")
sys.path.append(r"C:\Users\daemo\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")


from wpg.uti_io import *
from wpg.srwlib import SRWLOptT as SRWLOPT
from PIL import Image


import numpy as np

class np_to_sample:
    """The class for Samples from image file (e.g., .tif), NumPy array (.npy), etc.

    :param file_path: full path to the image or the saved NumPy array.
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_show_images: a flag to show the initial and processed images.
    :param is_save_images: a flag to save the initial and processed images.
    :param raw_image_name: the name of the raw file in case if it's saved.
    :param processed_image_name: the name of the processed file in case if it's saved.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    """

    def __init__(self, array,
                 area=None, rotate_angle=0, rotate_reshape=True, cutoff_background_noise=0.25, background_color=0,
                 tile=None, shift_x=None, shift_y=None, invert=None,
                 is_show_images=False):
        # Input parameters:
        self.array = array
        self.area = area
        self.rotate_angle = rotate_angle if rotate_angle is not None else 0
        self.rotate_reshape = True if rotate_reshape else False
        self.cutoff_background_noise = cutoff_background_noise if cutoff_background_noise is not None else 0
        self.background_color = background_color if background_color is not None else 0
        self.invert = invert
        self.tile = tile
        self.shift_x = shift_x
        self.shift_y = shift_y
        self.is_show_images = is_show_images
        

        # Output parameters:
        self.data = None
        self.raw_image = None
        self.processed_image = None
        self.nx = 128
        self.ny = 128
        self.limit_value = 1

        # Show the resulted images:
        if self.is_show_images:
            self.show_images()






        if self.tile:
            assert type(self.tile) in [list, tuple], 'The type of tile ({}) should be list/tuple'.format(
                type(self.tile).__name__)
            assert len(self.tile) == 2, 'The size ({}) of the list/tuple should be 2'.format(len(self.tile))
            self.array = np.tile(self.array, self.tile)

        if self.rotate_angle:
            assert -360 < self.rotate_angle < 360, 'The angle should be from -360 to 360 degrees.'
            self.array = rotate(self.array, self.rotate_angle, axes=(1, 0), reshape=self.rotate_reshape, output=None,
                               order=0, mode='constant', cval=self.background_color, prefilter=False)

        if self.shift_x:
            assert type(self.shift_x) is int, 'Type of shift_x ({}) should be int'.format(type(self.shift_x).__name__)
            if self.shift_x > 0:
                self.array = np.pad(self.array, ((0, 0), (self.shift_x, 0)), mode='constant',
                                   constant_values=(self.background_color))[:, :-self.shift_x]
            else:
                self.array = np.pad(self.array, ((0, 0), (0, -self.shift_x)), mode='constant',
                                   constant_values=(self.background_color))[:, -self.shift_x:]

        if self.shift_y:
            assert type(self.shift_y) is int, 'Type of shift_y ({}) should be int'.format(type(self.shift_y).__name__)
            if self.shift_y < 0:
                self.array = np.pad(self.array, ((-self.shift_y, 0), (0, 0)), mode='constant',
                                   constant_values=(self.background_color))[:self.shift_y, :]
            else:
                self.array = np.pad(self.array, ((0, self.shift_y), (0, 0)), mode='constant',
                                   constant_values=(self.background_color))[self.shift_y:, :]

        if self.invert:
            self.array = np.invert(self.array)


    def save_images(self):
        # self.raw_image.save(os.path.join(self.save_dir, self.raw_image_name))
        self.processed_image.save(os.path.join(self.save_dir, self.processed_image_name))

    def show_images(self):
        self.raw_image.show()
        self.processed_image.show()

    def _add_prefix(self, prefix, name, image_format):
        output_name = '{}_{}'.format(prefix, name) if prefix else name
        return '{}.{}'.format(output_name, image_format)

    def _check_files(self):
        if not os.path.isfile(self.file_path):
            raise ValueError('Provided file "{}" does not exist.'.format(self.file_path))

    def _check_input_type(self, file_path):
        self.possible_extensions = {
            'image': ['tif', 'tiff', 'png', 'bmp', 'gif', 'jpg', 'jpeg'],
            'npy': ['npy'],
        }
        extension = os.path.splitext(file_path)[1][1:].lower()
        for k in self.possible_extensions.keys():
            for e in self.possible_extensions[k]:
                if extension == e:
                    return k
        all_extensions = [x for x_list in self.possible_extensions.values() for x in x_list]
        all_extensions += [x.upper() for x in all_extensions]
        raise ValueError('Incorrect extension: {}. Possible values: {}.'.format(extension, ', '.join(all_extensions)))

# -*- coding: utf-8 -*-
#%%
def opt_from_array(
        array, res_x, res_y, thickness, delta, atten_len,
        arTr, extTr=0, fx=1e+23, fy=1e+23,
        xc=0, yc=0, ne=1, e_start=0, e_fin=0,
        area=None, rotate_angle=None, rotate_reshape=None, cutoff_background_noise=None,
        background_color=None, tile=None, shift_x=None, shift_y=None, invert=None,
        is_save_images=True, prefix='', output_image_format=None,
):
    """Setup Sample element.

    :param file_path: path to the input file (image or .npy).
    :param resolution: resolution of the image [m/pixel].
    :param thickness: thickness of the sample [m].
    :param delta: refractive index decrement.
    :param atten_len: attenuation length [m].
    :param arTr: complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse coordinates.
    :param extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1).
    :param fx: estimated focal length in the horizontal plane [m].
    :param fy: estimated focal length in the vertical plane [m].
    :param xc: horizontal coordinate of center [m].
    :param yc: vertical coordinate of center [m].
    :param ne: number of transmission data points vs photon energy.
    :param e_start: initial photon energy [eV].
    :param e_fin: final photon energy [eV].
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_save_images: a flag to save the initial and processed images.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    :return: transmission (SRWLOptT) type optical element which simulates the Sample.
    """

    # Input parameters to SRWLOptT:
    nx = array.shape[0]
    ny = array.shape[1]
    rx = nx * res_x
    ry = ny * res_y

    opT = SRWLOPT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
                          _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
                          _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    return opT

# ********************** Create transmission element from the data from an image file:
def srwl_opt_setup_transm_from_file(
        inputData, resolution, thickness, delta, atten_len,
        arTr=[], extTr=0, fx=1e+23, fy=1e+23,
        xc=0, yc=0, ne=1, e_start=0, e_fin=0,
        area=None, rotate_angle=None, rotate_reshape=None, cutoff_background_noise=None,
        background_color=None, tile=None, shift_x=None, shift_y=None, invert=None,
        is_save_images=True, prefix='', output_image_format=None, pad=None
):
    """Setup Sample element.

    :param inputData: path to the input file (image or .npy).
    :param resolution: resolution of the image [m/pixel].
    :param thickness: thickness of the sample [m].
    :param delta: refractive index decrement.
    :param atten_len: attenuation length [m].
    :param arTr: complex C-aligned data array (of 2*ne*nx*ny length) storing amplitude transmission and optical path difference as function of transverse coordinates.
    :param extTr: transmission outside the grid/mesh is zero (0), or it is same as on boundary (1).
    :param fx: estimated focal length in the horizontal plane [m].
    :param fy: estimated focal length in the vertical plane [m].
    :param xc: horizontal coordinate of center [m].
    :param yc: vertical coordinate of center [m].
    :param ne: number of transmission data points vs photon energy.
    :param e_start: initial photon energy [eV].
    :param e_fin: final photon energy [eV].
    :param area: the coordinates of the rectangle area listed in the following order: x_start, x_end, y_start, y_end.
    :param rotate_angle: the angle [deg] to rotate the read image counterclockwise. See scipy.ndimage.interpolation.rotate() for details.
    :param rotate_reshape: if reshape is true, the output shape is adapted so that the input array is contained completely in the output.
    :param cutoff_background_noise: the ratio for cutoff the background noise (between 0 and 1).
    :param background_color: the background color code to use instead of the background noise (0=black, 255=white).
    :param tile: the list/tuple (rows, columns) to tile the cut area of the image. See numpy.tile() for details.
    :param shift_x: shift the whole image horizontally. Positive value shifts the image to the right, negative - to the left. See numpy.pad() for details.
    :param shift_y: shift the whole image vertically. Positive value shifts the image to the top, negative - to the bottom. See numpy.pad() for details.
    :param invert: invert the image. See numpy.invert() for details.
    :param is_save_images: a flag to save the initial and processed images.
    :param prefix: the prefix to add to the names of the saved image files.
    :param output_image_format: the format of the output file. If not specified, the input format is used.
    :return: transmission (SRWLOptT) type optical element which simulates the Sample.
    """

    input_parms = {
        "type": "sample",
        "resolution": resolution,
        "thickness": thickness,
        "refractiveIndex": delta,
        "attenuationLength": atten_len,
        "horizontalCenterCoordinate": xc,
        "verticalCenterCoordinate": yc,
        "initialPhotonEnergy": e_start,
        "finalPhotonPnergy": e_fin,
        'area': area,
        'rotateAngle': rotate_angle,
        'rotateReshape': rotate_reshape,
        'cutoffBackgroundNoise': cutoff_background_noise,
        'backgroundColor': background_color,
        'tile': tile,
        'shiftX': shift_x,
        'shiftY': shift_y,
        'invert': invert,
        'outputImageFormat': output_image_format,
        'padding' : pad
    }

    s = np_to_sample(
        inputData=inputData,
        area=area,
        rotate_angle=rotate_angle,
        rotate_reshape=rotate_reshape,
        cutoff_background_noise=cutoff_background_noise,
        background_color=background_color,
        tile=tile,
        shift_x=shift_x,
        shift_y=shift_y,
        invert=invert,
        is_show_images=False,
        is_save_images=is_save_images,
        prefix=prefix,
        output_image_format=output_image_format,
        pad=pad
    )

    # Input parameters to SRWLOptT:
    nx = s.nx
    ny = s.ny

    print('image dimensions: %s x %s' % (str(nx),str(ny)))
    rx = nx * resolution
    ry = ny * resolution

    opT = SRWLOPT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
                          _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
                          _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    data = s.data

    # Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
    offset = 0
    for iy in range(ny):
        for ix in range(nx):
            for ie in range(ne):
                # In images Y=0 corresponds from upper-left corner, in SRW it's lower-left corner:
                pathInBody = thickness * data[ny - iy - 1, ix] / s.limit_value
                opT.arTr[offset] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
                opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
                offset += 2

    opT.input_parms = input_parms

    return opT