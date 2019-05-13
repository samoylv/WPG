import numpy as np
import math
import scipy

import os
import sys

sys.path.insert(0,os.path.join('..','..'))
sys.path.append(r"C:\Users\daemo\OneDrive\Documents\Github\Chromosome\src")
sys.path.append(r"C:\Users\daemo\OneDrive\Documents\Github\WavefrontPropogation")


from wpg.srwlib import SRWLOptT
from wpg.srwl_uti_smp import *
import wpg.uti_io
def SRWLOptMask(
        file_path, resolution_x, resolution_y, thickness, delta, atten_len,
        arTr=None, extTr=0, fx=1e+23, fy=1e+23,
        xc=0, yc=0, ne=1, e_start=0, e_fin=0,
        area=None, rotate_angle=None, rotate_reshape=None, cutoff_background_noise=None,
        background_color=None, tile=None, shift_x=None, shift_y=None, invert=None,
        is_save_images=False, prefix='', output_image_format=None,
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
    
    
    IN FUTURE FILE_PATH SHOULD BE REPLACE WITH AN EXPLICIT REF TO THE ARRAY IN
    WORKSPACE AS TO AVOID READ/WRITE DELAYS. SHOULD BE COOL FOR NOW THO. 01/05/19
    """

    input_parms = {
        "type": "sample",
        "resolution_x": resolution_x,
        "resolution_y": resolution_y,
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
    }

    s = SRWLUtiSmp(
        file_path=file_path,
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
    )

    # Input parameters to SRWLOptT:
    nx = s.nx
    ny = s.ny
    rx = nx * resolution_x
    ry = ny * resolution_y

    #opT = srwlib.SRWLOptT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
    #                      _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
    #                      _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    data = s.data

    #OC10112018
    specPropAreDef = False
    if(ne > 1):
        if((isinstance(delta, list) or isinstance(delta, array)) and (isinstance(atten_len, list) or isinstance(atten_len, array))):
            lenDelta = len(delta)
            if((lenDelta == len(atten_len)) and (lenDelta == ne)): specPropAreDef = True
            else: raise Exception("Inconsistent spectral refractive index decrement and/or attenuation length data")

    #OC10112018
    useNumPy = False
    try:
        import numpy as np
        useNumPy = True
    except:
        print('NumPy can not be loaded, native Python arrays / lists will be used instead, impacting performance')

    if(useNumPy): #RC161018
        
        thickByLim = thickness/s.limit_value
        nxny = nx*ny
        miHalfThickByLimByAttenLen = None
        miThickDeltaByLim = None

        if(ne <= 1):
            miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len
            miThickDeltaByLim = -thickByLim*delta

            #amplTransm = np.exp(-0.5 * data * thickness / (s.limit_value * atten_len))
            amplTransm = np.exp(miHalfThickByLimByAttenLen*data)

            #optPathDiff =  -1 * data * thickness * delta / s.limit_value
            optPathDiff =  miThickDeltaByLim*data

            arTr = np.empty((2*nxny), dtype=float)
            arTr[0::2] = np.reshape(amplTransm, nxny)
            arTr[1::2] = np.reshape(optPathDiff, nxny)
            #opT.arTr = arTr
            
        else:
            two_ne = 2*ne
            arTr = np.empty((nxny*two_ne), dtype=float)

            if(specPropAreDef == False):
                miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len
                miThickDeltaByLim = -thickByLim*delta
            
            for ie in range(ne):
                if(specPropAreDef):
                    miHalfThickByLimByAttenLen = -0.5*thickByLim/atten_len[ie]
                    miThickDeltaByLim = -thickByLim*delta[ie]

                amplTransm = np.exp(miHalfThickByLimByAttenLen*data)
                optPathDiff = miThickDeltaByLim*data

                two_ie = 2*ie
                arTr[two_ie::two_ne] = np.reshape(amplTransm, nxny) #to check!
                arTr[(two_ie + 1)::two_ne] = np.reshape(optPathDiff, nxny)
    else:
        #Same data alignment as for wavefront: outmost loop vs y, inmost loop vs e
        nTot = 2*ne*nx*ny
        arTr = array('d', [0]*nTot)
    
# =============================================================================
#         offset = 0
#         for iy in range(ny):
#             for ix in range(nx):
#                 #In images Y=0 corresponds from upper-left corner, in SRW it's lower-left corner:
#                 pathInBody = thickness * data[ny - iy - 1, ix] / s.limit_value
#             
#                 #OC10112018
#                 miHalfPathInBody = -0.5*pathInBody
#             
#                 if(specPropAreDef):
#                     for ie in range(ne):
#                         #opT.arTr[offset] = math.exp(-0.5 * pathInBody / atten_len)  # amplitude transmission
#                         #opT.arTr[offset + 1] = -delta * pathInBody  # optical path difference
#                         arTr[offset] = math.exp(miHalfPathInBody / atten_len[ie]) #amplitude transmission
#                         arTr[offset + 1] = -delta[ie] * pathInBody #optical path difference
#                         offset += 2
#                 else:
#                     for ie in range(ne):
#                         arTr[offset] = math.exp(miHalfPathInBody / atten_len) #amplitude transmission
#                         arTr[offset + 1] = -delta * pathInBody  #optical path difference
#                         offset += 2
# =============================================================================

    opT = SRWLOptT(_nx=nx, _ny=ny, _rx=rx, _ry=ry,
                          _arTr=arTr, _extTr=extTr, _Fx=fx, _Fy=fy,
                          _x=xc, _y=yc, _ne=ne, _eStart=e_start, _eFin=e_fin)

    opT.input_parms = input_parms

    return opT