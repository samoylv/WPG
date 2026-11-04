# -*- coding: utf-8 -*-
__author__ = 'O. Chubar, A. Buzmakov'

import os
import math

#**********************Auxiliary Functions
# Read data comumns from ASCII file:


def AuxReadInDataColumns(filePath, nCol, strSep):
    f = open(filePath, 'r')
    resCols = []
    for iCol in range(nCol):
        resCols.append([])

    curLine = f.readline()
    while len(curLine) > 0:
        curLineParts = curLine.split(strSep)
        for iCol in range(nCol):
            if(iCol < len(curLineParts)):
                resCols[iCol].append(float(curLineParts[iCol]))
        curLine = f.readline()
    f.close()
    return resCols  # attn: returns lists, not arrays!


# Write tabulated resulting Intensity data to ASCII file:


def AuxSaveIntData(arI, wfr, filePath):
    f = open(filePath, 'w')
    f.write('#C-aligned Intensity (inner loop is vs photon energy, outer loop vs vertical position)\n')
    f.write('#' + repr(wfr.eStart) + ' #Initial Photon Energy [eV]\n')
    f.write('#' + repr(wfr.eFin) + ' #Final Photon Energy [eV]\n')
    f.write('#' + repr(wfr.ne) + ' #Number of points vs Photon Energy\n')
    f.write('#' + repr(wfr.xStart) + ' #Initial Horizontal Position [m]\n')
    f.write('#' + repr(wfr.xFin) + ' #Final Horizontal Position [m]\n')
    f.write('#' + repr(wfr.nx) + ' #Number of points vs Horizontal Position\n')
    f.write('#' + repr(wfr.yStart) + ' #Initial Vertical Position [m]\n')
    f.write('#' + repr(wfr.yFin) + ' #Final Vertical Position [m]\n')
    f.write('#' + repr(wfr.ny) + ' #Number of points vs Vertical Position\n')
    for i in range(wfr.ne * wfr.nx * wfr.ny):  # write all data into one column using "C-alignment" as a "flat" 1D array
        f.write(' ' + repr(arI[i]) + '\n')
    f.close()


# NEEDED??
# Write tabulated resulting Wavefront  data to HDF5 file:

#TODO: delete ?
# def AuxSaveWfrData(arI, wfr, filePath):
# # wavefront structure based on reduced glossary:
#     wf_struct = {'Version': (0.1, 'f')}
#     wf_struct['header'] = {
#         'photonEnergy': ((wfr.avgPhotEn), 'f'),
#         'nx': ((wfr.mesh.nx), 'i'),
#         'ny': ((wfr.mesh.ny), 'i'),
#         'xMin': ((wfr.mesh.xStart), 'f'),
#         'xMax': ((wfr.mesh.xFin), 'f'),
#         'yMin': ((wfr.mesh.yStart), 'f'),
#         'yMax': ((wfr.mesh.yFin), 'f'),
#         'nSlices': (wf_data.shape[0], 'i'),
#         'sliceMin': (slMin, 'f'),
#         'sliceMax': (slMax, 'f'),
#     }
#     wf_struct['data'] = {
#         'arrEver': (wf_data, 'f'),
#         'arrEhor': (numpy.zeros(shape=wf_data.shape, dtype='float32'), 'f')
#     }
#     wf_struct['misc'] = {
#     }


# Write Optical Transmission characteristic data to ASCII file:


def AuxSaveOpTransmData(optTr, t, filePath):
    f = open(filePath, 'w')
    f.write(
        '#C-aligned optical Transmission characteristic (inner loop is vs horizontal position, outer loop vs vertical position)\n')
    f.write('#' + repr(1) + ' #Reserved for Initial Photon Energy [eV]\n')
    f.write('#' + repr(1) + ' #Reserved for Final Photon Energy [eV]\n')
    f.write('#' + repr(1) + ' #Reserved for Number of points vs Photon Energy\n')
    auxMesh = optTr.mesh
    f.write('#' + repr(auxMesh.xStart) + ' #Initial Horizontal Position [m]\n')
    f.write('#' + repr(auxMesh.xFin) + ' #Final Horizontal Position [m]\n')
    f.write('#' + repr(auxMesh.nx) + ' #Number of points vs Horizontal Position\n')
    f.write('#' + repr(auxMesh.yStart) + ' #Initial Vertical Position [m]\n')
    f.write('#' + repr(auxMesh.yFin) + ' #Final Vertical Position [m]\n')
    f.write('#' + repr(auxMesh.ny) + ' #Number of points vs Vertical Position\n')
    ntot = auxMesh.nx * auxMesh.ny

    for i in range(ntot):
        tr = 0
        if((t == 1) or (t == 2)):  # amplitude or intensity transmission
            tr = optTr.arTr[i * 2]
            if(t == 2):  # intensity transmission
                tr *= tr
        else:  # optical path difference
            tr = optTr.arTr[i * 2 + 1]
        f.write(' ' + repr(tr) + '\n')
    f.close()


# Setup Transmission optical element with 1D heght profile data


def AuxTransmAddSurfHeightProfileScaled(optSlopeErr, heightProfData, dim, ang, scale):
    """
    Fill a transmission element (SRWLOptT) with the optical path difference of a 1D mirror height
    profile: OPD = -2 sin(ang) h scale, amplitude 1, linear interpolation, 0 outside the profile.

    Vectorized version of the former loops (the same result, tests/test_surface.py); the profile is
    sorted by position (the loops gave zeros for a profile with decreasing positions).

    :param optSlopeErr: SRWLOptT, its mesh is the beam cross-section at the mirror
    :param heightProfData: [positions along the mirror [m], heights [m]]
    :param dim: plane of deflection, 'x' or 'y'
    :param ang: grazing angle [rad]
    :param scale: factor for the heights
    """
    import numpy as np
    from wpg.surface import _interp_inside, _set_transmission, _sorted_profile
    arg, val = _sorted_profile(heightProfData[0], heightProfData[1])
    m = optSlopeErr.mesh
    sin_ang = math.sin(ang)
    x = np.linspace(m.xStart, m.xFin, m.nx)
    y = np.linspace(m.yStart, m.yFin, m.ny)
    if 'x' in dim:
        h = _interp_inside(x, arg * sin_ang, val)[None, :] * np.ones((m.ny, 1))
    elif 'y' in dim:
        h = _interp_inside(y, arg * sin_ang, val)[:, None] * np.ones((1, m.nx))
    else:
        h = np.zeros((m.ny, m.nx))
    _set_transmission(optSlopeErr, -2 * sin_ang * h * scale)

