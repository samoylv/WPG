# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import os
import math
# Import standart libraries and addnig "../wavefront" directory to python search path
import sys
sys.path.insert(0, os.path.join('..', '..'))
# <codecell>

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

# <codecell>

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

# <codecell>

# NEEDED??
# Write tabulated resulting Wavefront  data to HDF5 file:


def AuxSaveWfrData(arI, wfr, filePath):
# wavefront structure based on reduced glossary:
    wf_struct = {'Version': (0.1, 'f')}
    wf_struct['header'] = {
        'photonEnergy': ((wfr.avgPhotEn), 'f'),
        'nx': ((wfr.mesh.nx), 'i'),
        'ny': ((wfr.mesh.ny), 'i'),
        'xMin': ((wfr.mesh.xStart), 'f'),
        'xMax': ((wfr.mesh.xFin), 'f'),
        'yMin': ((wfr.mesh.yStart), 'f'),
        'yMax': ((wfr.mesh.yFin), 'f'),
        'nSlices': (wf_data.shape[0], 'i'),
        'sliceMin': (slMin, 'f'),
        'sliceMax': (slMax, 'f'),
    }
    wf_struct['data'] = {
        'arrEver': (wf_data, 'f'),
        'arrEhor': (numpy.zeros(shape=wf_data.shape, dtype='float32'), 'f')
    }
    wf_struct['misc'] = {
    }

# <codecell>

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

# <codecell>

# Setup Transmission optical element with 1D heght profile data


def AuxTransmAddSurfHeightProfileScaled(optSlopeErr, heightProfData, dim, ang, scale):
    argHeightProfData = heightProfData[0]
    valHeightProfData = heightProfData[1]
    sinAng = math.sin(ang)
    npData = len(heightProfData[0])
    auxMesh = optSlopeErr.mesh
    xStep = (auxMesh.xFin - auxMesh.xStart)/(auxMesh.nx - 1)
    yStep = (auxMesh.yFin - auxMesh.yStart)/(auxMesh.ny - 1)
    y = auxMesh.yStart
    nyy = auxMesh.ny

    hApprox = 0
    ipStart = 0
    for iy in range(nyy):
        if('y' in dim):
            hApprox = 0
            y1 = argHeightProfData[ipStart] * sinAng
            for i in range(ipStart + 1, npData):
                y2 = argHeightProfData[i] * sinAng
                if((y1 <= y) and (y < y2)):
                    hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / ((argHeightProfData[
                               i] - argHeightProfData[i - 1]) * sinAng)) * (y - y1) + valHeightProfData[i - 1]
                    # print(ipStart, i, iy, y1, y, y2, argHeightProfData[i-1],
                    # argHeightProfData[i], valHeightProfData[i-1], valHeightProfData[i],
                    # hApprox)
                    ipStart = i - 1
                    break
                y1 = y2

        x = auxMesh.xStart
        nxx = auxMesh.nx
        for ix in range(nxx):
            if('x' in dim):
                if(ix == 0):
                    ipStart = 0
                hApprox = 0
                x1 = argHeightProfData[ipStart] * sinAng
                for i in range(ipStart + 1, npData):
                    x2 = argHeightProfData[i] * sinAng
                    if((x1 <= x) and (x < x2)):
                        hApprox = ((valHeightProfData[i] - valHeightProfData[i - 1]) / ((argHeightProfData[
                                   i] - argHeightProfData[i - 1]) * sinAng)) * (x - x1) + valHeightProfData[i - 1]
                        ipStart = i - 1
                        break
                    x1 = x2
            ofst = 2*ix + (2*auxMesh.nx)*iy
            optSlopeErr.arTr[ofst] = 1.  # Amplitude Transmission
            optSlopeErr.arTr[ofst + 1] = 0.  # Optical Path Difference
            if(hApprox != 0):
                optSlopeErr.arTr[ofst + 1] = -2 * sinAng * hApprox * scale  # Optical Path Difference (to check sign!)
            x += xStep
        y += yStep

# <codecell>
