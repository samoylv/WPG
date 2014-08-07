# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov'

import warnings
import srwlib
import srwlpy


def build_gauss_wavefront(nx, ny, nz, ekev, xMin, xMax, yMin, yMax, tau, sigX, sigY, d2waist, pulseEn=None,_mx=None,_my=None):
    """
    Build 3D Gaussian beam.

    :param nx: Number of point along x-axis
    :param ny: Number of point along y-axis
    :param nz: Number of point along z-axis (slices)
    :param ekev: Energy in kEv
    :param xMin: Initial Horizontal Position [m]
    :param xMax: Final Horizontal Position [m]
    :param yMin: Initial Vertical Position [m]
    :param yMax: Final Vertical Position [m]
    :param tau: Pulse duration [s]
    :param sigX: Horiz. RMS size at Waist [m]
    :param sigY:  Vert. RMS size at Waist [m]
    :param d2waist: Distance to Gaussian waist
    :param pulseEn: Energy per Pulse [J]
    :param _mx: transverse Gauss-Hermite mode order in horizontal direction
    :param _my: transverse Gauss-Hermite mode order in vertical direction
    :return: wpg.Wavefront structure
    """
    # TODO: fix comment

    GsnBm = srwlib.SRWLGsnBm()  # Gaussian Beam structure (just parameters)
    GsnBm.x = 0  # Transverse Coordinates of Gaussian Beam Center at Waist [m]
    GsnBm.y = 0
    GsnBm.z = 0  # Longitudinal Coordinate of Waist [m]
    GsnBm.xp = 0  # Average Angles of Gaussian Beam at Waist [rad]
    GsnBm.yp = 0

    GsnBm.avgPhotEn = ekev * 1.e3  # 15000. #Photon Energy [eV]
    if pulseEn is not None:
        GsnBm.pulseEn = pulseEn 
    else:
        GsnBm.pulseEn = 0.001 # was 1 mJ in the Tutorial exampes as well
    GsnBm.repRate = 1  # Rep. Rate [Hz] - to be corrected
    GsnBm.polar = 2  # 1- linear hoirizontal
    # Far field angular divergence: 14.1e-6 ./ (ekev) .^0.75
    # 0.17712e-09/(4*Pi)/(14.1e-06/((7)^0.75)) for 7 keV, 3.55561e-06 =
    # 0.0826561e-09/(4*Pi)/(14.1e-06/((15)^0.75)) for 15 keV #Horiz. RMS size
    # at Waist [m]
    GsnBm.sigX = sigX
    GsnBm.sigY = sigY  # Vert. RMS size at Waist [m]
    # Coherence time (~ Gaussian pulse duration)           0.12 fs @ 15 keV
    # and 0.17 fs @ 7 keV
    # 0.12e-15 #0.17e-15 for 15 keV #Pulse duration [s] #To check: Is it 0.12
    # fs or 12 fs ?
    GsnBm.sigT = tau
    if _mx is not None:
        GsnBm.mx = _mx 
    else:
        GsnBm.mx = 0  # Transverse Gauss-Hermite Mode Orders
    if _mx is not None:
        GsnBm.my = _my 
    else:
        GsnBm.my = 0

    wfr = srwlib.SRWLWfr()  # Initial Electric Field Wavefront
    wfr.allocate(nz, nx, ny)
     # Numbers of points vs Photon Energy (1), Horizontal and
     # Vertical Positions (dummy)
    wfr.presFT = 1  # Defining Initial Wavefront in Time Domain
    #wfr.presFT = 0 #Defining Initial Wavefront in Frequency Domain

    wfr.avgPhotEn = GsnBm.avgPhotEn
    wfr.mesh.eStart = -3 * GsnBm.sigT  # Initial Time [s]
    wfr.mesh.eFin = 3 * GsnBm.sigT  # Final Time [s]

    # Longitudinal Position [m] at which Electric Field has to be calculated,
    # i.e. the position of the first optical element
    wfr.mesh.zStart = d2waist
    wfr.mesh.xStart = xMin  # Initial Horizontal Position [m]
    wfr.mesh.xFin = xMax  # Final Horizontal Position [m]
    wfr.mesh.yStart = yMin  # Initial Vertical Position [m]
    wfr.mesh.yFin = yMax  # Final Vertical Position [m]

    wfr.mesh.ne = nz

    # Some information about the source in the Wavefront structure
    wfr.partBeam.partStatMom1.x = GsnBm.x
    wfr.partBeam.partStatMom1.y = GsnBm.y
    wfr.partBeam.partStatMom1.z = GsnBm.z
    wfr.partBeam.partStatMom1.xp = GsnBm.xp
    wfr.partBeam.partStatMom1.yp = GsnBm.yp

    sampFactNxNyForProp = - \
        1  # 5 #sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [sampFactNxNyForProp]
    #**********************Calculating Initial Wavefront
    srwlpy.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)

    return wfr


def build_gauss_wavefront_xy(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, d2waist,
                             xoff=0., yoff=0., tiltX=0., tiltY=0., 
                             pulseEn=None, pulseTau=None,repRate=None,_mx=None,_my=None):
    """
    Build 2D Gaussian beam.
    
    :param nx: Number of point along x-axis
    :param ny: Number of point along y-axis
    :param nz: Number of point along z-axis (slices)
    :param ekev: Energy in kEv
    :param xMin: Initial Horizontal Position [m]
    :param xMax: Final Horizontal Position [m]
    :param yMin: Initial Vertical Position [m]
    :param yMax: Final Vertical Position [m]
    :param sigX: Horiz. RMS size at Waist [m]
    :param sigY:  Vert. RMS size at Waist [m]
    :param d2waist: distance to Gaussian waist
    :param xoff:    Horizonal Coordinate of Gaussian Beam Center at Waist [m]
    :param yoff:    Vertical  Coordinate of Gaussian Beam Center at Waist [m]
    :param tiltX:   Average Angle of Gaussian Beam at Waist in Horizontal plane [rad] 
    :param tiltY:   Average Angle of Gaussian Beam at Waist in Vertical plane [rad]
    :param pulseEn:  Energy per Pulse [J]
    :param pulseTau: Coherence time [s] to get proper BW
    :param _mx: transverse Gauss-Hermite mode order in horizontal direction
    :param _my: transverse Gauss-Hermite mode order in vertical direction
    :return: wpg.Wavefront structure

    """
    GsnBm = srwlib.SRWLGsnBm()  # Gaussian Beam structure (just parameters)
    # Transverse Coordinates of Gaussian Beam Center at Waist [m]
    GsnBm.x = xoff
    GsnBm.y = yoff
    GsnBm.z = 0  # Longitudinal Coordinate of Waist [m]
    GsnBm.xp = tiltX  # Average Angles of Gaussian Beam at Waist [rad]
    GsnBm.yp = tiltY
    GsnBm.avgPhotEn = ekev * 1e3  # 5000 #Photon Energy [eV]
    if pulseEn is not None:
        GsnBm.pulseEn = pulseEn 
    else:
        GsnBm.pulseEn = 0.001  # Energy per Pulse [J] - to be corrected
    if repRate is not None:
        GsnBm.repRate = repRate 
    else:
        GsnBm.repRate = 1  # Rep. Rate [Hz] - to be corrected
    GsnBm.polar = 2  # 1- linear hoirizontal
    GsnBm.sigX = sigX  # Horiz. RMS size at Waist [m]
    GsnBm.sigY = sigY  # Vert. RMS size at Waist [m]
    if pulseTau is not None:
        GsnBm.sigT = pulseTau 
    else:
        GsnBm.sigT = 0.2e-15  # should be about coherence time to get proper BW
    if _mx is not None:
        GsnBm.mx = _mx 
    else:
        GsnBm.mx = 0  # Transverse Gauss-Hermite Mode Orders
    if _mx is not None:
        GsnBm.my = _my 
    else:
        GsnBm.my = 0

    wfr = srwlib.SRWLWfr()  # Initial Electric Field Wavefront
    wfr.allocate(1, nx, ny)
     # Numbers of points vs Photon Energy (1), Horizontal and
     # Vertical Positions (dummy)
    wfr.mesh.eStart = GsnBm.avgPhotEn  # Initial Photon Energy [eV]
    wfr.mesh.eFin = GsnBm.avgPhotEn  # Final Photon Energy [eV]
    wfr.avgPhotEn = (wfr.mesh.eStart + wfr.mesh.eFin) / 2
    wfr.mesh.zStart = d2waist
    wfr.mesh.xStart = xMin  # Initial Horizontal Position [m]
    wfr.mesh.xFin = xMax  # Final Horizontal Position [m]
    wfr.mesh.yStart = yMin  # Initial Vertical Position [m]
    wfr.mesh.yFin = yMax  # Final Vertical Position [m]

    #wfr.presFT = 1  # Defining Initial Wavefront in Time Domain
    wfr.presFT = 0  # Defining Initial Wavefront in Freq Domain

    # Some information about the source in the Wavefront structure
    wfr.partBeam.partStatMom1.x = GsnBm.x
    wfr.partBeam.partStatMom1.y = GsnBm.y
    wfr.partBeam.partStatMom1.z = GsnBm.z
    wfr.partBeam.partStatMom1.xp = GsnBm.xp
    wfr.partBeam.partStatMom1.yp = GsnBm.yp

    if (pulseEn is None) and (pulseTau is None):
        wfr.unitElFld = 0 # set to 'arbitrary'

    sampFactNxNyForProp = -1  # sampling factor for adjusting nx, ny (effective if > 0)
    arPrecPar = [sampFactNxNyForProp]
    srwlpy.CalcElecFieldGaussian(wfr, GsnBm, arPrecPar)
    return wfr


def build_gauss_wavefront_xy_(xoff, yoff, tiltX, tiltY, nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, d2waist):
    """
    This function depricated and will removed in next releases, use build_gauss_wavefront_xy instead. 
    """
    warnings.warn('This function depricated and will removed in next releases, use build_gauss_wavefront_xy instead.', DeprecationWarning)
    return build_gauss_wavefront_xy(nx, ny, ekev, xMin, xMax, yMin, yMax, sigX, sigY, d2waist, xoff, yoff, tiltX, tiltY)