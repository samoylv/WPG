"""
This module contains definitions of custom optical elements.

Described mapping (or aliases) of some of SRW optical elements (SRWLOpt* <-> wpg)

.. module:: wpg.optical_elements
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import errno
from wpg.srwlib import SRWLOptD as Drift
from wpg.srwlib import SRWLOptL as Lens
from wpg.srwlib import srwl, srwl_opt_setup_CRL
import wpg.srwlib
import numpy as np

import sys
if sys.version_info[0] == 3:
    import pickle
else:
    import cPickle as pickle


class WPGOpticalElement(object):

    """Base class for optical elements."""

    def __init__(self):
        pass

class Empty(WPGOpticalElement):

    """Optical element: Empty.
    This is empty propagator used for sampling and zooming wavefront
    """

    def __init__(self):
        super(Empty, self).__init__()

    def __str__(self):
        return ''

    def propagate(self, wfr, propagation_parameters):
        """
        Propagate wavefront through empty propagator,
        used for sampling and resizing wavefront
        """
        beamline = wpg.srwlib.SRWLOptC([], propagation_parameters)
        srwl.PropagElecField(wfr._srwl_wf, beamline)

class Screen(Empty):
    """
    class: Implements the Screen optical element
    """

    def __init__(self, filename=None):
        """ Constructor for the Screen class.

        :param filename: Name of file to store wavefront data.
        :type filename: str
        :raise IOError: File exists.

        """

        # Initialize base class.
        super(Screen, self).__init__()

        # Store filename for output.
        # Handle default.
        if filename is None:
            filename="screen.h5"
        # Check type.
        if not isinstance(filename, (str, unicode)):
            raise TypeError('The parameter "filename" must be str, received %s.' % (type(filename)))
        # Check if parent dir exists.
        filename = os.path.abspath(filename)
        if not os.path.isdir(os.path.dirname(filename)):
            raise IOError('%s is not a directory.' % (os.path.dirname(filename)))
        # Check if file exists. Don't overwrite.
        if os.path.isfile(filename):
            raise IOError('%s already exists. Cowardly refusing to overwrite.')

        self.__filename = filename

    def propagate(self, wfr, propagation_parameters):
        """ Overloaded propagation for this element. """

        super(Screen, self).propagate(wfr, propagation_parameters)
        wfr.store_hdf5(filename)


class Use_PP(object):

    """Short version of propagation parameters. Should be used with `wpg.beamline.Beamline`"""

    def __init__(self, auto_resize_before=None, auto_resize_after=None,
                 releative_precision=None, semi_analytical_treatment=None,
                 fft_resizing=None, zoom=None, zoom_h=None, zoom_v=None,
                 sampling=None, sampling_h=None, sampling_v=None, srw_pp=None
                 ):
        """
        :params srw_pp: propagation parameters in srw style
        """

        super(Use_PP, self).__init__()
        # default values fo propagation parameters
        # Wavefront Propagation Parameters:
        # [0]:  Auto-Resize (1) or not (0) Before propagation
        # [1]:  Auto-Resize (1) or not (0) After propagation
        # [2]:  Relative Precision for propagation with Auto-Resizing (1. is nominal)
        # [3]:  Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation
        # [4]:  Do any Resizing on Fourier side, using FFT, (1) or not (0)
        # [5]:  Horizontal Range modification factor at Resizing (1. means no modification)
        # [6]:  Horizontal Resolution modification factor at Resizing
        # [7]:  Vertical Range modification factor at Resizing
        # [8]:  Vertical Resolution modification factor at Resizing
        # [9]:  Type of wavefront Shift before Resizing (not yet implemented)
        # [10]: New Horizontal wavefront Center position after Shift (not yet implemented)
        # [11]: New Vertical wavefront Center position after Shift (not yet implemented)
        # [12]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Horizontal Coordinate
        # [13]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Vertical Coordinate
        # [14]: Optional: Orientation of the Output Optical Axis vector in the Incident Beam Frame: Longitudinal Coordinate
        # [15]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Horizontal Coordinate
        # [16]: Optional: Orientation of the Horizontal Base vector of the Output Frame in the Incident Beam Frame: Vertical Coordinate

        #         [0][1] [2] [3][4] [5]  [6]  [7]  [8] [9][10][11]
        self.pp = [0, 0, 1.0, 0, 0, 1.0, 1.0, 1.0, 1.0, 0, 0, 0]

        if srw_pp is not None:
            if isinstance(srw_pp, list) and (len(srw_pp) in [12, 17]):
                self.pp = srw_pp

        if auto_resize_before is not None:
            self.auto_resize_before = auto_resize_before

        if auto_resize_after is not None:
            self.auto_resize_after = auto_resize_after

        if releative_precision is not None:
            self.releative_precision = releative_precision

        if semi_analytical_treatment is not None:
            self.semi_analytical_treatment = semi_analytical_treatment

        if fft_resizing is not None:
            self.fft_resizing = fft_resizing

        if zoom is not None:
            self.zoom = zoom

        if zoom_h is not None:
            self.zoom_h = zoom_h

        if zoom_v is not None:
            self.zoom_v = zoom_v

        if sampling is not None:
            self.sampling = sampling

        if sampling_h is not None:
            self.sampling_h = sampling_h

        if sampling_v is not None:
            self.sampling_v = sampling_v

    @property
    def auto_resize_before(self):
        """Auto-Resize (1) or not (0) Before propagation"""
        return self.pp[0]

    @auto_resize_before.setter
    def auto_resize_before(self, value):
        if value in [0, 1]:
            self.pp[0] = value

    @property
    def auto_resize_after(self):
        """Auto-Resize (1) or not (0) After propagation"""
        return self.pp[1]

    @auto_resize_after.setter
    def auto_resize_after(self, value):
        if value in [0, 1]:
            self.pp[1] = value

    @property
    def releative_precision(self):
        """Relative Precision for propagation with Auto-Resizing (1. is nominal)"""
        return self.pp[2]

    @releative_precision.setter
    def releative_precision(self, value):
        self.pp[2] = float(value)

    @property
    def semi_analytical_treatment(self):
        """Allow (1) or not (0) for semi-analytical treatment of quadratic phase terms at propagation"""
        return self.pp[3]

    @semi_analytical_treatment.setter
    def semi_analytical_treatment(self, value):
        if value in [0, 1]:
            self.pp[3] = value

    @property
    def fft_resizing(self):
        """Do any Resizing on Fourier side, using FFT, (1) or not (0)"""
        return self.pp[4]

    @fft_resizing.setter
    def fft_resizing(self, value):
        if value in [0, 1]:
            self.pp[4] = value

    @property
    def zoom_h(self):
        """Horizontal Range modification factor at Resizing (1. means no modification)"""
        return self.pp[5]

    @zoom_h.setter
    def zoom_h(self, value):
        self.pp[5] = float(value)

    @property
    def sampling_h(self):
        """Horizontal Resolution modification factor at Resizing"""
        return self.pp[6]

    @sampling_h.setter
    def sampling_h(self, value):
        self.pp[6] = float(value)

    @property
    def zoom_v(self):
        """Vertical Range modification factor at Resizing"""
        return self.pp[7]

    @zoom_v.setter
    def zoom_v(self, value):
        self.pp[7] = float(value)

    @property
    def sampling_v(self):
        """Vertical Resolution modification factor at Resizing"""
        return self.pp[8]

    @sampling_v.setter
    def sampling_v(self, value):
        self.pp[8] = float(value)

    @property
    def zoom(self):
        """Range modification factor at Resizing (1. means no modification)"""
        if self.zoom_h == self.zoom_v:
            return self.zoom_h
        else:
            return "zoom_h != zoom_v"

    @zoom.setter
    def zoom(self, value):
        z = float(value)
        self.zoom_h = z
        self.zoom_v = z

    @property
    def sampling(self):
        """Resolution modification factor at Resizing"""
        if self.sampling_h == self.sampling_v:
            return self.sampling_h
        else:
            return "sampling_h != sampling_v"

    @sampling.setter
    def sampling(self, value):
        s = float(value)
        self.sampling_h = s
        self.sampling_v = s

    def get_srw_pp(self):
        """
        Return SRW propagation parameters list. Useful for interoperations with SRW tools.

        :return: list of floats (propagation parameters)
        """
        return self.pp

    def __str__(self):
        """
        Print propagation parameters in human readable format.
        """
        return '\n'.join([
            "zoom_h = {}".format(self.zoom_h),
            "zoom_v = {}".format(self.zoom_v),
            "sampling_h = {}".format(self.sampling_h),
            "sampling_v = {}".format(self.sampling_v),
            "semi_analytical_treatment = {}".format(
                self.semi_analytical_treatment),
            "auto_resize_before = {}".format(self.auto_resize_before),
            "auto_resize_after = {}".format(self.auto_resize_after),
            "releative_precision = {}".format(self.releative_precision),
            "fft_resizing = {}".format(self.fft_resizing),
            '\n'
        ])


def Aperture(shape, ap_or_ob, Dx, Dy=1e23, x=0, y=0):
    """
    Defining an aperture/obstacle propagator: A wrapper to a SRWL function SRWLOptA()

    :param shape:    'r' for rectangular, 'c' for circular
    :param ap_or_ob:  'a' for aperture, 'o' for obstacle
    :param Dx, Dy:   transverse dimensions [m]; in case of circular aperture, only Dx is used for diameter
    :param x, y:     transverse coordinates of center [m]
    :return: opAp  - aperture propagator, ``struct SRWLOptA``
    """
    from wpg.srwlib import SRWLOptA

    opAp = SRWLOptA(shape, ap_or_ob, Dx, Dy, x, y)
    return opAp


def Mirror_elliptical(orient, p, q, thetaE, theta0, length):
    """
    Defining a plane elliptical focusing mirror propagator: A wrapper to a SRWL function SRWLOptMirEl()

    :param orient:    mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param p:  distance to one ellipsis center (source), [m]
    :param q:  distance to the other ellipsis center (focus), [m]
    :param thetaE:  design incidence angle in the center of mirror, [rad]
    :param theta0:  "real" incidence angle in the center of mirror, [rad]
    :param length:  mirror length, [m]
    :return: opEFM  - elliptical mirror propagator, ``struct SRWLOptMirEl``
    """
    from wpg.srwlib import SRWLOptMirEl

    if orient == 'x':  # horizontal plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(_p=p, _q=q, _ang_graz=thetaE,
                             _r_sag=1.e+40, _size_tang=length,
                             _nvx=np.cos(theta0), _nvy=0, _nvz=-np.sin(theta0),
                             _tvx=-np.sin(theta0), _tvy=0, _x=0, _y=0, _treat_in_out=1)
    elif orient == 'y':  # vertical plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(_p=p, _q=q, _ang_graz=thetaE,
                             _r_sag=1.e+40, _size_tang=length,
                             _nvx=0, _nvy=np.cos(theta0), _nvz=-np.sin(theta0),
                             _tvx=0, _tvy=-np.sin(theta0), _x=0, _y=0, _treat_in_out=1)
    else:
        raise TypeError('orient should be "x" or "y"')
    return opEFM


def WF_dist(nx, ny, Dx, Dy):
    """
    Create a 'phase screen' propagator for wavefront distortions:   A wrapper to SRWL struct SRWLOptT

    :params nx: number of points in horizontal direction
    :params ny: number of points in vertical   direction
    :params Dx: size in m
    :params Dy: size in
    """
    from wpg.srwlib import SRWLOptT
    return SRWLOptT(nx, ny, Dx, Dy)


def Mirror_plane(orient, theta, length, range_xy, filename, scale=1, delim=' ', xscale = 1, x0 = 0., bPlot=False):
    """
    Defining a plane mirror propagator with taking into account surface height errors

    :param orient:  mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param theta:   incidence angle [rad]
    :param length:  mirror length, [m]
    :range_xy: range in which the incident WF defined [m]
    :filename: full file name with mirror profile of two columns, x and h(x) - heigh errors [m]
    :scale: - height errors scale factor, optical path difference OPD = 2*h*scale*sin(theta)
    :delim: delimiter between data columns
    :param xscale: scaling factor for the mirror profile x-axis (for taking an arbitrary scaled axis, i.e. im mm)
    :x0: shift of mirror longitudinal position [m]
    :return: opIPM  - imperfect plane mirror propagator
    """
    if orient == 'x':  # horizontal plane mirror
        opIPM = WF_dist(1500, 100, length*theta, range_xy)
    elif orient == 'y':  # vertical plane mirror
        opIPM = WF_dist(100, 1500, range_xy, length*theta)
    else:
        raise TypeError('orient should be "x" or "y"')

    calculateOPD(opIPM, mdatafile=filename, ncol=2, delim=delim,
                 Orient=orient, theta=theta, scale=scale,
                 length = length, xscale = xscale, x0=x0, bPlot=bPlot)
    return opIPM


def Mirror_plane_2d(orient, theta, length, range_xy, filename, scale=1, x0=0., y0=0., xscale=1., yscale=1.,bPlot=False):
    """
    Defining a plane mirror propagator with taking into account 2D surface height errors

    :param orient:  mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param theta:   incidence angle [rad]
    :param length:  mirror length, [m]
    :range_xy: range in which the incident WF defined [m]
    :filename: full file name with 2d mirror profile of three columns, x, y, and h(x, y) - heigh errors [m]
    :scale: scale factor, optical path difference OPD = 2*h*scale*sin(theta)
    :x0: shift of mirror longitudinal position [m]
    :y0: shift of mirror transverse position [m]
    :xscale: units of 1st column of filename,  x[m]=x[nits]*xscale  [m]
    :yscale: units of 1st column of filename,  y[m]=y[nits]*yscale  [m]
    :return: opIPM  - imperfect plane mirror propagator
    """
    from scipy import interpolate
    from wpg.srwlib import SRWLOptT
    import os

    sinTheta = np.sin(theta)
    _height_prof_data = np.loadtxt(filename)
    dim = np.shape(_height_prof_data)
    ntotal = dim[0]
    nx = np.size(np.where(_height_prof_data[:, 1] == _height_prof_data[0, 1]))
    ny = int(ntotal/nx)
    print('nx,ny:', nx, ny)

    xax = _height_prof_data[0:nx, 0]*xscale
    xmin = min(xax)
    xmax = max(xax)
    xc = (xmin+xmax)/2
    yax = _height_prof_data[0:ntotal:nx, 1]*yscale
    ymin = min(yax)
    ymax = max(yax)
    yc = (ymin+ymax)/2
    xax = xax - x0 - xc
    xmin = min(xax)
    xmax = max(xax)
    yax = yax - y0 - yc
    ymin = min(yax)
    ymax = max(yax)

    print('length: {:.1f} mm, width: {:.1f} mm'.format(
        (xmax-xmin)*1e3, (ymax-ymin)*1e3))
    if (xmin <= -length/2.) and (xmax >= length/2):
        xmin = -length/2
        xmax = length/2
    else:
        raise ValueError(
            'specified length -{0:.0f}:{0:.0f} mm exceeds \'{1:s}\' mirror limits {2:.0f}:{3:.0f} mm'.format(
                length*1e3/2, os.path.basename(filename), xmin*1e3, xmax*1e3)
        )
    if (ymin <= -range_xy/2) and (ymax >= range_xy/2):
        ymin = -range_xy/2
        ymax = range_xy/2
    else:
        raise ValueError(
            'specified width -{0:.0f}:{0:.0f} mm exceeds \'{1:s}\' mirror limits {2:.0f}:{3:.0f} mm'.format(
                range_xy*1e3/2, os.path.basename(filename), ymin*1e3, ymax*1e3)
        )

    # plt.figure();plt.plot(xax,'bx');plt.plot(yax,'ro');plt.title('xax(blue) yax(red)');
    _height_prof_data_val = np.reshape(_height_prof_data[:, 2], (ny, nx))
    # plt.figure();plt.imshow(_height_prof_data_val);plt.colorbar(orientation='horizontal')
    if orient == 'y':
        opIPM = SRWLOptT(100, 1500, (ymax-ymin), (xmax-xmin)*sinTheta)
    elif orient == 'x':
        opIPM = SRWLOptT(1500, 100, (xmax-xmin)*sinTheta, (ymax-ymin))
    else:
        raise TypeError('orient should be "x" or "y"')
    xnew, ynew = np.mgrid[xmin:xmax:1500j, ymin:ymax:100j]
    f = interpolate.RectBivariateSpline(xax, yax, _height_prof_data_val.T)
    h_new = f(xnew[:, 0], ynew[0, :])
    if bPlot:
        import pylab as plt
        plt.figure();plt.pcolor(xnew, ynew, h_new*scale*1e9);
        plt.axis([xnew.min(), xnew.max(), ynew.min(), ynew.max()])
        plt.colorbar(orientation='horizontal');
        plt.title('surface height errors map, nm');plt.show()
    # print('len:',len(_height_prof_data[2,:]))

    auxMesh = opIPM.mesh
    from array import array
    foo = array(str(u'd'), [])
    # for i in range(150000):
    #     foo.append(1.)
    foo = array(str(u'd'), [1.]*150000)
    opIPM.arTr[::2] = foo  # Amplitude Transmission
    foo = array(str(u'd'), [])
    if orient == 'y':
        for ix in range(1500):
            for iy in range(100):
                foo.append(-2 * sinTheta * h_new[ix, iy] * scale)
    elif orient == 'x':
        for iy in range(100):
            for ix in range(1500):
                foo.append(-2 * sinTheta * h_new[ix, iy] * scale)
    opIPM.arTr[1::2] = foo    # Optical Path Difference (to check sign!)
    return opIPM


def VLS_grating(_mirSub, _m=1, _grDen=100, _grDen1=0, _grDen2=0, _grDen3=0, _grDen4=0, _grAng=0):
    """
    Optical Element: Grating.

    param _mirSub: SRWLOptMir (or derived) type object defining substrate of the grating
    :param _m: output (diffraction) order
    :param _grDen: groove density [lines/mm] (coefficient a0 in the polynomial groove density: a0 + a1*y + a2*y^2 + a3*y^3 + a4*y^4)
    :param _grDen1: groove density polynomial coefficient a1 [lines/mm^2]
    :param _grDen2: groove density polynomial coefficient a2 [lines/mm^3]
    :param _grDen3: groove density polynomial coefficient a3 [lines/mm^4]
    :param _grDen4: groove density polynomial coefficient a4 [lines/mm^5]
    :param _grAng: angle between the grove direction and the saggital direction of the substrate [rad] (by default, groves are made along saggital direction (_grAng=0))
    """

    from .srwlib import SRWLOptG
    return SRWLOptG(_mirSub, _m, _grDen, _grDen1, _grDen2, _grDen3, _grDen4, _grAng)


# def Xtal(_d_sp=5.4309, _psi0r=-0.1446E-4, _psi0i=0.3202E-6, _psi_hr=0.88004E-5, _psi_hi=0.30808E-06,
#          _psi_hbr=0.88004E-5, _psi_hbi=0.30808E-06, _tc=100e-6, _ang_as=0):
#     """
#     Optical Element: Crystal.


#     :param _d_sp: (_d_space) crystal reflecting planes d-spacing (John's dA) [A]
#     :param _psi0r: real part of 0-th Fourier component of crystal polarizability (John's psi0c.real) (units?)
#     :param _psi0i: imaginary part of 0-th Fourier component of crystal polarizability (John's psi0c.imag) (units?)
#     :param _psi_hr: (_psiHr) real part of H-th Fourier component of crystal polarizability (John's psihc.real) (units?)
#     :param _psi_hi: (_psiHi) imaginary part of H-th Fourier component of crystal polarizability (John's psihc.imag) (units?)
#     :param _psi_hbr: (_psiHBr:) real part of -H-th Fourier component of crystal polarizability (John's psimhc.real) (units?)
#     :param _psi_hbi: (_psiHBi:) imaginary part of -H-th Fourier component of crystal polarizability (John's psimhc.imag) (units?)
#     :param _tc: crystal thickness [m] (John's thicum)
#     :param _ang_as: (_Tasym) asymmetry angle [rad] (John's alphdg)
#     :param _nvx: horizontal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
#     :param _nvy: vertical coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
#     :param _nvz: longitudinal coordinate of outward normal to crystal surface (John's angles: thdg, chidg, phidg)
#     :param _tvx: horizontal coordinate of central tangential vector (John's angles: thdg, chidg, phidg)
#     :param _tvy: vertical coordinate of central tangential vector (John's angles: thdg, chidg, phidg)

#     """

#     from .srwlib import SRWLOptC
# return SRWLOptC(_d_sp, _psi0r, _psi0i, _psi_hr, _psi_h, _psi_hbr,
# _psi_hbi, _tc, _ang_as)


# def define_Xtal(xtal='C', h=4, k=0, l=0, tc=100.e-6, idx=0, _dE=0., doPrint=False):
#     """
#     :param xtal: crystal type (now only 'Si' and 'C' are suppoted)
#     :param h,k,l: Miller indices
#     :param tc: crystal thickness, [m]
#     :param idx: the given photon energy line number in crystal parameters tables (ascii files '<xtal>_xih_<hkl>.dat' and '<xtal>_xi0.dat')
#     :param _dE: photon energy offset, eV
#     :param doPrint: if True additional printouts for crystal orientation matrixes
#     :return
#     """
#     # double dSp; /* crystal reflecting planes d-spacing (Angstroems) */
#     # double psi0r, psi0i; /* real and imaginary parts of 0-th Fourier component of crystal polarizability (units?) */
#     # double psiHr, psiHi; /* real and imaginary parts of H-th Fourier component of crystal polarizability (units?) */
#     # double psiHbr, psiHbi; /* real and imaginary parts of -H-th Fourier component of crystal polarizability (units?) */
#     # double tc; /* crystal thickness [m] */
#     # double angAs; /* asymmetry angle [rad] */
#     # double nvx, nvy, nvz; /* horizontal, vertical and longitudinal coordinates of outward normal
#     #                         to crystal surface in the frame of incident beam */
#     # double tvx, tvy; /* horizontal and vertical coordinates of central tangential vector [m] in the frame of incident beam */
#     # C(400)
#     a_Si = 5.4309e-10
#     a_C = 3.5590e-10
#     if (xtal == 'Si'):
#         a = a_Si
#     elif (xtal == 'C'):
#         a = a_C
#     else:
#         print('Unknown Xtal type ', xtal)
#         return
#     aa = np.loadtxt(data_dir+'/%s_xih_%1d%1d%1d.dat' % (xtal, h, k, l))
#     ekev0 = aa[:, 0]
#     xhr = aa[:, 1]
#     xhi = aa[:, 2]
#     Lex = aa[:, 3]
#     d = a/np.sqrt(h**2 + k**2 + l**2)
#     dSp = d*1e10  # 0.88975#e-10
#     idx = 4  # 8.23 keV;
#     print('ekev0[%d] {.4f}'.format(idx, ekev0[idx]))
#     if (ekev0[idx] != wf.params.photonEnergy*1e-3):
#         print('Warning: Central photon energy of the beam {:.4f} keV differs from  ekev0[{:d}] {:.4f}'.format(
#             wf.params.photonEnergy*1e-3, idx, ekev0[idx]))
#     aa = np.loadtxt(data_dir+'/%s_xi0.dat' % (xtal))
#     ekev0 = aa[:, 0]
#     x0r = aa[:, 1]
#     x0i = aa[:, 2]
#     if (ekev0[idx] != wf.params.photonEnergy*1e-3):
#         print('Warning:Central photon energy of the beam {:.4f}keV differs from xi0.dat ekev0[{:d}] {:.4f}'.format(
#             wf.params.photonEnergy*1e-3, idx, ekev0[idx]))
#     psi0r = x0r[idx]
#     psi0i = x0i[idx]
#     psiHr = -xhr[idx]
#     psiHi = xhi[idx]
#     thetaB = np.arcsin(12.39e-10/ekev0[idx]/(2*d))
#     angAs = 0.
#     b = -1
#     DeltaTheta = - psi0r * (1.-1./b)/(2*np.sin(2*thetaB))
#     DeltaE = -DeltaTheta/np.tan(thetaB)*ekev0[idx]
#     print('dTheta_0 {:.1f} urad, dE_0 {:.2f} eV'.format(
#         DeltaTheta*1e6, DeltaE))
#     print('{}({:1d}{:1d}{:1d}) Lex {:.2f} um, thetaB {:.2f} deg'.format(
#         xtal, h, k, l, Lex[idx], thetaB*180/np.pi))
#     Xtal = SRWLOptCryst(_d_sp=dSp, _psi0r=psi0r, _psi0i=psi0i,
#                         _psi_hr=psiHr, _psi_hi=psiHi, _psi_hbr=psiHr, _psi_hbi=psiHi,
#                         _tc=tc, _ang_as=angAs)
#     #Xtal.set_orient(_tc=tc,_ang_as=angAs, _nvx=nvx, _nvy=nvy, _nvz=nvz,_tvx=tvx, _tvy=tvy)
#     # Find appropriate orientation of the Crystal and the corresponding output beam frame:
#     #    """Finds optimal crystal orientation in the input beam frame (i.e. surface normal and tangential vectors) and the orientation of the output beam frame (i.e. coordinates of the longitudinal and horizontal vectors in the input beam frame)
#     #        :param _en: photon energy [eV]
#     #        :param _ang_dif_pl: diffraction plane angle (0 corresponds to the vertical deflection; pi/2 to the horizontal deflection; any value in between is allowed)
#     orientDataXtal = Xtal.find_orient(
#         wf.params.photonEnergy+_dE, _ang_dif_pl=pi/2)
#     orientXtal = orientDataXtal[0]  # Crystal orientation found
#     tXtal = orientXtal[0]
#     nXtal = orientXtal[2]  # Tangential and Normal vectors to crystal surface
#     if doPrint:
#         print('sin(thetaB) {:.4f} \ncos(thetaB) {:.4f}'.format(
#             np.sin(thetaB), np.cos(thetaB)))
#         print('sin(2thetaB) {:.4f} \ncos(2thetaB) {:.4f}'.format(
#             np.sin(2*thetaB), np.cos(2*thetaB)))
#         print('Xtal orientation:')
#         print('tangential vector:\t({:.4f} {:.4f} {:.4f})'.format(
#             tXtal[0], tXtal[1], tXtal[2]))
#         print('normal vector:\t\t({:.4f} {:.4f} {:.4f})'.format(
#             nXtal[0], nXtal[1], nXtal[2]))
#         print(
#             's-vector:\t\t({:.4f} {:.4f} {:.4f})'.format(orientXtal[0][0], orientXtal[0][1], orientXtal[0][2]))
#     # Set orientation of the Crystal:
#     Xtal.set_orient(nXtal[0], nXtal[1], nXtal[2], tXtal[0], tXtal[1])
#     # Orientation of the Outgoing beam frame being found
#     orientOutFrXtal = orientDataXtal[1]
#     # Horizontal, Vertical and Longitudinal base vectors of the Output beam
#     # frame
#     rxXtal = orientOutFrXtal[0]
#     ryXtal = orientOutFrXtal[1]
#     rzXtal = orientOutFrXtal[2]
#     if doPrint:
#         print('Orientation of the Outgoing beam frame:')
#         print('Horizontal base vector:\t\t({:.4f} {:.4f} {:.4f})'.format(
#             orientOutFrXtal[0][0], orientOutFrXtal[0][1], orientOutFrXtal[0][2]))
#         print('Vertical base vector:\t\t({:.4f}  {:.4f} {:.4f})'.format(
#             orientOutFrXtal[1][0], orientOutFrXtal[1][1], orientOutFrXtal[1][2]))
#         print('Longitudinal base vector:\t({:.4f} {:.4f} {:.4f})'.format(
# orientOutFrXtal[2][0], orientOutFrXtal[2][1], orientOutFrXtal[2][2]))

#     return Xtal, rxXtal, ryXtal, rzXtal


# def append_Xtal(bl, xtal='C', h=1, k=1, l=1, tc=100e-6, _dE=0., doPrint=False):
#     """
#     Append to a beamline Xtal propagator

#     :param bl:    Beamline() structure
#     """
#     Xtal, rxXtal, ryXtal, rzXtal = defineXtal(
#         xtal, h=_h, k=_k, l=_l, tc=tc, _dE=_dE, doPrint=False)


def CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n,
        _wall_thick, _xc, _yc, _void_cen_rad=None,
        _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    Setup Transmission type Optical Element which simulates Compound Refractive Lens (CRL).

    :param _foc_plane: plane of focusing: 1- horizontal, 2- vertical, 3- both
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _shape: 1- parabolic, 2- circular (spherical)
    :param _apert_h: horizontal aperture size [m]
    :param _apert_v: vertical aperture size [m]
    :param _r_min: radius (on tip of parabola for parabolic shape) [m]
    :param _n: number of lenses (/"holes")
    :param _wall_thick: min. wall thickness between "holes" [m]
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _void_cen_rad: flat array/list of void center coordinates and radii: [x1, y1, r1, x2, y2, r2,...]
    :param _e_start: initial photon energy
    :param _e_fin: final photon energy
    :return: transmission (SRWLOptT) type optical element which simulates CRL
    """

    return srwl_opt_setup_CRL(_foc_plane, _delta, _atten_len, _shape,
                              _apert_h, _apert_v, _r_min, _n, _wall_thick,
                              _xc, _yc, _void_cen_rad, _e_start, _e_fin, _nx, _ny)


def calculateOPD(wf_dist, mdatafile, ncol, delim, Orient, theta, scale=1., length=1., xscale=1., x0=0., bPlot=False):
    """
    Calculates optical path difference (OPD) from mirror profile and
    fills the struct wf_dist (``struct SRWLOptT``) for wavefront distortions


    :params wf_dist: struct SRWLOptT
    :params mdatafile: an ascii file with mirror profile data
    :params ncol: number of columns in the file
    :params delim: delimiter between numbers in an row, can be space (' '), tab '\t', etc
    :params orient: mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :params theta: incidence angle
    :params scale: scaling factor for the mirror profile
    :param xscale: scaling factor for the mirror profile x-axis (for taking an arbitrary scaled axis, i.e. im mm)
    :param length: mirror length, m, default value 1 m
    :return filled
    """
    from numpy import loadtxt
    # import SRW helpers functions
    from wpg.useful_code.srwutils import AuxTransmAddSurfHeightProfileScaled

    heightProfData = loadtxt(mdatafile).T
    heightProfData[0, :] = heightProfData[0, :] * xscale

    if bPlot:
        import pylab as plt
        plt.figure()
        plt.plot(heightProfData[0, :]*1e3, heightProfData[1, :]*scale*1.e9)
        plt.xlim([(-length+x0)*0.5e3,(length-x0)*0.5e3])
        plt.xlabel('mm'); plt.ylabel('nm')
        plt.title('Height error profile {:s}'.format(mdatafile))
        plt.show()
    AuxTransmAddSurfHeightProfileScaled(
        wf_dist, heightProfData, Orient, theta, scale)
    return wf_dist


def _save_object(obj, file_name):
    """
    Save any python object to file.

    :param: obj : - python objest to be saved
    :param: file_name : - output file, wil be overwrite if exists
    """
    with open(file_name, 'wb') as f:
        pickle.dump(obj, f)


def _load_object(file_name):
    """
    Save any python object to file.

    :param: file_name : - output file, wil be overwrite if exists
    :return: obj : - loaded pthon object
    """
    res = None
    with open(file_name, 'rb') as f:
        res = pickle.load(f)

    return res


def mkdir_p(path):
    """
    Create directory with subfolders (like Linux mkdir -p)

    :param path: Path to be created
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def create_CRL(directory, voids_params,
               _foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n,
               _wall_thick, _xc, _yc, _void_cen_rad=None,
               _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    This function build CLR or load it from file if it was created beforehand.
    Out/input filename builded as sequence of function parameters.

    Adiitinal parameters (*args) passed to srwlib.srwl_opt_setup_CRL function

    :param directory: output directory to save file.
    :param voids_params: void params to build CRL and construct unique file name
    :param _foc_plane: plane of focusing: 1- horizontal, 2- vertical, 3- both
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _shape: 1- parabolic, 2- circular (spherical)
    :param _apert_h: horizontal aperture size [m]
    :param _apert_v: vertical aperture size [m]
    :param _r_min: radius (on tip of parabola for parabolic shape) [m]
    :param _n: number of lenses (/"holes")
    :param _wall_thick: min. wall thickness between "holes" [m]
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _void_cen_rad: flat array/list of void center coordinates and radii: [x1, y1, r1, x2, y2, r2,...]
    :param _e_start: initial photon energy
    :param _e_fin: final photon energy
    :return: SRWL CRL object
    """
    if not isinstance(voids_params, tuple):
        raise TypeError('Voids_params must be tuple')

    file_name = '_'.join([str(a) for a in args[:-1]])
    subdir_name = '_'.join([str(v) for v in voids_params])
    if directory is None:
        full_path = os.path.join(subdir_name, file_name + '.pkl')
    else:
        full_path = os.path.join(directory, subdir_name, file_name + '.pkl')

    if os.path.isfile(full_path):
        print('Found file {}. CLR will be loaded from file'.format(full_path))
        res = _load_object(full_path)
        return res
    else:
        print('CLR file NOT found. CLR will be recalculated and saved in file {}'.format(
            full_path))
        res = CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n,
                  _wall_thick, _xc, _yc, _void_cen_rad,
                  _e_start, _e_fin, _nx, _ny)
        mkdir_p(os.path.dirname(full_path))
        _save_object(res, full_path)
        return res


def create_CRL_from_file(directory, file_name,
                         _foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n,
                         _wall_thick, _xc, _yc, _void_cen_rad=None,
                         _e_start=0, _e_fin=0, _nx=1001, _ny=1001):
    """
    This function build CLR or load it from file.
    Out/input filename builded as sequence of function parameters.
    Adiitinal parameters (*args) passed to srwlib.srwl_opt_setup_CRL function

    :param directory: output directory
    :param fiel_name: CRL file name
    :param _foc_plane: plane of focusing: 1- horizontal, 2- vertical, 3- both
    :param _delta: refractive index decrement (can be one number of array vs photon energy)
    :param _atten_len: attenuation length [m] (can be one number of array vs photon energy)
    :param _shape: 1- parabolic, 2- circular (spherical)
    :param _apert_h: horizontal aperture size [m]
    :param _apert_v: vertical aperture size [m]
    :param _r_min: radius (on tip of parabola for parabolic shape) [m]
    :param _n: number of lenses (/"holes")
    :param _wall_thick: min. wall thickness between "holes" [m]
    :param _xc: horizontal coordinate of center [m]
    :param _yc: vertical coordinate of center [m]
    :param _void_cen_rad: flat array/list of void center coordinates and radii: [x1, y1, r1, x2, y2, r2,...]
    :param _e_start: initial photon energy
    :param _e_fin: final photon energy
    :return: SRWL CRL object
    """

    full_path = os.path.join(directory, file_name + '.pkl')

    if os.path.isfile(full_path):
        print('Found file {}. CLR will be loaded from file'.format(full_path))
        res = _load_object(full_path)
        return res
    else:
        print('CLR file NOT found. CLR will be recalculated and saved in file {}'.format(
            full_path))
        res = CRL(_foc_plane, _delta, _atten_len, _shape, _apert_h, _apert_v, _r_min, _n,
                  _wall_thick, _xc, _yc, _void_cen_rad,
                  _e_start, _e_fin, _nx, _ny)
        mkdir_p(os.path.dirname(full_path))
        _save_object(res, full_path)
        return res
