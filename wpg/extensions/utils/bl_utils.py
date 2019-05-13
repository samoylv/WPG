import math
import numpy as np


def mirror_rot(beam_diameter, mir_width):
    """
    Determine the appropriate angle of rotation of plane mirrors
    
    :param beam_diameter: width of the incident beam
    :param mir_width: width of the plane mirror
    """
    theta1 = math.asin(beam_diameter/mir_width)
    theta2 = np.pi/2 - math.acos(beam_diameter/mir_width)
    
    return theta1, theta2


def bragg_angle(wavelength, d):
    """
    Calculates bragg angle of a crystal
    
    :param wavelength: wavelength of incident light 
    :param d: lattice spacing of crystal
    """
    
    theta = math.asin((wavelength)/(2*d))
    return theta

def cyl_focus(f, alpha0):
    """
    
    TO BE DEVELOPED 
    
    Calculates the sag and tang radius of a toroidal mirror
    
    :param f: desired focal length (m)
    :param alpha0: grazing angle (deg)
    """
    
    alpha0 = math.radians(alpha0)
    alpha = np.pi/2-alpha0
    
    r = 2*f*math.cos(alpha)
    R = (2*f)/(math.cos(alpha))
    
    return r, R

def toroid_focus(f, alpha0):
    """
    Calculates the sag and tang radius of a toroidal mirror
    
    :param f: desired focal length (m)
    :param alpha0: grazing angle (deg)
    """
    
    alpha0 = math.radians(alpha0)
    alpha = np.pi/2-alpha0
    
    r = 2*f*math.cos(alpha)
    R = (2*f)/(math.cos(alpha))
    
    return r, R

def ZP_focu_gvr(D, rn, wavelength, order = 1):
    
    """
    Calculates the focal point of a zone plate
    
    :param D: outer diameter [m]
    :param oz: width of the outer zone [m]
    :param wavelength: wavelength of radiation [m]
    :param order: diffraction order
    """
    
    f = (D*(rn)) / ( order*wavelength )
    
    return f

def ZP_focus(r_n, n, wav = 6.7e-09):
    
    """
    Calculates the focal point of a zone plate
    
    :param r_n: radius of outermost zone [m]
    :param n: number of zones
    :param wav: wavelength of light [m]
    """
    a = 1/(wav*n)
    b = r_n**2
    c = (n*wav)/4
    
    s = a*b-c
    
    return s

def PGM_delta_y(width = 1e-03, thetai = 1.5, wav = 6.7e-09, deltawav = 1e-09, d = 8.334e-07):
    
    """
    Calculates the ideal plane mirror grating seperation to seperate
    polychromatic light to a bandwidth deltawav.
    
    Note that this formula assumes that the size of the incident wavefield is
    negligible in comparison to the size of the PGM and strikes at the centre
    
    :param width: width of the PGM [m]
    

    :param thetai: angle of incidence [degrees]
    :param wav: central wavelength [m]
    :param deltawav: desired spectral bandwidth [m]
    :param d: grating separation [m]        
    """

    
    a = (2*width)/d
    b = d*math.sin(thetai)-(wav+deltawav)
    delta_y = a*b

    return delta_y


def wfr_center(wfr):
    [xc, yc] = calculate_peak_pos(wfr)
    wfr.params.Mesh.xMin -= xc
    wfr.params.Mesh.xMax -= xc
    wfr.params.Mesh.yMin -= yc
    wfr.params.Mesh.yMax -= yc
    
def get_mesh(mwf):
    wf_mesh = mwf.params.Mesh
    nx = wf_mesh.nx
    ny = wf_mesh.ny
    [xmin, xmax, ymin, ymax] = [wf_mesh.xMin,
                                wf_mesh.xMax, wf_mesh.yMin, wf_mesh.yMax]
    return [nx, ny, xmin, xmax, ymin, ymax]

    
def calculate_peak_pos(mwf):
    # irradiance
    irr = mwf.get_intensity(slice_number=0, polarization='horizontal')
    irr_max = np.max(irr)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    x_axis = np.linspace(xmin, xmax, nx)
    y_axis = np.linspace(ymin, ymax, ny)
    nc = np.where(irr == irr_max)
    irr_x = irr[ny // 2, :]
    irr_y = irr[:, nx // 2]
    x0 = np.max(x_axis[np.where(irr_x == np.max(irr_x))])
    y0 = np.max(y_axis[np.where(irr_y == np.max(irr_y))])
    return [x0, y0]



