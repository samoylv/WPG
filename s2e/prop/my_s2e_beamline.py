import numpy as np
import os
from wpg import Wavefront
import wpg.optical_elements
from wpg.optical_elements import Use_PP

#mirror_data_dir = '/data/S2E/modules/prop/data_common'
mirror_data_dir = 'data_common'

def defineOPD(opTrErMirr, mdatafile, ncol, delim, orient, theta, scale=1., stretching=1.):
    """
    Define optical path difference (OPD) from mirror profile
    
    :param opTrErMirr: the struct with wave front distortions from mirror susrface errors
    :param mdatafile: an ascii file with mirror profile data
    :param ncol: number of columns in the mirror profile file
    :param delim: delimiter between numbers in an row, can be space (' '), tab '\t', etc
    :param orient: mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param theta: mirror incidence angle
    :param scale: scaling factor for the mirror profile height errors
    :param stretching: scaling factor for the mirror profile x-axis (a hack, should be removed ASAP) 
    """

    heightProfData = np.loadtxt(mdatafile).T
    heightProfData[0,:] = heightProfData[0,:] * stretching
    wpg.useful_code.srwutils.AuxTransmAddSurfHeightProfileScaled(opTrErMirr, heightProfData, orient, theta, scale)
    # if isIpynb:
    #     pylab.figure(); pylab.plot(heightProfData[0],heightProfData[ncol-1]*1e9)
    #     pylab.title('profile from %s' %mdatafile);pylab.xlabel('x (m)');pylab.ylabel('h (nm)')

def get_beamline():
    """
    This function will called in propagation script.
    
    :return: Beamline.
    """

    distance0 = 246.5
    distance1 = 683.5
    distance = distance0 + distance1
    f_hfm    = 3.0       # nominal focal length for HFM KB
    f_vfm    = 1.9       # nominal focal length for VFM KB
    distance_hfm_vfm = f_hfm - f_vfm
    distance_foc =  1. /(1./f_vfm + 1. / (distance + distance_hfm_vfm))
    theta_om = 3.5e-3 # offset mirrors incidence angle 
    theta_kb = 3.5e-3 # KB mirrors incidence angle 
    
    drift0 = wpg.optical_elements.Drift(distance0)
    drift1 = wpg.optical_elements.Drift(distance1)
    drift_in_kb = wpg.optical_elements.Drift(distance_hfm_vfm)
    drift_to_foc = wpg.optical_elements.Drift(distance_foc)
    
    om_mirror_length = 0.8; om_clear_ap = om_mirror_length*theta_om
    kb_mirror_length = 0.9; kb_clear_ap = kb_mirror_length*theta_kb
    ap0   = wpg.optical_elements.Aperture('r','a', 120.e-6, 120.e-6)
    ap1   = wpg.optical_elements.Aperture('r','a', om_clear_ap, 2*om_clear_ap)
    ap_kb = wpg.optical_elements.Aperture('r','a', kb_clear_ap, kb_clear_ap)
    hfm    = wpg.optical_elements.Mirror_elliptical(
                    _p=distance, _q=(distance_hfm_vfm+distance_foc), _ang_graz=theta_kb,
                    _r_sag=1.e+40, _size_tang=0.9, _nvx=np.cos(theta_kb),  _nvy=0,
                    _nvz=-np.sin(theta_kb), _tvx=-np.sin(theta_kb), _tvy=0, _x=0, _y=0,
                    _treat_in_out=1)     
    vfm    = wpg.optical_elements.Mirror_elliptical(
                    _p=(distance+distance_hfm_vfm), _q=distance_foc, _ang_graz=theta_kb,
                    _r_sag=1.e+40, _size_tang=0.9, _nvx=0, _nvy=np.cos(theta_kb),
                    _nvz=-np.sin(theta_kb), _tvx=0, _tvy=-np.sin(theta_kb), _x=0, _y=0,
                    _treat_in_out=1) 
    wf_dist_om = wpg.optical_elements.WF_dist(1500, 100, om_clear_ap, 2*om_clear_ap)
    defineOPD(wf_dist_om, os.path.join(mirror_data_dir,'mirror2.dat'), 2, '\t', 'x',
              theta_kb, scale=2)
    # if isIpynb:
    #     meshT = wf_dist_om.mesh
    #     opdTmp=np.array(wf_dist_om.arTr)[1::2].reshape(meshT.ny,meshT.nx)
    #     figure(); pylab.imshow(opdTmp,extent=[meshT.xStart,meshT.xFin,meshT.yStart,meshT.yFin])
    #     pylab.title('OPD [m]');pylab.xlabel('x (m)'); pylab.ylabel('y (m)')
        
    wf_dist_hfm = wpg.optical_elements.WF_dist(1500, 100, kb_clear_ap, kb_clear_ap)
    defineOPD(wf_dist_hfm, os.path.join(mirror_data_dir,'mirror1.dat'), 2, '\t', 'x',
              theta_kb, scale=2, stretching=kb_mirror_length/0.8)
    # if isIpynb:
    #     meshT = wf_dist_hfm.mesh
    #     opdTmp=np.array(wf_dist_hfm.arTr)[1::2].reshape(meshT.ny,meshT.nx)
    #     figure(); pylab.imshow(opdTmp,extent=[meshT.xStart,meshT.xFin,meshT.yStart,meshT.yFin])
    #     pylab.title('OPD [m]');pylab.xlabel('x (m)'); pylab.ylabel('y (m)')  
    
    wf_dist_vfm = wpg.optical_elements.WF_dist(1100, 1500, kb_clear_ap, kb_clear_ap)
    defineOPD(wf_dist_vfm, os.path.join(mirror_data_dir,'mirror2.dat'), 2, ' ', 'y',
              theta_kb, scale=2, stretching=kb_mirror_length/0.8)
    
    # if isIpynb:
    #     meshT = wf_dist_vfm.mesh
    #     opdTmp=np.array(wf_dist_vfm.arTr)[1::2].reshape(meshT.ny,meshT.nx)
    #     figure(); pylab.imshow(opdTmp,extent=[meshT.xStart,meshT.xFin,meshT.yStart,meshT.yFin])
    #     pylab.title('OPD [m]');pylab.xlabel('x (m)'); pylab.ylabel('y (m)')
        
    bl0 = wpg.Beamline()
    bl0.append(ap0,   Use_PP(semi_analytical_treatment=0, zoom=14.4, sampling=1/1.6))
    bl0.append(drift0,Use_PP(semi_analytical_treatment=0))
    bl0.append(ap1,    Use_PP(zoom=0.8))   #bl0.append(ap1,    Use_PP(zoom=1.6, sampling=1/1.5))
    bl0.append(wf_dist_om, Use_PP())
    bl0.append(drift1, Use_PP(semi_analytical_treatment=1))
    bl0.append(ap_kb,  Use_PP(zoom = 6.4, sampling = 1/16.))#bl0.append(ap_kb,    Use_PP(zoom=5.4, sampling=1/6.4))
    bl0.append(hfm, Use_PP())
    bl0.append(wf_dist_hfm, Use_PP())
    bl0.append(drift_in_kb, Use_PP(semi_analytical_treatment=1))
    bl0.append(vfm, Use_PP())
    bl0.append(wf_dist_vfm, Use_PP())
    bl0.append(drift_to_foc, Use_PP(semi_analytical_treatment=1))

    return bl0
