# Module holding WPG beamlines at European XFEL.

from opd import defineOPD
from wpg import optical_elements, Beamline
from wpg.optical_elements import Use_PP
import numpy
import os

from propagate_s2e import MIRROR_DATA_DIR as mirror_data_dir


def get_beamline():
    """ Setup and return the WPG.Beamline object representing the SPB/SFX nanofocus beamline (KB mirrors).

    :return: beamline
    :rtype: wpg.Beamline
    """

    # Distances
    distance0 = 246.5
    distance1 = 683.5
    distance = distance0 + distance1

    # Focal lengths.
    f_hfm    = 3.0          # nominal focal length for HFM KB
    f_vfm    = 1.9          # nominal focal length for VFM KB
    distance_hfm_vfm = f_hfm - f_vfm
    distance_foc =  1. /(1./f_vfm + 1. / (distance + distance_hfm_vfm))

    # Mirror incidence angles
    theta_om = 3.5e-3       # offset mirrors incidence angle
    theta_kb = 3.5e-3       # KB mirrors incidence angle

    # Mirror lengths
    om_mirror_length = 0.8; om_clear_ap = om_mirror_length*theta_om
    kb_mirror_length = 0.9; kb_clear_ap = kb_mirror_length*theta_kb

    # Drifts.
    drift0 = optical_elements.Drift(distance0)
    drift1 = optical_elements.Drift(distance1)
    drift_in_kb = optical_elements.Drift(distance_hfm_vfm)
    drift_to_foc = optical_elements.Drift(distance_foc)

    # Mirror apertures.
    ap0   = optical_elements.Aperture('r','a', 120.e-6, 120.e-6)
    ap1   = optical_elements.Aperture('r','a', om_clear_ap, 2*om_clear_ap)
    ap_kb = optical_elements.Aperture('r','a', kb_clear_ap, kb_clear_ap)

    # Mirror definitions.
    hfm = optical_elements.Mirror_elliptical(
                    orient='x',
                    p=distance,
                    q=(distance_hfm_vfm+distance_foc),
                    thetaE=theta_kb,
                    theta0=theta_kb,
                    length=0.9
                    )
    vfm = optical_elements.Mirror_elliptical(
                    orient='y',
                    p=(distance+distance_hfm_vfm),
                    q=distance_foc,
                    thetaE=theta_kb,
                    theta0=theta_kb,
                    length=0.9
                    )


    # Wavefront distortions due to mirror profile.
    wf_dist_om = optical_elements.WF_dist(1500, 100, om_clear_ap, 2*om_clear_ap)
    defineOPD(wf_dist_om, os.path.join(mirror_data_dir,'mirror2.dat'), 2, '\t', 'x',
              theta_kb, scale=2)

    wf_dist_hfm = optical_elements.WF_dist(1500, 100, kb_clear_ap, kb_clear_ap)
    defineOPD(wf_dist_hfm, os.path.join(mirror_data_dir,'mirror1.dat'), 2, '\t', 'x',
              theta_kb, scale=2, stretching=kb_mirror_length/0.8)

    wf_dist_vfm = optical_elements.WF_dist(1100, 1500, kb_clear_ap, kb_clear_ap)
    defineOPD(wf_dist_vfm, os.path.join(mirror_data_dir,'mirror2.dat'), 2, ' ', 'y',
              theta_kb, scale=2, stretching=kb_mirror_length/0.8)


    # Assemble the beamline with PP parameters.
    bl0 = Beamline()
    bl0.append(ap0,   Use_PP(semi_analytical_treatment=0,
                             zoom=14.4,
                             sampling=1/1.6))
    bl0.append(drift0,Use_PP(semi_analytical_treatment=0))
    bl0.append(ap1, Use_PP(zoom=0.8))
    bl0.append(wf_dist_om, Use_PP())
    bl0.append(drift1, Use_PP(semi_analytical_treatment=1))
    bl0.append(ap_kb,  Use_PP(zoom = 6.4, sampling = 1/16.))
    bl0.append(hfm, Use_PP())
    bl0.append(wf_dist_hfm, Use_PP())
    bl0.append(drift_in_kb, Use_PP(semi_analytical_treatment=1))
    bl0.append(vfm, Use_PP())
    bl0.append(wf_dist_vfm, Use_PP())
    bl0.append(drift_to_foc, Use_PP(semi_analytical_treatment=1))

    # All done, return.
    return bl0
