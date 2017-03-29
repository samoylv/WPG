# Module holding the SPB/SFX beamline (nanofocus, KB mirrors) at European XFEL.

from wpg import optical_elements, Beamline
from wpg.optical_elements import Use_PP, Aperture, Drift, CRL
from opd import defineOPD
import os


def get_beamline():
    """ Setup and return the WPG.Beamline object representing the SPB/SFX nanofocus beamline (KB mirrors).

    :return: beamline
    :rtype beamline: wpg.Beamline
    """

    ### Geometry ###
    src_to_hom1 = 257.8 # Distance source to HOM 1 [m]
    src_to_hom2 = 267.8 # Distance source to HOM 2 [m]
    src_to_crl = 887.8  # Distance source to CRL [m]
    src_to_exp = 920.42 # Distance source to experiment [m]

    #Incidence angle at HOM
    theta_om = 3.6e-3       # [rad]

    om_mirror_length = 0.8 # [m]
    om_clear_ap = om_mirror_length*theta_om


    #define the beamline:
    beamline = Beamline()
    zoom=1

    # Define HOM1 = Aperture + Wavefront distortion.
    aperture_x_to_y_ratio = 1
    hom1_aperture = Aperture(shape='r',ap_or_ob='a',Dx=om_clear_ap,Dy=om_clear_ap/aperture_x_to_y_ratio)

    # Append to beamline.
    beamline.append( hom1_aperture, Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom) )

    # Free space propagation from hom1 to hom2
    hom1_to_hom2_drift = Drift(src_to_hom2 - src_to_hom1)
    beamline.append( hom1_to_hom2_drift, Use_PP(semi_analytical_treatment=0))

    # Define HOM2.
    zoom = 1.0
    hom2_aperture = Aperture('r','a', om_clear_ap, om_clear_ap/aperture_x_to_y_ratio)
    beamline.append( hom2_aperture,  Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom))

    #drift to CRL aperture
    hom2_to_crl_drift = Drift( src_to_crl - src_to_hom2 )
    beamline.append( hom2_to_crl_drift, Use_PP(semi_analytical_treatment=1))

    # Circular Aperture before CRL.
    crl_front_aperture_diameter = 2.8e-3
    crl_front_aperture = Aperture('c','a', crl_front_aperture_diameter, crl_front_aperture_diameter)

    ### Define CRL
    crl_focussing_plane = 3 # Both horizontal and vertical.
    crl_delta = 4.8308e-06 # Refractive index decrement (n = 1- delta - i*beta) @ 8.4 keV
    crl_attenuation_length  = 6.053e-3    # Attenuation length [m], Henke data.
    crl_shape = 1         # Parabolic lenses
    crl_aperture = 3.0e-3 # [m]
    crl_curvature_radius = 5.8e-3 # [m]
    crl_number_of_lenses = 19
    crl_wall_thickness = 8.0e-5 # Thickness
    crl_center_horizontal_coordinate = 0.0
    crl_center_vertical_coordinate = 0.0
    crl_initial_photon_energy = 8.48e3 # [eV]
    crl_final_photon_energy = 8.52e3 # [eV]

    crl = CRL( _foc_plane=crl_focussing_plane,
               _delta=crl_delta,
               _atten_len=crl_attenuation_length,
               _shape=crl_shape,
               _apert_h=crl_aperture,
               _apert_v=crl_aperture,
               _r_min=crl_curvature_radius,
               _n=crl_number_of_lenses,
               _wall_thick=crl_wall_thickness,
               _xc=crl_center_horizontal_coordinate,
               _yc=crl_center_vertical_coordinate,
               _void_cen_rad=None,
               _e_start=crl_initial_photon_energy,
               _e_fin=crl_final_photon_energy,
               )

    zoom = 0.6
    beamline.append( crl_front_aperture, Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom/0.1) )
    beamline.append( crl, Use_PP(semi_analytical_treatment=0, zoom=1, sampling=1) )

    # Drift to focus aperture
    crl_to_exp_drift = Drift( src_to_exp - src_to_crl )
    beamline.append( crl_to_exp_drift, Use_PP(semi_analytical_treatment=1, zoom=1, sampling=1) )

    return beamline




