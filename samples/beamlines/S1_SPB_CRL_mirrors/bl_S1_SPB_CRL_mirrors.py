
def get_beamline():
    import os
    import wpg
    from wpg import Beamline
    from wpg.optical_elements import Aperture, Drift, CRL, Empty, Use_PP, Mirror_plane

    wpg_path = os.path.abspath(os.path.dirname(wpg.__file__))

    # S1 beamline layout
    # Geometry ###
    src_to_hom1 = 257.8  # Distance source to HOM 1 [m]
    src_to_hom2 = 267.8  # Distance source to HOM 2 [m]
    src_to_crl = 887.8  # Distance source to CRL [m]
    #     src_to_exp = 920.42 # Distance source to experiment [m]

    # Drift to focus aperture
    # crl_to_exp_drift = Drift( src_to_exp - src_to_crl )
    z = 34.0
    # define distances, angles, etc
    # ...
    # Incidence angle at HOM

      # should be checked for other beams !!!

    theta_om = 3.6e-3  # [rad]

    om_mirror_length = 0.8  # [m]
    om_clear_ap = om_mirror_length * theta_om

    # define the beamline:
    bl0 = Beamline()
    zoom = 1

    # Define HOM1.
    aperture_x_to_y_ratio = 1
    hom1 = Aperture(
        shape='r', ap_or_ob='a', Dx=om_clear_ap, Dy=om_clear_ap / aperture_x_to_y_ratio)
    bl0.append(
        hom1, Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom))

    # Define and apply mirror distortion profile.
    mirrors_path = os.path.join(wpg_path, '..', 'samples', 'data_common')
    hom1_wavefront_distortion = Mirror_plane(orient='x', 
                                             theta=theta_om, 
                                             length=om_mirror_length, 
                                             range_xy=om_clear_ap/aperture_x_to_y_ratio, 
                                             filename=os.path.join(
                                             mirrors_path, 'mirror1.dat'), 
                                             scale=1.,
                                             bPlot=True)
    bl0.append(hom1_wavefront_distortion,
               Use_PP(semi_analytical_treatment=0, zoom=zoom, sampling=zoom))

    # Free space propagation from hom1 to hom2
    hom1_to_hom2_drift = Drift(src_to_hom2 - src_to_hom1)
    bl0.append(hom1_to_hom2_drift, Use_PP(semi_analytical_treatment=0))

    # Define HOM2.
    zoom = 1.0
    hom2 = Aperture('r', 'a', om_clear_ap, om_clear_ap / aperture_x_to_y_ratio)
    bl0.append(hom2, Use_PP(semi_analytical_treatment=0,
                            zoom=zoom, sampling=zoom / 0.75))

    # define and apply mirror distortion.
    hom2_wavefront_distortion = Mirror_plane(orient='x', 
                                             theta=theta_om, 
                                             length=om_mirror_length, 
                                             range_xy=om_clear_ap/aperture_x_to_y_ratio, 
                                             filename=os.path.join(
                                             mirrors_path, 'mirror2.dat'), 
                                             scale=1.,
                                             xscale=1.,
                                             bPlot=True)
    bl0.append(hom2_wavefront_distortion, Use_PP(
        semi_analytical_treatment=0, zoom=zoom, sampling=zoom))

    # drift to CRL aperture
    hom2_to_crl_drift = Drift(src_to_crl - src_to_hom2)

    bl0.append(hom2_to_crl_drift, Use_PP(semi_analytical_treatment=1))

    # Define CRL
    crl_focussing_plane = 3  # Both horizontal and vertical.
    # Refractive index decrement (n = 1- delta - i*beta)
    crl_delta = 4.7177e-06
    crl_attenuation_length = 6.3e-3    # Attenuation length [m], Henke data.
    crl_shape = 1         # Parabolic lenses
    crl_aperture = 5.0e-3  # [m]
    crl_curvature_radius = 5.8e-3  # [m]
    crl_number_of_lenses = 19
    crl_wall_thickness = 8.0e-5  # Thickness
    crl_center_horizontal_coordinate = 0.0
    crl_center_vertical_coordinate = 0.0
    crl_initial_photon_energy = 8.48e3  # [eV] ### OK ???
    crl_final_photon_energy = 8.52e3  # [eV]   ### OK ???

    crl = CRL(_foc_plane=crl_focussing_plane,
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

    bl0.append(
        crl, Use_PP(semi_analytical_treatment=1, zoom=zoom, sampling=zoom/0.1))

    crl_to_exp_drift = Drift(z)
    bl0.append(crl_to_exp_drift, Use_PP(
        semi_analytical_treatment=1, zoom=1, sampling=1))
    #     bl0.append(Empty(),Use_PP(zoom=0.25, sampling=0.25))

    return bl0