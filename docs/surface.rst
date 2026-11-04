:mod:`wpg.surface` module
-------------------------

Surface height errors of mirrors as transmission elements, in the model of SRW
(``srwl_opt_setup_surf_height_1d/2d``), and random surfaces with a power-law spectrum.

.. code-block:: python

    from wpg.optical_elements import Mirror_elliptical, Use_PP
    from wpg.surface import Surface_errors, random_profile

    theta = 3e-3
    x, h = random_profile(length=0.8, n=4001, sigma=2e-9, exponent=-2.5, seed=1)
    bl.append(Mirror_elliptical('x', p=50., q=5., thetaE=theta, theta0=theta, length=0.8), Use_PP())
    bl.append(Surface_errors('x', theta, (x, h), width=1.2e-3), Use_PP())

.. automodule:: wpg.surface
    :members:
