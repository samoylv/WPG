:mod:`wpg.materials` module
---------------------------

Optical constants from xraylib (``pip install wpg[materials]``) and reflectivity of mirrors
for ``reflectivity=`` of ``Mirror_elliptical``, ``Mirror_hyperbolic``, ``Mirror_parabolic``,
``Mirror_spherical`` and ``Mirror_toroidal``.

.. code-block:: python

    import numpy as np
    from wpg.materials import Reflectivity
    from wpg.optical_elements import Mirror_elliptical

    theta = 3e-3
    refl = Reflectivity.from_material('B4C', ekev=np.linspace(8.9, 9.1, 41),
                                      angles=np.linspace(0.8 * theta, 1.2 * theta, 201), rho=2.52)
    mirror = Mirror_elliptical('x', p=50., q=5., thetaE=theta, theta0=theta, length=0.5,
                               reflectivity=refl)

.. automodule:: wpg.materials
    :members:
    :undoc-members:
