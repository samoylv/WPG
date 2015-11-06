:mod:`wpg.optical_elements` module
----------------------------------

.. automodule:: wpg.optical_elements
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: wpg.optical_elements.Drift
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: wpg.optical_elements.Lens
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: wpg.optical_elements.Aperture
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: wpg.optical_elements.Mirror_elliptical
    :members:
    :undoc-members:
    :show-inheritance:

.. autoclass:: wpg.optical_elements.WF_dist
    :members:
    :undoc-members:
    :show-inheritance:

Some useful redefined optical elements
--------------------------------------

.. function:: defineEFM(orient,p,q,thetaEFM,theta0,lengthEFM):

    A wrapper to a SRWL function SRWLOptMirEl() for defining a plane elliptical focusing mirror propagator
    
    :param Orient:    mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param p:  the distance to two ellipsis centers
    :param q:  the distance to two ellipsis centers
    :param thetaEFM:  the design incidence angle in the center of the mirror
    :param theta0:    the "real" incidence angle in the center of the mirror
    :param lengthEFM: mirror length, [m]
    :return: the struct opEFM

.. code-block:: python

    def defineEFM(orient,p,q,thetaEFM,theta0,lengthEFM):
    """
    A wrapper to a SRWL function SRWLOptMirEl() for defining a plane elliptical focusing mirror propagator
    
    :param Orient:    mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param p:  the distance to two ellipsis centers
    :param q:  the distance to two ellipsis centers
    :param thetaEFM:  the design incidence angle in the center of the mirror
    :param theta0:    the "real" incidence angle in the center of the mirror
    :param lengthEFM: mirror length, [m]
    :return: the struct opEFM
    """

    if orient == 'x':     #horizontal plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(_p=p, _q=q, _ang_graz=thetaEFM, _r_sag=1.e+40, _size_tang=lengthEFM, 
            _nvx=cos(theta0), _nvy=0, _nvz=-sin(theta0), _tvx=-sin(theta0), _tvy=0, _x=0, _y=0, _treat_in_out=1) 
    elif orient == 'y': #vertical plane ellipsoidal mirror
        opEFM = SRWLOptMirEl(_p=p, _q=q, _ang_graz=thetaEFM, _r_sag=1.e+40, _size_tang=lengthEFM, 
            _nvx=0, _nvy=cos(theta0), _nvz=-sin(theta0), _tvx=0, _tvy=-sin(theta0), _x=0, _y=0, _treat_in_out=1)
    else:
        raise TypeError('orient should be "x" or "y"')
    return opEFM

