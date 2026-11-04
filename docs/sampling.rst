:mod:`wpg.sampling` module
--------------------------

Sampling of a wavefront and the choice of ``Use_PP`` zoom and sampling for an optical element.

.. code-block:: python

    from wpg.sampling import check_sampling, set_radius_from_field, suggest_propagation_parameters

    print(check_sampling(wf))       # phase step between pixels, window border, SRW radius vs field
    set_radius_from_field(wf)       # if the report says the SRW radius disagrees with the field
    pp, report = suggest_propagation_parameters(wf, 50.)   # drift of 50 m, trial propagations of one slice
    bl.append(Drift(50.), pp)

.. automodule:: wpg.sampling
    :members:
