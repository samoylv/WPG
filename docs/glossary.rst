Glossary definition
-------------------

The glossary (:mod:`wpg.glossary`) defines the fields of a wavefront: mapping of ``srwlib.SRWLWfr``
to :class:`wpg.Wavefront` and the names of the datasets in HDF5 files.
Each field is declared once in a node class: ``wf.params`` (:class:`~wpg.glossary.WavefrontParams`),
``wf.params.Mesh`` (:class:`~wpg.glossary.WavefrontMesh`) and ``wf.data`` (:class:`~wpg.glossary.WavefrontData`).
The fields read and write the SRW wavefront ``wf._srwl_wf`` directly.
The HDF5 dataset name of a field is its glossary path (shown in brackets),
the description, units, limits and SRW alias of the field are stored as attributes of the dataset.
The list of all fields with units is printed by ``wpg.print_glossary()``.

.. autoclass:: wpg.glossary.WavefrontParams
    :members:
    :undoc-members:

.. autoclass:: wpg.glossary.WavefrontMesh
    :members:
    :undoc-members:

.. autoclass:: wpg.glossary.WavefrontData
    :members:
    :undoc-members:

Besides the fields of the nodes, the HDF5 file contains ``version`` (format version, ``wf.version``)
and the user fields of ``wf.custom_fields`` (e.g. ``history``, ``misc``).

.. automodule:: wpg.glossary
    :members: Field, get_wf_fields, print_glossary
