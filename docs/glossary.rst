Glossary definition
-------------------

Current gloaasry contains the follaing fields:


**params/wSpace** - Real space or  q-space WF presentation [string] - ***string***

**params/nval** - complex electric field nval==2 - ******

**params/Mesh/ny** - Numbers of points, vertical - ******

**params/Mesh/nx** - Numbers of points, horizontal - ******

**params/Mesh/nSlices** - Numbers of points vs photon energy/time for the pulse - ******

**data/arrEhor** -      EM field (Re, Im) pairs written in 3D array, slice number changes first.     Horizontal polarization      - ******

**data/arrEver** -      EM field (Re, Im) pairs written in 3D array, slice number changes first.     Vertical polarization      - ******

**params/dRx** - Error of wavefront horizontal radius [m] - ***m***

**params/dRy** - Error of wavefront horizontal radius [m] - ***m***

**params/Mesh/qxMax** - Maximum of horizontal frequency [1/m] - ***1/m***

**params/Mesh/qxMin** - Minimum of horizontal frequency [1/m] - ***1/m***

**params/Mesh/qyMax** - Maximum of vertical frequency [1/m] - ***1/m***

**params/Mesh/qyMin** - Minimum of vertical frequency [1/m] - ***1/m***

**params/Mesh/sliceMax** - Max value of time [s] or energy [ev] for pulse (fragment) - ***s or ev***

**params/Mesh/sliceMin** - Min value of time [s] or energy [ev] for pulse (fragment) - ***s or ev***

**params/Mesh/xMax** - Maximum of horizontal range [m] - ***m***

**params/Mesh/xMin** - Minimum of horizontal range [m] - ***m***

**params/Mesh/yMax** - Maximum of vertical range [m] - ***m***

**params/Mesh/yMin** - Minimum of vertical range [m] - ***m***

**params/Mesh/zCoord** - Longitudinal position [m], Fast data: length of  active undulator, Gaussian source: distance to waist  - ***m***

**params/photonEnergy** - Average photon energy [ev] - ***ev***

**params/Rx** - Instantaneous horizontal wavefront radius [m] - ***m***

**params/Ry** - Instantaneous vertical wavefront radius [m] - ***m***

**params/wDomain** - WF in time or frequency (photon energy) domain [string] - ***string***

**params/wEFieldUnit** - Electric field units [string] - ***string***

**params/wFloatType** - Electric field numerical type [string] - ***string***

**params/xCentre** - Horizontal transverse coordinates of wavefront instant 'source center' [m] - ***m***

**params/yCentre** - Vertical transverse coordinates of wavefront instant 'source center' [m] - ***m***

**version** - Hdf5 format version (glossary) - ******

:mod:`wpg.glossary` module
----------------------------

.. automodule:: wpg.glossary
    :members:
    :undoc-members:
    :show-inheritance: