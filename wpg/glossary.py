# -*- coding: utf-8 -*-

"""
This module contains definitions (glossary) of Wavefront fields: mapping of fields SRWLWfr <-> wpg.Wavefront
and their description in HDF5 files.

Each field is declared once as a class attribute (:class:`Field` descriptor) of a wavefront node:
``wf.params`` (:class:`WavefrontParams`), ``wf.params.Mesh`` (:class:`WavefrontMesh`) and
``wf.data`` (:class:`WavefrontData`). The descriptor reads and writes the SRW structure ``wf._srwl_wf``
directly, so the fields are a live view of the SRW wavefront, visible to IDEs and static analysis.
The glossary path of the field (e.g. ``params/Mesh/nx``) is the dataset name in the HDF5 file,
the field metadata (description, units, limits, alias) are stored as the dataset attributes.

.. module:: wpg.glossary
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

import array
import operator
import warnings
from typing import Any, Callable, Generic, Iterator, List, NoReturn, Optional, TypeVar, overload

import numpy as np

T = TypeVar('T')


def _units_from_description(description):
    """Units written in square brackets in the description: 'time [s] or energy [ev]' -> 's or ev'."""
    units = []
    start = stop = 0
    while start >= 0 and stop >= 0:
        start = description.find('[', stop)
        stop = description.find(']', start)
        if start >= 0 and stop >= 0:
            units.append(description[start + 1: stop])
    return ' or '.join(units)


class Field(Generic[T]):
    """
    Wavefront field: descriptor of a node class mapped to an attribute of SRWLWfr and to an HDF5 dataset.

    :param path: glossary path, also the dataset name in HDF5 files, e.g. ``'params/Mesh/nx'``
    :param description: field description (HDF5 attribute ``description``)
    :param srw: attribute of SRWLWfr, e.g. ``'mesh.nx'``, for fields directly mapped to SRW
    :param convert: type conversion applied on assignment to ``srw`` (``int``, ``float``)
    :param get: ``get(wf)`` for fields not directly mapped to SRW
    :param set: ``set(wf, value)`` for fields not directly mapped to SRW
    :param units: units (HDF5 attribute ``units``), by default taken from ``[...]`` in the description
    :param limits: allowed values (HDF5 attribute ``limits``)
    :param alias: name in SRW (HDF5 attribute ``alias``)
    """

    def __init__(self, path: str, description: str, *, srw: Optional[str] = None,
                 convert: Optional[Callable[[Any], Any]] = None,
                 get: Optional[Callable[[Any], T]] = None,
                 set: Optional[Callable[[Any, T], None]] = None,
                 units: Optional[str] = None, limits: Optional[str] = None, alias: Optional[str] = None):
        if (srw is None) == (get is None):
            raise ValueError(f'{path}: exactly one of srw= and get= must be given')
        self.path = path
        self.description = description
        self.srw = srw
        self._convert = convert
        self._get = get
        self._set = set
        self.attributes = {'description': description,
                           'units': _units_from_description(description) if units is None else units}
        if limits is not None:
            self.attributes['limits'] = limits
        if alias is not None:
            self.attributes['alias'] = alias
        self.__doc__ = f"{description} (HDF5: ``{path}``)"
        self.name = path.split('/')[-1]

    @property
    def glossary_name(self) -> str:
        """Glossary path of the field (name used in older versions of WPG)."""
        return self.path

    def __set_name__(self, owner, name):
        self.name = name

    @overload
    def __get__(self, node: None, owner: Any = None) -> 'Field[T]': ...

    @overload
    def __get__(self, node: Any, owner: Any = None) -> T: ...

    def __get__(self, node, owner=None):
        if node is None:
            return self
        wf = node._wf
        if self._get is not None:
            return self._get(wf)
        return operator.attrgetter(self.srw)(wf._srwl_wf)

    def __set__(self, node, value: T) -> None:
        wf = node._wf
        if self._get is not None:
            if self._set is None:
                raise AttributeError(f'{self.path} is read-only')
            self._set(wf, value)
            return
        parent, _, attr = self.srw.rpartition('.')
        obj = operator.attrgetter(parent)(wf._srwl_wf) if parent else wf._srwl_wf
        setattr(obj, attr, self._convert(value) if self._convert is not None else value)

    def __repr__(self):
        return f'Field({self.path!r})'


class ReadOnlyField(Field[T]):
    """Field that can not be assigned (static analyzers report an assignment as an error)."""

    def __set__(self, node, value: NoReturn) -> None:
        raise AttributeError(f'{self.name} is a read-only shortcut, assign wf.params... instead')


class _Node:
    """Base class of wavefront nodes (``wf.params``, ``wf.params.Mesh``, ``wf.data``)."""

    __slots__ = ('_wf',)

    def __init__(self, wf):
        self._wf = wf

    @classmethod
    def fields(cls) -> List[Field]:
        """Fields of this node in the declaration order."""
        res = []
        for klass in reversed(cls.__mro__):
            res += [v for v in vars(klass).values() if isinstance(v, Field)]
        return res

    def __repr__(self):
        values = ', '.join(f'{f.name}={getattr(self, f.name)!r}' for f in self.fields())
        return f'{type(self).__name__}({values})'


# ---------------------------------------------------------------- helpers for non-trivial fields

def _enum_field(srw_attr, values, name):
    """get/set for a string field mapped to an integer (or char) SRW flag: values = {srw_value: string}."""
    def get(wf):
        v = getattr(wf._srwl_wf, srw_attr)
        if isinstance(v, bytes):
            v = v.decode()
        try:
            return values[v]
        except KeyError:
            raise ValueError(f'internal error, wrong wavefront field value {name}={v!r}') from None

    def set(wf, val):
        if isinstance(val, bytes):
            val = val.decode()
        for k, v in values.items():
            if v == val:
                setattr(wf._srwl_wf, srw_attr, k)
                return
        raise ValueError(f'{name}: value must be one of {list(values.values())}, given: {val!r}')

    return get, set


def _get_w_float_type(wf):
    v = wf._srwl_wf.numTypeElFld
    if v in ('f', b'f'):
        return 'float'
    elif v in ('d', b'd'):
        return 'double'
    raise ValueError(f'internal error, wrong wavefront field value numTypeElFld={v!r}')


def _set_w_float_type(wf, val):
    if val in ('float', b'float'):
        wf._srwl_wf.numTypeElFld = 'f'
    elif val in ('double', b'double'):
        raise ValueError('"double" type not supported yet')
    else:
        raise ValueError('value must be "float" or "double"')


_E_FIELD_UNITS = {
    'arbitrary': 0,
    'sqrt(Phot/s/0.1%bw/mm^2)': 1,
    'sqrt(J/eV/mm^2)': 2,  # frequency domain
    'sqrt(W/mm^2)': 2,     # time domain
}


def _get_w_efield_unit(wf):
    unit = wf._srwl_wf.unitElFld
    if unit == 0:
        return 'arbitrary'
    elif unit == 1:
        return 'sqrt(Phot/s/0.1%bw/mm^2)'
    elif unit == 2:
        if wf._srwl_wf.presFT == 0:
            return 'sqrt(J/eV/mm^2)'
        elif wf._srwl_wf.presFT == 1:
            return 'sqrt(W/mm^2)'
    raise ValueError(f'internal error, wrong wavefront field value unitElFld={unit!r}')


def _set_w_efield_unit(wf, val):
    if isinstance(val, bytes):
        val = val.decode()
    if val not in _E_FIELD_UNITS:
        raise ValueError(f'value must be one of {list(_E_FIELD_UNITS)}, given: {val!r}')
    wf._srwl_wf.unitElFld = _E_FIELD_UNITS[val]


def _space_dependent_range(srw_attr, space, path):
    """get=, set= for mesh range defined only in R-space (xMin, ...) or Q-space (qxMin, ...)."""
    def get(wf):
        if wf.params.wSpace == space:
            return getattr(wf._srwl_wf.mesh, srw_attr)
        warnings.warn(f'{path} not defined if NOT params/wSpace=={space}')
        return None

    def set(wf, val):
        if not wf.params.wSpace == space:
            warnings.warn(f'{path} not defined if NOT params/wSpace=={space}')
        setattr(wf._srwl_wf.mesh, srw_attr, float(val))

    return {'get': get, 'set': set}


def _set_n_slices(wf, val):
    wf._srwl_wf.mesh.ne = int(val)
    wf._allocate_srw_moments()


def _set_nval(wf, val):
    wf._nval = int(val)


def _set_version(wf, val):
    wf._version = val


def _efield_array(srw_attr):
    """get/set of EM field array (Re, Im pairs) of one polarization as a 4D numpy array (ny, nx, nSlices, nval)."""
    def get(wf):
        res = np.asarray(getattr(wf._srwl_wf, srw_attr), dtype='float32')
        mesh = wf.params.Mesh
        correct_shape = (mesh.ny, mesh.nx, mesh.nSlices, wf.params.nval)
        if not res.shape:
            pass
        elif res.size == np.prod(correct_shape):
            res = res.reshape(correct_shape)
        else:
            raise ValueError(f'size of SRW field array {res.size} does not match the wavefront mesh {correct_shape}')
        if not res.flags['C_CONTIGUOUS']:
            res = np.ascontiguousarray(res)
        return res

    def set(wf, val):
        """val: 4D array (ny, nx, nSlices, nval) or flat array / array.array of the mesh size"""
        n_total = wf._get_total_elements() * wf.params.nval
        if isinstance(val, array.array) and len(val) != n_total:
            warnings.warn('New array size not equal to wavefront size. You must set it by yourself.')
        val = np.ascontiguousarray(np.array(val, dtype='float32'))
        if val.size != n_total:
            if val.ndim == 4:
                warnings.warn('New array size not equal to wavefront size. '
                              'Mesh (ny, nx, nSlices) is set from the array shape.')
                mesh = wf.params.Mesh
                mesh.ny, mesh.nx, mesh.nSlices = val.shape[:3]
                other = 'arEy' if srw_attr == 'arEx' else 'arEx'
                other_arr = getattr(wf._srwl_wf, other)  # None in a fresh SRWLWfr
                if other_arr is None or len(other_arr) != val.size:
                    warnings.warn(f'Field of the other polarization ({other}) does not match the new mesh, '
                                  'it is set to zeros.')
                    setattr(wf._srwl_wf, other, np.zeros(val.size, dtype='float32'))
            else:
                warnings.warn('New array size not equal to wavefront size. You must set it by yourself.')
        wf._allocate_srw_moments()
        setattr(wf._srwl_wf, srw_attr, val.reshape(-1))

    return get, set


# ---------------------------------------------------------------- wavefront nodes

class WavefrontMesh(_Node):
    """Mesh of the wavefront (``wf.params.Mesh``), view of ``SRWLWfr.mesh``."""

    __slots__ = ()

    # mesh size first: it is needed to load field arrays and ranges from files
    nx = Field[int]('params/Mesh/nx', 'Numbers of points, horizontal', srw='mesh.nx', convert=int,
                    units='-', limits='[2:LONG_MAX]', alias='')
    ny = Field[int]('params/Mesh/ny', 'Numbers of points, vertical', srw='mesh.ny', convert=int,
                    units='-', limits='[2:LONG_MAX]', alias='')
    nSlices = Field[int]('params/Mesh/nSlices', 'Numbers of points vs photon energy/time for the pulse',
                         get=lambda wf: wf._srwl_wf.mesh.ne, set=_set_n_slices,
                         units='', limits='[1:LONG_MAX]', alias='mesh.ne')
    zCoord = Field[float]('params/Mesh/zCoord',
                          'Longitudinal position [m], Fast data: length of active undulator, '
                          'Gaussian source: distance to waist',
                          srw='mesh.zStart', convert=float, limits='{???}', alias='mesh.zStart')
    xMin = Field[Optional[float]]('params/Mesh/xMin', 'Minimum of horizontal range [m]',
                                  **_space_dependent_range('xStart', 'R-space', 'params/Mesh/xMin'))
    xMax = Field[Optional[float]]('params/Mesh/xMax', 'Maximum of horizontal range [m]',
                                  **_space_dependent_range('xFin', 'R-space', 'params/Mesh/xMax'))
    yMin = Field[Optional[float]]('params/Mesh/yMin', 'Minimum of vertical range [m]',
                                  **_space_dependent_range('yStart', 'R-space', 'params/Mesh/yMin'))
    yMax = Field[Optional[float]]('params/Mesh/yMax', 'Maximum of vertical range [m]',
                                  **_space_dependent_range('yFin', 'R-space', 'params/Mesh/yMax'))
    qxMin = Field[Optional[float]]('params/Mesh/qxMin', 'Minimum of horizontal frequency [1/m]',
                                   **_space_dependent_range('xStart', 'Q-space', 'params/Mesh/qxMin'))
    qxMax = Field[Optional[float]]('params/Mesh/qxMax', 'Maximum of horizontal frequency [1/m]',
                                   **_space_dependent_range('xFin', 'Q-space', 'params/Mesh/qxMax'))
    qyMin = Field[Optional[float]]('params/Mesh/qyMin', 'Minimum of vertical frequency [1/m]',
                                   **_space_dependent_range('yStart', 'Q-space', 'params/Mesh/qyMin'))
    qyMax = Field[Optional[float]]('params/Mesh/qyMax', 'Maximum of vertical frequency [1/m]',
                                   **_space_dependent_range('yFin', 'Q-space', 'params/Mesh/qyMax'))
    sliceMin = Field[float]('params/Mesh/sliceMin', 'Min value of time [s] or energy [ev] for pulse (fragment)',
                            srw='mesh.eStart', convert=float, limits='[1:LONG_MAX]', alias='mesh.eStart')
    sliceMax = Field[float]('params/Mesh/sliceMax', 'Max value of time [s] or energy [ev] for pulse (fragment)',
                            srw='mesh.eFin', convert=float, limits='[1:LONG_MAX]', alias='mesh.eFin')
    # components of unit vectors of the observation plane (float; before 2026.11 they were truncated to int)
    hvx = Field[float]('params/Mesh/hvx', 'Lab-frame horizontal base vector of the observation plane '
                       '(/ surface in its center)', srw='mesh.hvx', convert=float, units='-', limits='[-1:1]', alias='')
    hvy = Field[float]('params/Mesh/hvy', 'Lab-frame horizontal base vector of the observation plane '
                       '(/ surface in its center)', srw='mesh.hvy', convert=float, units='-', limits='[-1:1]', alias='')
    hvz = Field[float]('params/Mesh/hvz', 'Lab-frame horizontal base vector of the observation plane '
                       '(/ surface in its center)', srw='mesh.hvz', convert=float, units='-', limits='[-1:1]', alias='')
    nvx = Field[float]('params/Mesh/nvx', 'Lab-frame coordinate of the inner normal to observation plane '
                       '(/ surface in its center)', srw='mesh.nvx', convert=float, units='-', limits='[-1:1]', alias='')
    nvy = Field[float]('params/Mesh/nvy', 'Lab-frame coordinate of the inner normal to observation plane '
                       '(/ surface in its center)', srw='mesh.nvy', convert=float, units='-', limits='[-1:1]', alias='')
    nvz = Field[float]('params/Mesh/nvz', 'Lab-frame coordinate of the inner normal to observation plane '
                       '(/ surface in its center)', srw='mesh.nvz', convert=float, units='-', limits='[-1:1]', alias='')


_get_w_domain, _set_w_domain = _enum_field('presFT', {0: 'frequency', 1: 'time'}, 'params/wDomain')
_get_w_space, _set_w_space = _enum_field('presCA', {0: 'R-space', 1: 'Q-space'}, 'params/wSpace')


class WavefrontParams(_Node):
    """Wavefront parameters (``wf.params``), view of SRWLWfr."""

    __slots__ = ('Mesh',)

    Mesh: WavefrontMesh

    def __init__(self, wf):
        super().__init__(wf)
        self.Mesh = WavefrontMesh(wf)

    # nval and wSpace first: they are needed to load field arrays and ranges from files
    nval = Field[int]('params/nval', 'complex electric field nval==2',
                      get=lambda wf: wf._nval, set=_set_nval, units='', limits='{2}', alias='-')
    wSpace = Field[str]('params/wSpace', 'Real space or q-space WF presentation [string]',
                        get=_get_w_space, set=_set_w_space, limits='{"R-space", "Q-space"}', alias='presCA')
    wDomain = Field[str]('params/wDomain', 'WF in time or frequency (photon energy) domain [string]',
                         get=_get_w_domain, set=_set_w_domain, limits='{"time", "frequency"}', alias='presFT')
    wFloatType = Field[str]('params/wFloatType', 'Electric field numerical type [string]',
                            get=_get_w_float_type, set=_set_w_float_type,
                            limits='{"float", "double"}', alias='numTypeElFld')
    wEFieldUnit = Field[str]('params/wEFieldUnit', 'Electric field units [string]',
                             get=_get_w_efield_unit, set=_set_w_efield_unit,
                             limits='{"sqrt(Phot/s/0.1%bw/mm^2)", "sqrt(W/mm^2)", "sqrt(J/eV/mm^2)", "arbitrary"}',
                             alias='unitElFld')
    photonEnergy = Field[float]('params/photonEnergy', 'Average photon energy [ev]',
                                srw='avgPhotEn', convert=float, limits='[1:1e6]', alias='avgPhotEn')
    Rx = Field[float]('params/Rx', 'Instantaneous horizontal wavefront radius [m]',
                      srw='Rx', convert=float, limits='[FLOAT_MIN:FLOAT_MAX]', alias='Rx')
    Ry = Field[float]('params/Ry', 'Instantaneous vertical wavefront radius [m]',
                      srw='Ry', convert=float, limits='[FLOAT_MIN:FLOAT_MAX]', alias='Ry')
    dRx = Field[float]('params/dRx', 'Error of wavefront horizontal radius [m]',
                       srw='dRx', convert=float, limits='[FLOAT_MIN:FLOAT_MAX]', alias='dRx')
    dRy = Field[float]('params/dRy', 'Error of wavefront vertical radius [m]',
                       srw='dRy', convert=float, limits='[FLOAT_MIN:FLOAT_MAX]', alias='dRy')
    xCentre = Field[float]('params/xCentre',
                           "Horizontal transverse coordinates of wavefront instant 'source center' [m]",
                           srw='xc', convert=float)
    yCentre = Field[float]('params/yCentre',
                           "Vertical transverse coordinates of wavefront instant 'source center' [m]",
                           srw='yc', convert=float)


_get_arr_ehor, _set_arr_ehor = _efield_array('arEx')
_get_arr_ever, _set_arr_ever = _efield_array('arEy')


class WavefrontData(_Node):
    """EM field arrays of the wavefront (``wf.data``)."""

    __slots__ = ()

    arrEhor = Field[np.ndarray]('data/arrEhor', 'EM field (Re, Im) pairs written in 3D array, '
                                'slice number changes first. Horizontal polarization',
                                get=_get_arr_ehor, set=_set_arr_ehor,
                                units='see params/wEFieldUnit', limits='[FLOAT_MIN:FLOAT_MAX]', alias='arEx')
    arrEver = Field[np.ndarray]('data/arrEver', 'EM field (Re, Im) pairs written in 3D array, '
                                'slice number changes first. Vertical polarization',
                                get=_get_arr_ever, set=_set_arr_ever,
                                units='see params/wEFieldUnit', limits='[FLOAT_MIN:FLOAT_MAX]', alias='arEy')


version_field = Field[float]('version', 'Hdf5 format version (glossary)',
                             get=lambda wf: wf._version, set=_set_version,
                             units='-', limits='[???]', alias='VERSION')


def get_wf_fields() -> Iterator[Field]:
    """
    All wavefront fields in the order of loading from files: mesh size and space first, field arrays last.

    :return: iterator over :class:`Field`
    """
    yield from WavefrontParams.fields()[:2]  # nval, wSpace
    yield from WavefrontMesh.fields()
    yield from WavefrontParams.fields()[2:]
    yield version_field
    yield from WavefrontData.fields()


def field_value(wf, field: Field):
    """Value of the field in the wavefront wf."""
    return field.__get__(field_node(wf, field))


def shortcut(field: Field[T]) -> ReadOnlyField[T]:
    """
    Read-only alias of a field of ``wf.params`` or ``wf.params.Mesh`` at the top level of the wavefront:
    ``wf.nx`` is ``wf.params.Mesh.nx``. Not stored in HDF5 files (not in :func:`get_wf_fields`).
    """
    node = 'params.Mesh' if field.path.startswith('params/Mesh/') else 'params'
    res = ReadOnlyField[T](field.path, field.description, units=field.attributes['units'],
                           get=lambda wf: field.__get__(field_node(wf, field)))
    res.__doc__ = f"Read-only shortcut of ``wf.{node}.{field.name}``: {field.description}"
    return res


def field_node(wf, field: Field):
    """Node of the wavefront wf the field belongs to."""
    if field is version_field:
        return wf
    if field.path.startswith('data/'):
        return wf.data
    if field.path.startswith('params/Mesh/'):
        return wf.params.Mesh
    return wf.params


def get_glosary_info():
    """
    :return: dictionary field_name -> (description, units)
    """
    return {f.path: (f.description, f.attributes['units']) for f in get_wf_fields()}


def print_glossary():
    """
    Print glossary docs
    """
    for f in get_wf_fields():
        print('**{name}** - {description} - ***{units}***'.format(
            name=f.path, description=f.description, units=f.attributes['units']))


def print_glossary_html():
    """
    Print glossary docs as html table. Used to build alfresco documentation page.
    """
    res = '<table border="2">\n'
    res += '<tr><th>Filed name</th><th>Description</th><th>Units</th></tr>\n'
    for name, descr, units in sorted([f.path, f.description, f.attributes['units']] for f in get_wf_fields()):
        res += '<tr>\n'
        for gg in (name, descr, units):
            res += f'\t<td>{gg}</td>\n'
        res += '</tr>\n'
    res += '</table>\n'
    return res
