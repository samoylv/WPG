# -*- coding: utf-8 -*-
"""
This module contains definitions (glossary) of Wavefront fields. Described mapping fields SRWLWfr <-> wpg.Wavefront

.. module:: wpg.glossary
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import inspect
import sys
import wpg.utils as utils
import warnings
import numpy
import array


class RadiationField(object):

    """
    This is base class for all Wavefront fileds.
    """

    glossary_name = None  # used as path for mapping in wpg.Wavefront

    def __init__(self, wf):
        """
        Used for map values to Wavefront. Also map description string from docstrings and attributes.

        :param wf: Wavefront
        :type wf: wpg.Wavefront
        """

        if not wf.__class__.__name__ == 'Wavefront':
            raise TypeError
        self._wf = wf
        self._value = None
        self.attributes = {'description': self.__doc__,
                           'units': self.find_units_label()}

    def find_units_label(self):
        """Search [units] in field docstring

        :return: units string
        """

        descr = self.__doc__
        descr = descr.replace('\n', ' ')
        units = []
        start = stop = 0
        while start >= 0 and stop >= 0:
            start = descr.find(r'[', stop)
            stop = descr.find(r']', start)
            if start >= 0 and stop >= 0:
                units.append(descr[start+1: stop])

        units = ' or '.join(units)
        return units

    def _map_to_dict(self):
        """
        Map Radiation field to dictionary

        :return: dictionary
        """

        t = {}
        utils.set_value(t, self.keys_chain, self.value)
        return t

    def _map_from_dict(self, dic):
        """Map Radiation field from dictionary
        :param dic: dictionary
        """

        self.value = utils.get_value(dic, self.keys_chain)

    @property
    def value(self):
        """
        Property where value stored.
        """
        raise NotImplemented

    @value.setter
    def value(self, val):
        raise NotImplemented

    @property
    def keys_chain(self):
        """
        Split field name to the parts.

        :return: tpule of strings
        """
        return self.glossary_name.split('/')


class WFVersion(RadiationField):

    """Hdf5 format version (glossary)"""

    glossary_name = 'version'

    def __init__(self, wf):
        """
        Version field.

        :param wf: Wavefront
        :type wf: wpg.Wavefront
        """

        super(WFVersion, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[???]',
             'alias': 'VERSION'
             })

        self._value = 0.1

    @property
    def value(self):
        """Hdf5 format version (glossary)"""

        return self._value

    @value.setter
    def value(self, val):
        self._value = val


class WFRadiationPhotonEnergy(RadiationField):

    """Average photon energy [ev]"""

    glossary_name = 'params/photonEnergy'

    def __init__(self, wf):
        """
        params/photonEnergy field.

        :param wf: Wavefront
        :type wf: wpg.Wavefront
        """

        super(WFRadiationPhotonEnergy, self).__init__(wf)

        self.attributes.update(
            {'limits': '[1:1e6]',
             'alias': 'avgPhotEn'
             })

    @property
    def value(self):
        """Average photon energy [ev]"""

        return self._wf._srwl_wf.avgPhotEn

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.avgPhotEn = float(val)


class WFRadiationMeshZCoord(RadiationField):

    """Longitudinal position [m], Fast data: length of  active undulator,
     Gaussian source: distance to waist """

    glossary_name = 'params/Mesh/zCoord'

    def __init__(self, wf):
        """params/Mesh/zCoord field"""

        super(WFRadiationMeshZCoord, self).__init__(wf)

        self.attributes.update(
            {'limits': '{???}',
             'alias': 'mesh.zStart'
             })

    @property
    def value(self):
        """Longitudinal position, for fast data - length of  active undulator [m]"""

        return self._wf._srwl_wf.mesh.zStart

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.zStart = float(val)


class WFRadiationWDomain(RadiationField):

    """WF in time or frequency (photon energy) domain [string]"""

    glossary_name = 'params/wDomain'

    def __init__(self, wf):
        """params/wDomain field"""

        super(WFRadiationWDomain, self).__init__(wf)

        self.attributes.update(
            {'limits': '{"time", "frequency"}',
             'alias': 'presFT'
             })

    @property
    def value(self):
        """WF in time or frequency (photon energy) domain"""

        if self._wf._srwl_wf.presFT == 0:
            return 'frequency'
        elif self._wf._srwl_wf.presFT == 1:
            return 'time'
        else:
            raise ValueError('internal error, wrong wavefront field vaule')

    @value.setter
    def value(self, val):
        if val in ['frequency', b'frequency']:
            self._wf._srwl_wf.presFT = 0
        elif val in ['time', b'time']:
            self._wf._srwl_wf.presFT = 1
        else:
            raise ValueError('value must be "frequency" or "time"')


class WFRadiationWSpace(RadiationField):

    """Real space or  q-space WF presentation [string]"""

    glossary_name = 'params/wSpace'

    def __init__(self, wf):
        """params/wSpace field"""

        super(WFRadiationWSpace, self).__init__(wf)

        self.attributes.update(
            {'limits': '{"R-space", "Q-space"}',
             'alias': 'presCA'
             })

    @property
    def value(self):
        """Real space or  q-space WF presentation"""

        if self._wf._srwl_wf.presCA == 0:
            return 'R-space'
        elif self._wf._srwl_wf.presCA == 1:
            return 'Q-space'
        else:
            raise ValueError('internal error, wrong wavefront field vaule')

    @value.setter
    def value(self, val):
        if val in ['R-space', b'R-space']:
            self._wf._srwl_wf.presCA = 0
        elif val in ['Q-space', b'Q-space']:
            self._wf._srwl_wf.presCA = 1
        else:
            raise ValueError('value must be "R-space" or "Q-space"')


class WFRadiationWFloatType(RadiationField):

    """Electric field numerical type [string]"""

    glossary_name = 'params/wFloatType'

    def __init__(self, wf):
        """params/wFloatType field"""

        super(WFRadiationWFloatType, self).__init__(wf)

        self.attributes.update(
            {'limits': '{"float", "double"}',
             'alias': 'numTypeElFld'
             })

    @property
    def value(self):
        """Electric field numerical type"""
        if self._wf._srwl_wf.numTypeElFld in ['f', b'f']:
            return 'float'
        elif self._wf._srwl_wf.numTypeElFld in ['d', b'd']:
            return 'double'
        else:
            raise ValueError('internal error, wrong wavefront field value')

    @value.setter
    def value(self, val):
        if val in ['float', b'float']:
            self._wf._srwl_wf.numTypeElFld = 'f'
        elif val in ['double', b'double']:
            raise ValueError('"double" type not supprted yet')
        #            self._wf._srwl_wf.numTypeElFld = 'd'
        else:
            raise ValueError('value must be "float" or "double"')


class WFRadiationWEFieldUnit(RadiationField):

    """Electric field units [string]"""

    glossary_name = 'params/wEFieldUnit'

    def __init__(self, wf):
        """params/wEFieldUnit field"""

        super(WFRadiationWEFieldUnit, self).__init__(wf)

        self.attributes.update(
            {'limits': r'{"sqrt(Phot/s/0.1%bw/mm^2)","sqrt(W/mm^2)",\
             "sqrt(J/eV/mm^2)", "arbitrary"]',
             'alias': 'unitElFld'
             })

    @property
    def value(self):
        """Electric field units"""

        if self._wf._srwl_wf.unitElFld == 0:
            return 'arbitrary'
        elif self._wf._srwl_wf.unitElFld == 1:
            return r'sqrt(Phot/s/0.1%bw/mm^2)'
        elif self._wf._srwl_wf.unitElFld == 2:
            if self._wf._srwl_wf.presFT == 0:
                return "sqrt(J/eV/mm^2)"
            elif self._wf._srwl_wf.presFT == 1:
                return "sqrt(W/mm^2)"
            else:
                raise ValueError('internal error, wrong wavefront field value')

        else:
            raise ValueError('internal error, wrong wavefront field value')

    @value.setter
    def value(self, val):
        if val == 'arbitrary':
            self._wf._srwl_wf.unitElFld = 0
        elif val in ['sqrt(Phot/s/0.1%bw/mm^2)', b'sqrt(Phot/s/0.1%bw/mm^2)']:
            self._wf._srwl_wf.unitElFld = 1
        elif val in ['sqrt(J/eV/mm^2)', 'sqrt(W/mm^2)', b'sqrt(J/eV/mm^2)', b'sqrt(W/mm^2)']:
            self._wf._srwl_wf.unitElFld = 2
        else:
            raise ValueError(
                'value must be "arbitrary" or "sqrt(J/eV/mm^2)"' +
                'or "sqrt(W/mm^2)" or "sqrt(Phot/s/0.1%bw/mm^2)"')


class WFRadiationMeshNx(RadiationField):

    """Numbers of points, horizontal"""

    glossary_name = 'params/Mesh/nx'

    def __init__(self, wf):
        """params/Mesh/nx field"""

        super(WFRadiationMeshNx, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, horizontal"""

        return self._wf._srwl_wf.mesh.nx

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.nx = int(val)


class WFRadiationMeshNy(RadiationField):

    """Numbers of points, vertical"""

    glossary_name = 'params/Mesh/ny'

    def __init__(self, wf):
        """params/Mesh/ny field"""

        super(WFRadiationMeshNy, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.ny

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.ny = int(val)


class WFRadiationMeshNSlices(RadiationField):

    """Numbers of points vs photon energy/time for the pulse"""

    glossary_name = 'params/Mesh/nSlices'

    def __init__(self, wf):
        """params/Mesh/nSlices field"""

        super(WFRadiationMeshNSlices, self).__init__(wf)

        self.attributes.update(
            {'units': '',
             'limits': '[1:LONG_MAX]',
             'alias': 'mesh.ne'
             })

    @property
    def value(self):
        """Numbers of points vs photon energy/time for the pulse"""

        return self._wf._srwl_wf.mesh.ne

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.ne = int(val)
        self._wf._allocate_srw_moments()


class WFRadiationMeshNvx(RadiationField):

    """Lab-frame coordinate of the inner normal to observation plane
    (/ surface in its center)"""

    glossary_name = 'params/Mesh/nvx'

    def __init__(self, wf):
        """params/Mesh/nvx field"""

        super(WFRadiationMeshNvx, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.nvx

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.nvx = int(val)


class WFRadiationMeshNvy(RadiationField):

    """Lab-frame coordinate of the inner normal to observation plane
    (/ surface in its center)"""

    glossary_name = 'params/Mesh/nvy'

    def __init__(self, wf):
        """params/Mesh/nvy field"""

        super(WFRadiationMeshNvy, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.nvy

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.nvy = int(val)


class WFRadiationMeshNvz(RadiationField):

    """Lab-frame coordinate of the inner normal to observation plane
     (/ surface in its center)"""

    glossary_name = 'params/Mesh/nvz'

    def __init__(self, wf):
        """params/Mesh/nvz field"""

        super(WFRadiationMeshNvz, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.nvz

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.nvz = int(val)


class WFRadiationMeshHvx(RadiationField):

    """Lab-frame horizontal base vector of the observation plane
     (/ surface in its center) """

    glossary_name = 'params/Mesh/hvx'

    def __init__(self, wf):
        """params/Mesh/hvx field"""

        super(WFRadiationMeshHvx, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.hvx

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.hvx = int(val)


class WFRadiationMeshHvy(RadiationField):

    """Lab-frame horizontal base vector of the observation plane
     (/ surface in its center) """

    glossary_name = 'params/Mesh/hvy'

    def __init__(self, wf):
        """params/Mesh/hvy field"""

        super(WFRadiationMeshHvy, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.hvy

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.hvy = int(val)


class WFRadiationMeshHvz(RadiationField):

    """Lab-frame horizontal base vector of the observation plane
     (/ surface in its center) """

    glossary_name = 'params/Mesh/hvz'

    def __init__(self, wf):
        """params/Mesh/hvz field"""

        super(WFRadiationMeshHvz, self).__init__(wf)

        self.attributes.update(
            {'units': '-',
             'limits': '[2:LONG_MAX]',
             'alias': ''
             })

    @property
    def value(self):
        """Numbers of points, vertical"""

        return self._wf._srwl_wf.mesh.hvz

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.hvz = int(val)

# TODO: add wrapper for mesh._arSurf (we should know it size)


class WFRadiationMeshXMin(RadiationField):

    """Minimum of horizontal range [m]"""

    glossary_name = 'params/Mesh/xMin'

    def __init__(self, wf):
        """params/Mesh/xMin field"""

        super(WFRadiationMeshXMin, self).__init__(wf)

    @property
    def value(self):
        """Minimum of horizontal range [m]"""

        if self._wf.params.wSpace == 'R-space':
            return self._wf._srwl_wf.mesh.xStart
        else:
            warnings.warn(
                'params/Mesh/xMin not defined if NOT params/wSpace==R-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'R-space':
            self._wf._srwl_wf.mesh.xStart = float(val)
        else:
            warnings.warn(
                'params/Mesh/xMin not defined if NOT params/wSpace==R-space')
            return None


class WFRadiationMeshXMax(RadiationField):

    """Maximum of horizontal range [m]"""

    glossary_name = 'params/Mesh/xMax'

    def __init__(self, wf):
        """params/Mesh/xMax field"""

        super(WFRadiationMeshXMax, self).__init__(wf)

    @property
    def value(self):
        """Maximum of horizontal range [m]"""

        if self._wf.params.wSpace == 'R-space':
            return self._wf._srwl_wf.mesh.xFin
        else:
            warnings.warn(
                'params/Mesh/xMax not defined if NOT params/wSpace==R-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'R-space':
            self._wf._srwl_wf.mesh.xFin = float(val)
        else:
            warnings.warn(
                'params/Mesh/xMax not defined if NOT params/wSpace==Rspace')
            return None


class WFRadiationMeshYMin(RadiationField):

    """Minimum of vertical range [m]"""

    glossary_name = 'params/Mesh/yMin'

    def __init__(self, wf):
        """params/Mesh/yMin field"""

        super(WFRadiationMeshYMin, self).__init__(wf)

    @property
    def value(self):
        """Minimum of vertical range [m]"""

        if self._wf.params.wSpace == 'R-space':
            return self._wf._srwl_wf.mesh.yStart
        else:
            warnings.warn(
                'params/Mesh/yMin not defined if NOT params/wSpace==Rs-pace')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'R-space':
            self._wf._srwl_wf.mesh.yStart = float(val)
        else:
            warnings.warn(
                'params/Mesh/yMin not defined if NOT params/wSpace==R-space')
            return None


class WFRadiationMeshYMax(RadiationField):

    """Maximum of vertical range [m]"""

    glossary_name = 'params/Mesh/yMax'

    def __init__(self, wf):
        """params/Mesh/yMax field"""

        super(WFRadiationMeshYMax, self).__init__(wf)

    @property
    def value(self):
        """Maximum of vertical range [m]"""

        if self._wf.params.wSpace == 'R-space':
            return self._wf._srwl_wf.mesh.yFin
        else:
            warnings.warn(
                'params/Mesh/yMax not defined if NOT params/wSpace==R-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'R-space':
            self._wf._srwl_wf.mesh.yFin = float(val)
        else:
            warnings.warn(
                'params/Mesh/yMax not defined if NOT params/wSpace==R-space')
            return None


class WFRadiationMeshQxMin(RadiationField):

    """Minimum of horizontal frequency [1/m]"""

    glossary_name = 'params/Mesh/qxMin'

    def __init__(self, wf):
        """params/Mesh/qxMin field"""

        super(WFRadiationMeshQxMin, self).__init__(wf)

    @property
    def value(self):
        """Minimum of horizontal frequency [1/m]"""

        if self._wf.params.wSpace == 'Q-space':
            return self._wf._srwl_wf.mesh.xStart
        else:
            warnings.warn(
                'params/Mesh/qxMin not defined if NOT params/wSpace==Q-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'Q-space':
            self._wf._srwl_wf.mesh.xStart = float(val)
        else:
            warnings.warn(
                'params/Mesh/qxMin not defined if NOT params/wSpace==Q-space')
            return None


class WFRadiationMeshQxMax(RadiationField):

    """Maximum of horizontal frequency [1/m]"""

    glossary_name = 'params/Mesh/qxMax'

    def __init__(self, wf):
        """params/Mesh/qxMax field"""

        super(WFRadiationMeshQxMax, self).__init__(wf)

    @property
    def value(self):
        """Maximum of horizontal frequency [1/m]"""

        if self._wf.params.wSpace == 'Q-space':
            return self._wf._srwl_wf.mesh.xFin
        else:
            warnings.warn(
                'params/Mesh/qxMax not defined if NOT params/wSpace==Q-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'Q-space':
            self._wf._srwl_wf.mesh.xFin = float(val)
        else:
            warnings.warn(
                'params/Mesh/qxMax not defined if NOT params/wSpace==Q-space')
            return None


class WFRadiationMeshQyMin(RadiationField):

    """Minimum of vertical frequency [1/m]"""

    glossary_name = 'params/Mesh/qyMin'

    def __init__(self, wf):
        """params/qyMin field"""

        super(WFRadiationMeshQyMin, self).__init__(wf)

    @property
    def value(self):
        """Minimum of vertical frequency [1/m]"""

        if self._wf.params.wSpace == 'Q-space':
            return self._wf._srwl_wf.mesh.yStart
        else:
            warnings.warn(
                'params/Mesh/qyMin not defined if NOT params/wSpace==Q-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'Q-space':
            self._wf._srwl_wf.mesh.yStart = float(val)
        else:
            warnings.warn(
                'params/Mesh/qyMin not defined if NOT params/wSpace==Q-space')
            return None


class WFRadiationMeshQyMax(RadiationField):

    """Maximum of vertical frequency [1/m]"""

    glossary_name = 'params/Mesh/qyMax'

    def __init__(self, wf):
        """params/Mesh/qyMax field"""

        super(WFRadiationMeshQyMax, self).__init__(wf)

    @property
    def value(self):
        """Maximum of vertical frequency [1/m]"""

        if self._wf.params.wSpace == 'Q-space':
            return self._wf._srwl_wf.mesh.yFin
        else:
            warnings.warn(
                'params/Mesh/qyMax not defined if NOT params/wSpace==Q-space')
            return None

    @value.setter
    def value(self, val):
        if self._wf.params.wSpace == 'Q-space':
            self._wf._srwl_wf.mesh.yFin = float(val)
        else:
            warnings.warn(
                'params/Mesh/qyMax not defined if NOT params/wSpace==Q-space')
            return None


class WFRadiationMeshSliceMin(RadiationField):

    """Min value of time [s] or energy [ev] for pulse (fragment)"""

    glossary_name = 'params/Mesh/sliceMin'

    def __init__(self, wf):
        """params/Mesh/sliceMin field"""

        super(WFRadiationMeshSliceMin, self).__init__(wf)

        self.attributes.update(
            {'limits': '[1:LONG_MAX]',
             'alias': 'mesh.eFin'
             })

    @property
    def value(self):
        """Min value of time [s] or energy [ev] for pulse (fragment)"""

        return self._wf._srwl_wf.mesh.eStart

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.eStart = float(val)


class WFRadiationMeshSliceMax(RadiationField):

    """Max value of time [s] or energy [ev] for pulse (fragment)"""

    glossary_name = 'params/Mesh/sliceMax'

    def __init__(self, wf):
        """params/Mesh/sliceMax field"""

        super(WFRadiationMeshSliceMax, self).__init__(wf)

        self.attributes.update(
            {'limits': '[1:LONG_MAX]',
             'alias': 'mesh.eStart'
             })

    @property
    def value(self):
        """Max value of time [s] or energy [ev] for pulse (fragment)"""

        return self._wf._srwl_wf.mesh.eFin

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.mesh.eFin = float(val)


class WFRadiationRx(RadiationField):

    """Instantaneous horizontal wavefront radius [m]"""

    glossary_name = 'params/Rx'

    def __init__(self, wf):
        """params/Rx field"""

        super(WFRadiationRx, self).__init__(wf)
        self.attributes.update(
            {'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'Rx'
             })

    @property
    def value(self):
        """Instantaneous horizontal wavefront radius [m]"""

        return self._wf._srwl_wf.Rx

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.Rx = float(val)


class WFRadiationRy(RadiationField):

    """Instantaneous vertical wavefront radius [m]"""

    glossary_name = 'params/Ry'

    def __init__(self, wf):
        """params/Ry field"""

        super(WFRadiationRy, self).__init__(wf)

        self.attributes.update(
            {'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'Ry'
             })

    @property
    def value(self):
        """Instantaneous vertical wavefront radius [m]"""

        return self._wf._srwl_wf.Ry

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.Ry = float(val)


class WFRadiationDRx(RadiationField):

    """Error of wavefront horizontal radius [m]"""

    glossary_name = 'params/dRx'

    def __init__(self, wf):
        """params/dRx field"""

        super(WFRadiationDRx, self).__init__(wf)
        self.attributes.update(
            {'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'dRx'
             })

    @property
    def value(self):
        """Error of wavefront horizontal radius [m]"""

        return self._wf._srwl_wf.dRx

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.dRx = float(val)


class WFRadiationDRy(RadiationField):

    """Error of wavefront horizontal radius [m]"""

    glossary_name = 'params/dRy'

    def __init__(self, wf):
        """params/dRy field"""

        super(WFRadiationDRy, self).__init__(wf)

        self.attributes.update(
            {'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'dRy'
             })

    @property
    def value(self):
        """Error of wavefront horizontal radius [m]"""

        return self._wf._srwl_wf.dRy

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.dRy = float(val)


class WFRadiationNval(RadiationField):

    """complex electric field nval==2"""

    glossary_name = 'params/nval'

    def __init__(self, wf):
        """params/nval field"""

        super(WFRadiationNval, self).__init__(wf)
        self.attributes.update(
            {'units': '',
             'limits': '{2}',
             'alias': '-'
             })

        self._value = 2

    @property
    def value(self):
        """complex electric field nval==2"""

        return self._value

    @value.setter
    def value(self, val):
        self._value = int(val)


class WFRadiationXCentre(RadiationField):

    """Horizontal transverse coordinates of wavefront instant 'source center' [m]"""

    glossary_name = 'params/xCentre'

    def __init__(self, wf):
        """params/xCentre field"""

        super(WFRadiationXCentre, self).__init__(wf)

    @property
    def value(self):
        """Horizontal transverse coordinates of wavefront instant 'source center' [m]"""

        return self._wf._srwl_wf.xc

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.xc = float(val)


class WFRadiationYCentre(RadiationField):

    """Vertical transverse coordinates of wavefront instant 'source center' [m]"""

    glossary_name = 'params/yCentre'

    def __init__(self, wf):
        """params/yCentre field"""

        super(WFRadiationYCentre, self).__init__(wf)

    @property
    def value(self):
        """Vertical transverse coordinates of wavefront instant 'source center' [m]"""

        return self._wf._srwl_wf.yc

    @value.setter
    def value(self, val):
        self._wf._srwl_wf.yc = float(val)


class WFDataArrEhor(RadiationField):

    """
    EM field (Re, Im) pairs written in 3D array, slice number changes first.
    Horizontal polarization
    """

    glossary_name = 'data/arrEhor'

    def __init__(self, wf):
        """data/arrEhor field"""

        super(WFDataArrEhor, self).__init__(wf)
        self.attributes.update(
            {'units': 'see params/wEFieldUnit',
             'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'arEx'
             })

    @property
    def value(self):
        """
        EM field (Re, Im) pairs written in 3D array, slice number changes first.
        Horizontal polarization
        """

        res = numpy.array(self._wf._srwl_wf.arEx, dtype='float32', copy=False)
        if res.shape:
            res.shape = (self._wf.params.Mesh.ny, self._wf.params.Mesh.nx,
                         self._wf.params.Mesh.nSlices, self._wf.params.nval)
        return res

    @value.setter
    def value(self, val):
        """

        :param val: complex numpy 3D array or array.array. if array.array - just copy
        """
        n_total = self._wf._get_total_elements() * self._wf.params.nval
        self._wf._allocate_srw_moments()
        if type(val) == array.array:
            if not val.count() == n_total:
                warnings.warn(
                    'New array size not equal to wavefront size. You must set it by yourself.')
            self._wf._srwl_wf.arEx = array.array(str(u'f'), val)
        else:
            val = numpy.array(val, dtype='float32')
            if not numpy.prod(val.shape) == n_total:
                warnings.warn(
                    'New array size not equal to wavefront size. It will set automaticaly to array size.')
                self._wf.params.nx = val.shape[1]
                self._wf.params.ny = val.shape[0]
                self._wf.params.nSlices = val.shape[2]
            # self._wf._srwl_wf.arEx = array.array(str(u'f'), val.flat)
            self._wf._srwl_wf.arEx = array.array(str(u'f'))
            val_s = val.tostring()
            self._wf._srwl_wf.arEx.fromstring(val_s)


class WFDataArrEver(RadiationField):

    """
    EM field (Re, Im) pairs written in 3D array, slice number changes first.
    Vertical polarization
    """

    glossary_name = 'data/arrEver'

    def __init__(self, wf):
        """data/arrEver field"""

        super(WFDataArrEver, self).__init__(wf)
        self.attributes.update(
            {'units': 'see params/wEFieldUnit',
             'limits': '[FLOAT_MIN:FLOAT_MAX]',
             'alias': 'arEy'
             })

    @property
    def value(self):
        """
        EM field (Re, Im) pairs written in 3D array, slice number changes first.
        Vertical polarization
        """

        res = numpy.array(self._wf._srwl_wf.arEy, dtype='float32', copy=False)
        if res.shape:
            res.shape = (self._wf.params.Mesh.ny, self._wf.params.Mesh.nx,
                         self._wf.params.Mesh.nSlices, self._wf.params.nval)
        return res

    @value.setter
    def value(self, val):
        n_total = self._wf._get_total_elements() * self._wf.params.nval
        self._wf._allocate_srw_moments()
        if type(val) == array.array:
            if not val.count() == n_total:
                warnings.warn(
                    'New array size not equal to wavefront size. You must set it by yourself.')
            self._wf._srwl_wf.arEy = array.array(str(u'f'), val)
        else:
            val = numpy.array(val, dtype='float32')
            if not numpy.prod(val.shape) == n_total:
                warnings.warn(
                    'New array size not equal to wavefront size. It will set automaticaly to array size.')
                self._wf.params.nx = val.shape[1]
                self._wf.params.ny = val.shape[0]
                self._wf.params.nSlices = val.shape[2]
#             self._wf._srwl_wf.arEy = array.array(str(u'f'), val.flat)
            self._wf._srwl_wf.arEy = array.array(str(u'f'))
            val_s = val.tostring()
            self._wf._srwl_wf.arEy.fromstring(val_s)


# TODO: fix allocation in N(x,y,z)

# TODO: add history section

# TODO: add misc section

def get_wf_fields():
    """
    Return fields in proper order to map it in Wavefront

    :return: iterator over wavefront fields in glossary
    """

    clsmembers = inspect.getmembers(sys.modules[__name__],
                                    inspect.isclass)

    # hack for normalize initialization order

    for (ic, c) in enumerate(clsmembers[:]):
        if c[1].__name__ in ['WFRadiationWSpace', 'WFRadiationMeshNSlices',
                             'WFRadiationMeshNx', 'WFRadiationMeshNy',
                             'WFRadiationNval']:
            ct = clsmembers.pop(ic)
            clsmembers.insert(0, ct)

    for (key, val) in clsmembers:
        if inspect.isclass(val):
            if issubclass(val, RadiationField):
                if not val.__name__ == 'RadiationField':
                    yield val


def get_glosary_info():
    """
    :return: dictionary field_name -> doc
    """

    res = {}
    for wf in get_wf_fields():
        res[wf.glossary_name] = wf.__doc__, wf().find_units_label()
    return res


def print_glossary():
    """
    Print glossary docs
    """

    for wf in get_wf_fields():
        name = wf.glossary_name
        descr = wf.__doc__
        descr = descr.replace('\n', ' ')
        units = []
        start = stop = 0
        while start >= 0 and stop >= 0:
            start = descr.find(r'[', stop)
            stop = descr.find(r']', start)
            if start >= 0 and stop >= 0:
                units.append(descr[start+1: stop])

        units = ' or '.join(units)
        print('**{name}** - {decsription} - ***{units}***'.format(
            name=name, decsription=descr, units=units))
    # pprint.pprint(get_glosary_info())


def print_glossary_html():
    """
    Print glossry docsas html table. Used to build alfresco documentaion page.
    """

    gloss = []
    for wf in get_wf_fields():
        name = wf.glossary_name
        descr = wf.__doc__
        descr = descr.replace('\n', ' ')
        units = []
        start = stop = 0
        while start >= 0 and stop >= 0:
            start = descr.find(r'[', stop)
            stop = descr.find(r']', start)
            if start >= 0 and stop >= 0:
                units.append(descr[start+1: stop])

        units = ' or '.join(units)
        gloss.append([name, descr, units])
    res = r'<table border="2">' + '\n'
    res += r'<tr><th>Filed name</th><th>Description</th><th>Units</th></tr>' + \
        '\n'
    for g in sorted(gloss):
        res += r'<tr>' + '\n'
        for gg in g:
            res += '\t' + r'<td>{}</td>'.format(gg) + '\n'
        res += r'</tr>' + '\n'
    res += r'</table>' + '\n'
    return res

# class GlossaryField(dict):

#     def __init__(
#         self,
#         type=None,
#         units=None,
#         limits=None,
#         alias=None,
#         description=None,
#         ):
#         self['type'] = type
#         self['units'] = units
#         self['limits'] = limits
#         self['alias'] = alias
#         self['description'] = description


# class GlossaryDict(object):

#     def __init__(self):
#         self.definition = {
#             'VERSION': GlossaryField(type='float',
#                     description='Hdf5 format version (glossary)'),
#             'history': {
#                 'parentContact': GlossaryField(type='string',
#                         description=r'For fast data: section Contacts from readme'
#                         ),
#                 'parentData': GlossaryField(type='string',
#                         description=r'Source file name'),
#                 'parentMethods': GlossaryField(type='string',
#                         description=r'For fast data - first 13 lines + "Conversion from polar to xy: fast2xy" from readme.txt'
#                         ),
#                 'parentTestFiles': GlossaryField(type='string',
#                         description=r'For fast data – section test files of readme'
#                         ),
#                 },
#             'Radiation': {
#                 'photonEnergy': GlossaryField(type='float',
#                         alias='avgPhotEn', units='ev',
#                         description=r'average photon energy'),
#                 'zCoord': GlossaryField(type='double',
#                         alias='mesh.zStart', units='m',
#                         description=r'Longitudinal position, for fast data - length of  active undulator'
#                         ),
#                 'wDomain': GlossaryField(type='string', limits=('twEFieldUnitime',
#                         'frequency'), alias='presFT',
#                         description=r'WF in time || frequency domain '
#                         ),
#                 'wSpace': GlossaryField(type='string', limits=('Rspace', 'Qspace'), alias='presCA',
#                         description=r'real space || q-space WF presentation'
#                         ),
#                 'wFloatType': GlossaryField(type='string',
#                         limits=('Float', 'Double'), alias='numTypeElFld', description=r'electric field numerical type'
#                         ),
#                 'wEFieldUnit ': GlossaryField(type='string',
#                         limits=('sqrt(GW)',
#                         'sqrt(Nphoton/s/0.1%bw/mm^2)', 'sqrt(Nphotons)'
#                         ), alias='unitElFld',
#                         description=r'units of electric field intensity |E|^2'
#                         ),
#                 'nx': GlossaryField(type='long', alias='mesh.nx',
#                                     description=r'numbers of points, horizontal'
#                                     ),
#                 'ny': GlossaryField(type='long', alias='mesh.ny',
#                                     description=r'numbers of points, vertical'
#                                     ),
#                 'xMin': GlossaryField(type='double', units='m',
#                         alias='mesh.xStart',
#                         description=r'minimum of horizontal range'),
#                 'qxMin': GlossaryField(type='double', units='rad',
#                         alias='mesh.xStart',
#                         description=r'minimum of horizontal range'),
#                 'xMax': GlossaryField(type='double', units='m',
#                         alias='mesh.xFin',
#                         description=r'maximum of horizontal range'),
#                 'qxMax': GlossaryField(type='double', units='rad',
#                         alias='mesh.xFin',
#                         description=r'maximum of horizontal range'),
#                 'yMin': GlossaryField(type='double', units='m',
#                         alias='mesh.yStart',
#                         description=r'minimum of vertical range'),
#                 'qyMin': GlossaryField(type='double', units='rad',
#                         alias='mesh.yStart',
#                         description=r'minimum of vertical range'),
#                 'yMax': GlossaryField(type='double', units='m',
#                         alias='mesh.yFin',
#                         description=r'maximum of vertical range'),
#                 'qyMax': GlossaryField(type='double', units='rad',
#                         alias='mesh.yFin',
#                         description=r'maximum of vertical range'),
#                 'nSlices': GlossaryField(type='long', alias='mesh.ne',
#                         description=r'numbers of points vs photon energy/time for the pulse'
#                         ),
#                 'sliceMin': GlossaryField(type='double',
#                         alias='mesh.eStart',
#                         description=r'min value of time [s] or energy [ev] for pulse [fragment]'
#                         ),
#                 'sliceMax': GlossaryField(type='double',
#                         alias='mesh.eFin',
#                         description=r'max value of time [s] or energy [ev] for pulse [fragment]'
#                         ),
#                 'Rx': GlossaryField(type='double', units='m',
#                                     alias='Rx,',
#                                     description=r'instantaneous horizontal wavefront radius'
#                                     ),
#                 'Ry': GlossaryField(type='double', units='m', alias='Ry',
#                                     description=r'instantaneous vertical wavefront radius'
#                                     ),
#                 'dRx': GlossaryField(type='double', units='m',
#                         alias='dRx',
#                         description=r'error of wavefront horizontal radius'
#                         ),
#                 'dRy': GlossaryField(type='double', units='m',
#                         alias='dRy',
#                         description=r'error of wavefront vertical radius'
#                         ),
#                 'nval': GlossaryField(type='int',
#                         description='complex electric field nval==2'),
#                 'xCentre': GlossaryField(type='double', units='m',
#                         alias='xC',
#                         description=r'horizontal transverse coordinates of wavefront instant "source center"'
#                         ),
#                 'yCentre': GlossaryField(type='double', units='m',
#                         alias='yC',
#                         description=r'vertical transverse coordinates of wavefront instant "source center"'
#                         ),
#                 },
#             'Data': {'arrEhor': GlossaryField(type='float array',
#                      alias='arEx',
#                      description='EM field (Re, Im) pairs written in 3D array, slice number changes first. Horizontal polarization '
#                      ), 'arrEver': GlossaryField(type='float array',
#                      alias='arEy',
#                      description=r'EM field. Vertical polarization')},
#             'misc': {
#                 'WignerDistrMoments': GlossaryField(type='string',
#                         alias='arMomX,arMomY',
#                         description=r'h5 structure with statistical moments (of Wigner distribution)'
#                         ),
#                 'aux': GlossaryField(type='string', alias='arWfrAuxData',
#                         description=r'Name of file with array of auxiliary WF data'
#                         ),
#                 'SRWLPartBeam': GlossaryField(type='string',
#                         alias='partBeam',
#                         description=r'filename with source parameters (e.g. e-beam for SR or e-bunch for FEL ?)'
#                         ),
#                 'ElecPropMatr': GlossaryField(type='string',
#                         alias='arElecPropMatr', description=r''),
#                 },
#             }
