# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

"""
This module contains  wrapper for SRWLOptC (optical container) and propagation parameters.

.. module:: wpg.beamline
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""

import wpg.srwlib as srwlib
from wpg.srwlib import srwl
from wpg.utils import srw_obj2str
import wpg.optical_elements


class Beamline(object):
    """
    Set of optical elements and propagation parameters.
    """

    def __init__(self, srwl_beamline=None):
        """
        Init beamline.

        :params srwl_beamline: if present will used for initialization.
        :type srwl_wavefront: SRWLOptC
        """
        self.propagation_options = [{'optical_elements': [],
                                    'propagation_parameters':[]}]
        if srwl_beamline is not None:
            tolal_elements = max(
                len(srwl_beamline.arProp), len(srwl_beamline.arOpt))
            for ti in range(tolal_elements):
                try:
                    elem = srwl_beamline.arOpt[ti]
                except IndexError:
                    elem = wpg.optical_elements.Empty()

                try:
                    pp = srwl_beamline.arProp[ti]
                except IndexError:
                    pp = None
                self.append(elem, pp)

    def __str__(self):
        """
        String representaion of beamline (used with print function).

        :return: string
        """
        res = ''
        for po in self.propagation_options:
            tolal_elements = max(
                len(po['optical_elements']), len(po['propagation_parameters']))
            for ti in range(tolal_elements):
                try:
                    elem = po['optical_elements'][ti]
                except IndexError:
                    elem = wpg.optical_elements.Empty()

                try:
                    pp = po['propagation_parameters'][ti]
                except IndexError:
                    pp = None

                s1 = elem.__doc__
                s2 = 'Prop. parameters = {0}'.format(pp)
                if isinstance(elem, srwlib.SRWLOpt):
                    s3 = '\t' + '\n\t'.join(srw_obj2str(elem).split('\n'))
                else:
                    s3 = '\t' + str(elem)

                res += '{0}\n{1}\n{2}\n'.format(s1, s2, s3)
        return res

    def append(self, optical_element, propagation_parameters):
        """
        Appends optical element and propagation propagation parameters to the end of beamline

        :param optical_element: SRW or wpg optical element
        :param propagation_parameters: SRW propagation parameters list or wpg.optical_elements.UsePP object
        """

        # TODO: check types
        last_pp_opt = self.propagation_options[-1]
        # if numbers of propagation parameters and optical elements different
        # create new propagation option
        if not len(last_pp_opt['optical_elements']) == len(last_pp_opt['propagation_parameters']):
            self.propagation_options.append({'optical_elements': [],
                                             'propagation_parameters': []})
            last_pp_opt = self.propagation_options[-1]

        if not all([isinstance(o, srwlib.SRWLOpt) for o in last_pp_opt['optical_elements']]):
            self.propagation_options.append({'optical_elements': [],
                                             'propagation_parameters': []})

        # if current parameter is SRWOpt and last parameter was SRWOpt lets
        # stack it
        if isinstance(optical_element, srwlib.SRWLOpt):
            opt = optical_element
            pp = _get_srw_pp(propagation_parameters)

            last_pp_opt['optical_elements'].append(opt)
            last_pp_opt['propagation_parameters'].append(pp)

        # support resizing element
        if optical_element == [] or isinstance(optical_element, wpg.optical_elements.Empty):
            pp = _get_srw_pp(propagation_parameters)

            last_pp_opt['propagation_parameters'].append(pp)

        # self.srwl_beamline.arOpt.append(optical_element)
        # self.srwl_beamline.arProp.append(propagation_parameters)
    def propagate(self, wfr):
        """
        Propagate wavefront through beamline.

        :param wfr: Input wavefront (will be re-writed after propagation)
        :type wfr: wpg.wavefront.Wavefront
        """
        for propagation_option in self.propagation_options:
            if all([isinstance(o, srwlib.SRWLOpt) for o in propagation_option['optical_elements']]):
                srwl_beamline = srwlib.SRWLOptC(
                    propagation_option['optical_elements'],
                    propagation_option['propagation_parameters'])
                
                srwl.PropagElecField(wfr._srwl_wf, srwl_beamline)

                # fixing wf._srwl_wf.mesh.zStart bug for Drift
                # TODO: Try to fix it with _treat parameter
                for opt_element in propagation_option['optical_elements']:
                    if isinstance(opt_element, srwlib.SRWLOptD):
                        wfr.params.Mesh.zCoord = wfr.params.Mesh.zCoord + \
                            opt_element.L
            else:
                raise ValueError('Unknown type of propagators')


def _check_srw_pp(pp):
    """
    Check is propagation parameters valid SRW propagation parameters

    :param pp: propagation parameters
    :type pp: list of floats
    """
    return isinstance(pp, list) and len(pp) in [12, 17]


def _get_srw_pp(propagation_parameters):
    """ Try to get propagation parameters from object calling get_srw_pp() method

    :param propagation_parameters: propagation parameters
    :type propagation_parameters: list of floats
    :return: propagation SRW parameters
    """
    if _check_srw_pp(propagation_parameters):
        return propagation_parameters
    elif 'get_srw_pp' in dir(propagation_parameters):
        tmp_pp = propagation_parameters.get_srw_pp()
        if _check_srw_pp(tmp_pp):
            return tmp_pp

    raise TypeError(
        'Propagation parameters should be a list of 12 or 17 numbers or have valid "get_srw_pp" method')
