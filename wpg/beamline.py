# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov'

import wpg.srwlib as srwlib
from wpg.srwlib import srwl
from wpg.utils import srw_obj2str
import wpg.optical_elements


class Beamline(object):

    def __init__(self):
        self.propagation_options = [{'optical_elements': [],
                                    'propagation_parameters':[]}]

    def __str__(self):
        # TODO: fix beamline printing
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
            pp = get_srw_pp(propagation_parameters)

            last_pp_opt['optical_elements'].append(opt)
            last_pp_opt['propagation_parameters'].append(pp)

        # support resizing element
        if optical_element == [] or isinstance(optical_element, wpg.optical_elements.Empty):
            pp = get_srw_pp(propagation_parameters)

            last_pp_opt['propagation_parameters'].append(pp)

        # self.srwl_beamline.arOpt.append(optical_element)
        # self.srwl_beamline.arProp.append(propagation_parameters)
    def propagate(self, wfr):
        for propagation_option in self.propagation_options:
            if all([isinstance(o, srwlib.SRWLOpt) for o in propagation_option['optical_elements']]):
                srwl_beamline = srwlib.SRWLOptC(
                    propagation_option['optical_elements'],
                    propagation_option['propagation_parameters'])

                srwl.PropagElecField(wfr._srwl_wf, srwl_beamline)

                # fixing wf._srwl_wf.mesh.zStart bug for Drift
                for opt_element in propagation_option['optical_elements']:
                    if isinstance(opt_element, srwlib.SRWLOptD):
                        wfr.params.Mesh.zCoord = wfr.params.Mesh.zCoord + \
                            opt_element.L
            else:
                raise ValueError('Unknown type of propagators')


def check_srw_pp(pp):
    """ Check is propagation parameters valid SRW propagation parameters"""
    return isinstance(pp, list) and len(pp) == 12


def get_srw_pp(propagation_parameters):
    """ Try to get propagation parameters from object call in get_srw_pp() method"""
    if check_srw_pp(propagation_parameters):
        return propagation_parameters
    elif 'get_srw_pp' in dir(propagation_parameters):
        tmp_pp = propagation_parameters.get_srw_pp()
        if check_srw_pp(tmp_pp):
            return tmp_pp

    raise TypeError(
        'Propagation parameters should be a list of 12 numbers or have valid "get_srw_pp" method')
