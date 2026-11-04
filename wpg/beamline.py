# -*- coding: utf-8 -*-

"""
This module contains  wrapper for SRWLOptC (optical container) and propagation parameters.

.. module:: wpg.beamline
   :platform: Linux, Mac OSX, Windows

.. moduleauthor:: Alexey Buzmakov <buzmakov@gmail.com>
"""
import warnings

from srwpy import srwlib
from srwpy.srwlib import srwl
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
                    elem = po.get('screen') or wpg.optical_elements.Empty()

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

        # if current parameter is SRWOpt and last parameter was SRWOpt lets
        # stack it
        if isinstance(optical_element, srwlib.SRWLOpt):
            opt = optical_element
            pp = _get_srw_pp(propagation_parameters)
            _drift_length(opt)  # warns about interferometer branches of different length

            last_pp_opt['optical_elements'].append(opt)
            last_pp_opt['propagation_parameters'].append(pp)

        # support resizing element: its propagation parameters close the group
        # (SRW applies them as resizing after the last element), Screen also stores the wavefront
        elif isinstance(optical_element, wpg.optical_elements.Empty) or \
                (isinstance(optical_element, list) and not optical_element):
            pp = _get_srw_pp(propagation_parameters)

            last_pp_opt['propagation_parameters'].append(pp)
            if isinstance(optical_element, wpg.optical_elements.Screen):
                last_pp_opt['screen'] = optical_element

        else:
            raise TypeError(
                'Optical element should be an SRW optical element (SRWLOpt), wpg Empty or Screen, '
                'received {}'.format(type(optical_element)))
    def propagate(self, wfr, gpu=None):
        """
        Propagate wavefront through beamline.

        :param wfr: Input wavefront (will be re-writed after propagation)
        :type wfr: wpg.wavefront.Wavefront
        :param gpu: number of the GPU device (1 - first GPU) to propagate on,
            None - propagate on CPU. Requires srwpy built with CUDA
            (``make srw SRW_MODE=cuda``), otherwise SRW propagates on CPU.
            Experimental: SRW 4.2.1 on GPU agrees with CPU for a single slice and,
            for several slices, with semi-analytical treatment of the quadratic phase
            (``Use_PP(semi_analytical_treatment=1)``); for several slices without it
            the drift length is ignored (see tests/test_wpg.py::test_gpu_option).
        :type gpu: int or None
        """
        if wfr.params.wDomain == 'time':
            # SRW propagates every slice as a photon energy: in the time domain the result is noise
            # (5 keV gaussian, 20 slices, drift 5 m: RMS size 875 um instead of 452 um)
            warnings.warn('Wavefront is in the time domain: SRW propagates only in the frequency '
                          "domain, call srwl.SetRepresElecField(wfr._srwl_wf, 'f') before propagation")
        for propagation_option in self.propagation_options:
            if all([isinstance(o, srwlib.SRWLOpt) for o in propagation_option['optical_elements']]):
                optical_elements = propagation_option['optical_elements']
                propagation_parameters = propagation_option['propagation_parameters']
                if not propagation_parameters:
                    continue  # nothing to do (empty beamline or all elements of the group removed)
                if not optical_elements:
                    # Empty as the first element of the group: SRW 4.2.1 rejects SRWLOptC([], [pp])
                    # ("No objects were submitted for parsing"), resize with a drift of zero length
                    optical_elements = [srwlib.SRWLOptD(0.)]
                srwl_beamline = srwlib.SRWLOptC(optical_elements, propagation_parameters)

                if gpu is None:
                    srwl.PropagElecField(wfr._srwl_wf, srwl_beamline)
                else:
                    srwl.PropagElecField(wfr._srwl_wf, srwl_beamline, None, int(gpu))

                # fixing wf._srwl_wf.mesh.zStart bug for Drift (also inside nested containers)
                # TODO: Try to fix it with _treat parameter
                for opt_element in propagation_option['optical_elements']:
                    wfr.params.Mesh.zCoord = wfr.params.Mesh.zCoord + \
                        _drift_length(opt_element, warn=False)

                screen = propagation_option.get('screen')
                if screen is not None:
                    screen.store(wfr)
            else:
                raise ValueError('Unknown type of propagators')


    def _elements(self):
        """(group, number in the group, element) for optical elements of all propagation options."""
        for group, option in enumerate(self.propagation_options):
            for number, element in enumerate(option['optical_elements']):
                yield group, number, element

    @property
    def index(self):
        """
        Names of optical elements and their numbers along the whole beamline (all propagation
        options, elements without a name are 'Element <number>').

        :return: dict {name: number}
        """
        index = {}
        for itr, (_, _, a) in enumerate(self._elements()):
            try:
                index[a.name] = itr
            except AttributeError:
                index['Element {}'.format(itr)] = itr
        return index

    def _locate(self, oe_name):
        """
        (group, number in the group) of the optical element named oe_name.

        :raises KeyError: if there is no element with this name
        """
        positions = [(group, number) for group, number, _ in self._elements()]
        return positions[self.index[oe_name]]

    def get_element_properties(self, oe_name):
        group, number = self._locate(oe_name)
        return self.propagation_options[group]['optical_elements'][number].__dict__

    def get_propagation_paramters(self, oe_name):
        group, number = self._locate(oe_name)
        print(self.propagation_options[group]['propagation_parameters'][number])

    def replace_element(self, oe_name, el):
        """
        replace an optical element labelled oe_name w/ a pre-defined element oe.
        """
        group, number = self._locate(oe_name)
        self.propagation_options[group]['optical_elements'][number] = el

    def remove_element(self, oe_name):
        """
        remove an optical element from a beamline
        """
        group, number = self._locate(oe_name)
        del self.propagation_options[group]['propagation_parameters'][number]
        del self.propagation_options[group]['optical_elements'][number]

    def edit_element_property(self, oe_name, prop, value):
        group, number = self._locate(oe_name)
        self.propagation_options[group]['optical_elements'][number].__dict__[prop] = value

    def edit_propagation_parameters(self, oe_name, new_parameters):
        """
        edit the propagation parameters of an element by name
        """
        group, number = self._locate(oe_name)
        self.propagation_options[group]['propagation_parameters'][number] = _get_srw_pp(new_parameters)

    def crop_beamline(self, element1=None, element2=None):
        """
        crop the beamline after the optical element named element1 (element1 is kept)

        Elements and propagation parameters after element1 are removed, including the resizing
        (Empty) after it and all following propagation options.

        :param element1: name of the last optical element to keep, None or unknown name - do nothing
        :param element2: not used
        """
        if element1 is None or element1 not in self.index:
            return
        group, number = self._locate(element1)
        option = self.propagation_options[group]
        option['optical_elements'] = option['optical_elements'][:number + 1]
        option['propagation_parameters'] = option['propagation_parameters'][:number + 1]
        option.pop('screen', None)  # the resizing of the Screen after element1 is removed
        del self.propagation_options[group + 1:]


def _drift_length(element, warn=True):
    """
    Total length of drifts in an SRW optical element, including nested containers.

    SRW does not change mesh.zStart at a drift, Beamline adds this length to zCoord.
    For an interferometer (SRWLOptI) the length of the branch whose mesh SRW takes at the
    recombination (arPar[0]) is used.

    :param warn: warn if branches of an interferometer have different drift lengths
    """
    if isinstance(element, srwlib.SRWLOptD):
        return element.L
    if isinstance(element, srwlib.SRWLOptC):
        return sum(_drift_length(e, warn) for e in element.arOpt)
    opt_i = getattr(srwlib, 'SRWLOptI', None)  # interferometer, SRW after 4.2.1
    if opt_i is not None and isinstance(element, opt_i) and element.arOptC:
        lengths = [_drift_length(c, warn) for c in element.arOptC]
        if warn and max(lengths) - min(lengths) > 1e-9 * max(1., max(abs(v) for v in lengths)):
            warnings.warn('Drift lengths of interferometer branches differ ({}), zCoord is '
                          'shifted by the length of branch {}'.format(lengths, element.arPar[0]))
        return lengths[element.arPar[0]]
    return 0.


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
