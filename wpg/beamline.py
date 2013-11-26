#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__ = 'makov'

import srwlib
from utils import srw_obj2str


class Beamline(object):
    def __init__(self, srwl_beamline=None):
        if srwl_beamline is None:
            self.srwl_beamline = srwlib.SRWLOptC()
        else:
            self.srwl_beamline = srwl_beamline

    def __str__(self):
        res = ''
        for elem, pp in zip(self.srwl_beamline.arOpt, self.srwl_beamline.arProp):
            s1 = elem.__doc__
            s2 = 'Prop. parameters = {0}'.format(pp)
            s3 = '\t'+'\n\t'.join(srw_obj2str(elem).split('\n'))
            res += '{0}\n{1}\n{2}\n'.format(s1, s2, s3)
        return res


def flatten(lst):
    for x in lst:
        if isinstance(x, list):
            for y in flatten(x):
                yield y
        else:
            yield x
