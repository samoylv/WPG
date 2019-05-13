# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import sys
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/utils")
sys.path.append(r"C:\Users\twguest\OneDrive - LA TROBE UNIVERSITY\code\WPG_Wavefront_Simulations/root")

import copy
import numpy
import pylab

from wpg.beamline import Beamline
from wpg.srwlib import *
from wpg.wavefront import Wavefront
from wpg.wpg_uti_wf import *
from wpg.useful_code.wfrutils import plot_wfront
from bl_utils import wfr_center
from wpg.optical_elements import Empty# -*- coding: utf-8 -*-
################################# TO BE MOVED
def plot(wfr):
    plot_wfront(wfr, title_fig='at '+ str(14) +' m',
                isHlog=False, isVlog=False,
                i_x_min=0, i_y_min=0, orient='x', onePlot=False)
##############################################
def zoom_param(zoom_x, zoom_y):
    ppzoom = [0, 0, 1.0, 0, 0, 1/zoom_x, zoom_x, 1/zoom_y, zoom_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return ppzoom

def zoom_wfr(wfr, zoom_x, zoom_y = None, show = True):
    """
    Propagate wavefront through beamline.

    :param wfr: Input wavefront (will be re-writed after propagation)
    :type wfr: wpg.wavefront.Wavefront
    """
    if zoom_y == None:
        zoom_y = zoom_x

    if isinstance(wfr, Wavefront):
        wfr = Wavefront(srwl_wavefront=wfr._srwl_wf)
    

        
    bl = Beamline("Temporary Beamline")
    bl.append(SRWLOptD(0), zoom_param(zoom_x, zoom_y))
    bl.propagate(wfr)   
    
    if show == True:
        plot(wfr)
    del bl
    
def rescale_param(scalefactor_x, scalefactor_y = None):
    if scalefactor_y == None:
        scalefactor_y = scalefactor_x
    ppzoom = [0, 0, 1.0, 0, 0, 1, scalefactor_x, 1, scalefactor_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    return ppzoom

def rescale_wfr(wfr, rescale_x, rescale_y):
    """
    Propagate wavefront through beamline.

    :param wfr: Input wavefront (will be re-writed after propagation)
    :type wfr: wpg.wavefront.Wavefront
    """

    if isinstance(wfr, Wavefront):
        wfr = Wavefront(srwl_wavefront=wfr._srwl_wf)
    

        
    bl = Beamline("Temporary Beamline")
    bl.append(SRWLOptD(0), rescale_param(rescale_x, rescale_y))
    bl.propagate(wfr)   
    
    del bl