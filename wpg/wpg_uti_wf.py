# -*- coding: utf-8 -*-
__author__ = 'A. Buzmakov, L. Samoylova'

import time
import numpy
#import pylab

## Import standart libraries and addnig "../wavefront" directory to python
## search path
#import os
#import sys
#sys.path.insert(0, os.path.join('..','..'))

from wpg.wavefront import Wavefront
from wpg.beamline  import  Beamline
# from srwlib import *
import wpg.srwlib
from   wpg.srwlib import srwl

def print_mesh(wfr):
    """
    Print out wfr wavefront mesh.
    """    
 
    wf_mesh = wfr.params.Mesh
    print 'nx {:5d}  range_x [{:.1e}, {:.1e}]'.format(wf_mesh.nx,wf_mesh.xMin,wf_mesh.xMax)
    print 'ny {:5d}  range_y [{:.1e}, {:.1e}]'.format(wf_mesh.ny,wf_mesh.yMin,wf_mesh.yMax)   


def propagate_wavefront(wavefront, beamline, output_file = None):
    """
    Propagate wavefront and store it in output file.
    
    :param wavefront: Wavefront object or path to HDF5 file
    :param beamline: SRWLOptC container of beamline
    :param output_file: if parameter present - store propagaed wavefront to file
    :return: propagated wavefront object:
    """
    
    if not isinstance(beamline, Beamline):
        bl = Beamline(beamline)
    else:
        bl = beamline
    print bl
    if isinstance(wavefront, Wavefront):
        wfr = Wavefront(srwl_wavefront=wavefront._srw_wf)
    else:
        print '*****reading wavefront from h5 file...'
        wfr = Wavefront()
        wfr.load_hdf5(wavefront)
        
    
    print_mesh(wfr)
    print '*****propagating wavefront (with resizing)...'
    bl.propagate(wfr)

    if not output_file is None:
        print 'save hdf5:', output_file
        wfr.store_hdf5(output_file)
    print 'done'
    return wfr