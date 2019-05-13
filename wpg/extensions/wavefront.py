# -*- coding: utf-8 -*-
"""
Created on Tue May 14 09:17:08 2019

@author: twguest
"""

import array
import warnings

import numpy as np
import h5py
import pylab as plt
import wpg.srwlib as srwlib

import wpg.utils as utils
import wpg.glossary as glossary

from wpg.utils import srw_obj2str
import imageio
from wpg.beamline import *

from wpg.srwlib import srwl,SRWLOptD,SRWLOptA,SRWLOptC,SRWLOptT,SRWLOptL,SRWLOptMirEl

from coherence_1d import calc_1D_coherent_fraction as calc_coherence
from coherence_1d import plot_1D_degree_of_coherence as plot_coherence
from coherence_1d import calc_degree_of_transverse_coherence_PCA as TDOC

from root.wpg.wavefront import Wavefront as SRWLWavefront

class Wavefront(SRWLWavefront):
    
    def __init__(self, srwl_wavefront = None):
        super().__init__(srwl_wavefront)

   
    def rescale(self, scale_x, scale_y):
        """
    
        :param scale_x: pixel magnification along x-axis
        :param scale_y: pixel magnification along y_axis
        """
    
        
        def rescale_param(scale_x, scale_y = None):
            
            if scale_y == None:
                scale_y = scale_x
            ppscale = [0, 0, 1.0, 0, 0, 1, scale_x, 1, scale_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            return ppscale
    
    
            
        bl = Beamline("Temporary Beamline")
        bl.append(SRWLOptD(0), rescale_param(scale_x, scale_y))
        bl.propagate(self)   
        
        del bl
        

    
    def zoom(self, zoom_x, zoom_y = None):
        """
        :param scale_x: magnification along x-axis
        :param scale_y: magnification along y_axis
        """
        if zoom_y == None:
            zoom_y = zoom_x
    
        def zoom_param(zoom_x, zoom_y):
            ppzoom = [0, 0, 1.0, 0, 0, 1/zoom_x, zoom_x, 1/zoom_y, zoom_y, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            return ppzoom
    
        bl = Beamline("Temporary Beamline")
        bl.append(SRWLOptD(0), zoom_param(zoom_x, zoom_y))
        bl.propagate(self)   
        
        del bl
        
    def save_tif(self, filename):
        imgReal = self.get_intensity(slice_number=0,polarization='horizontal')
        imageio.imwrite(filename + ".tif", imgReal)
        
    def get_amp(self):
        
        amp = self.data.arrEhor[:,:,:,1]

        return amp
    
    def coherence(self, axisName = 'x', plot = True):
        axis = np.linspace(self.params.Mesh.xMin, self.params.Mesh.xMax, self.params.Mesh.nx)
        amp = self.get_amp()

        coh = calc_coherence(amp, axisName, axis)
        tdoc = TDOC(amp)
        
        
        print("Transverse Degree of Coherence: {}".format(tdoc))
        if plot == True:
            plot_coherence(coh, axisName, axis)
        
        return coh, tdoc
        
    def plot(self):
        

        try: 
            plt.figure()
            plt.imshow(self.II[:,:,0])
        except(NameError):
            print("Single or Non-Collapsible Wavefront")
        
    def pixelsize(self):
        
        px = (self.params.Mesh.xMax-self.params.Mesh.xMin)/self.params.Mesh.nx
        py = (self.params.Mesh.yMax-self.params.Mesh.yMin)/self.params.Mesh.ny
        
        return px, py
    

    
