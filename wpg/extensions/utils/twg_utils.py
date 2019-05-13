# -*- coding: utf-8 -*-
from wpg import Wavefront
from matplotlib import pyplot as plt
from math import *
from wpg.srwlib import srwl,SRWLOptD,SRWLOptA,SRWLOptC,SRWLOptT,SRWLOptL,SRWLOptMirEl
import numpy as np
import os
from datetime import datetime
import sys

def h5_name(name, itr, strOutput):
    fname0 = 'twg' + name + str(itr)
    print('*****saving hdf5: '+ fname0 + '.h5')
    file_name = os.path.join(strOutput,fname0+'.h5')
    itr+=1
    return file_name

def plot_intensity(wavefront):
    """IN 2D
    PLOT THE INTENSITY OF A WAVFRONT"""
    print("******plotting wavefront intensity")
    
    # calculate intensity
    wavefront_intensity = Wavefront.get_intensity(wavefront)
    wavefront_intensity = wavefront_intensity.sum(-1)
    
    plt.figure()
    plt.set_cmap('jet')
    plt.imshow(wavefront_intensity)
    
    # remove ticks
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])

    # color bar
    cbar = plt.colorbar()
    cbar.ax.set_ylabel('Intensity')
    



def wfr_dose(wavefront):
    """ IN 2D
    CALCULATE THE DOSE AND DOSE DENSITY OF A WAVEFRONT IN W AND W/ MM^2"""
    
    intensity = wavefront.get_intensity(slice_number = 0, polarization = 'horizontal')
    pE = wavefront.params.photonEnergy 
    J2EV = 6.24150934e18
    
    Power = intensity*wavefront.params.photonEnergy/J2EV # 2D power mesh across face of wfr # GW
    
    MeshDim = [xMin, xMax, yMin, yMax] = wavefront.get_limits()
    dx = MeshDim[1]-MeshDim[0]
    dy = MeshDim[3]-MeshDim[2]
    
    mwf_area = abs(dx*dy*(1e3)**2) # mm^2
    npoints = np.shape(intensity)[0]*np.shape(intensity)[1]
    apoint = mwf_area/npoints

    plt.figure()
    plt.imshow(Power, extent=[xMin*1e3, xMax*1e3, yMin*1e3, yMax*1e3])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Dose [W/$mm^2$]")
    plt.title("Dose Density")
    plt.ylabel("y [mm]")
    plt.xlabel("x [mm]")
    
    plt.figure()
    plt.imshow(Power*apoint, extent=[xMin*1e3, xMax*1e3, yMin*1e3, yMax*1e3])
    cbar = plt.colorbar()
    cbar.ax.set_ylabel("Dose [W/$mm^2$]")
    plt.title("Dose")
    plt.ylabel("y [mm]")
    plt.xlabel("x [mm]")
    
def simpleplot(mwf):
    
    print("*****Plotting Intensity and Phase at {} m".format(mwf.params.Mesh.zCoord))
    
    plt.subplot(1,2,1)
    plt.title("Intensity")
    plt.imshow(mwf.get_intensity(slice_number=0))
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])
    
    plt.subplot(1,2,2)
    plt.title("Phase")
    plt.imshow(mwf.get_phase(slice_number=0,polarization='horizontal'))
    plt.gca().set_xticks([])
    plt.gca().set_yticks([])

import math
import os
import scipy
try:
    from wpg import srwlpy as srwl
except ImportError:
    import srwlpy as srwl  #  Hack for read the docs


from wpg.uti_io import *
from wpg.srwlib import SRWLOptT
from PIL import Image


import numpy as np
#%%
class log:
    
    def __init__(self, time):
        self.time = time
    def start(self):
        self.f = open("log/log" + self.time + ".txt", 'w')
        sys.stdout = self.f
    def stop(self):
        self.f.close()

class Experiment:

    def __init__(self, islog = True):
        self.time = datetime.now().strftime('_%Y-%m-%d_%H_%M_%S')
        self.info = {"Experiment Start": self.time}
        self.islog = islog
        self.L = log(self.time)
   
    def start(self):
        if self.islog == True:
            self.L.start()
            self.info["log"] = "out/log" + self.time + ".txt"
        elif self.islog == False:
            self.info["log"] = "False"
            
    def stop(self):
        "Stop Logging"
    
    def info(self):
        return self.info 

    
def plot_treyfront(mwf, title_fig, i_x_min, i_y_min, width, scale, X1, X2, Fineplot = None):
    from wpg.useful_code.wfrutils import calculate_peak_pos, get_mesh, calculate_fwhm_x, calculate_fwhm_y
    import pylab
    
    J2EV = 6.24150934e18
    """
    Plot 2D wavefront (a slice).
    
    :param mwf: 2D wavefront structure , 
    :param title_fig: Figure title
    :param isHlog: if True, plot the horizontal cut in logarithmic scale
    :param isVlog: if True, plot the vertical cut in logarithmic scale
    :param i_x_min: Intensity threshold for horizontral cut, i.e. x-axis  limits are [min(where(i_x<i_x_min):max(where(i_x<i_x_min)]
    :param i_y_min: Intensity threshold for vertical cut,
    :param orient: 'x' for returning horizontal cut, 'y' for vertical cut
    :param onePlot: if True, put intensity map and plot of cuts on one plot, as  subplots
    :param bPlotPha: if True, plot the cuts of WF phase
    :return: 2-column array containing horizontal or vertical cut data in dependence of 'orient' parameter
    """
    
    [xc, yc] = calculate_peak_pos(mwf) # get centre

    if scale == "m":
        S = 1
        print('Coordinates of center, [m]:', xc * S, yc * S)
    elif scale == "mm":
        S = 1e3
        print('Coordinates of center, [mm]:', xc * S, yc * S)        
    elif scale == "um":
        S = 1e6
        print('Coordinates of center, [um]:', xc * S, yc * S)
    elif scale == "nm":
        S = 1e9
        print('Coordinates of center, [nm]:', xc * S, yc * S)
    else:
        assert("Scale should be m, mm, um or nm")


    

    ii = mwf.get_intensity(slice_number=0, polarization='horizontal')


    ii = ii*mwf.params.photonEnergy/J2EV#*1e3
    imax = np.max(ii)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    ph = mwf.get_phase(slice_number=0, polarization='horizontal')
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    
    
    
    xa = np.linspace(xmin, xmax, nx); 
    ya = np.linspace(ymin, ymax, ny); 
    
    if scale == "m":
        print('stepX, stepY [m]:', dx * S, dy * S, '\n')
    elif scale == "mm":
        print('stepX, stepY [mm]:', dx * S, dy * S, '\n')    
    elif scale == "um":\
        print('stepX, stepY [um]:', dx * S, dy * S, '\n')
    elif scale == "nm":
        print('stepX, stepY [nm]:', dx * S, dy * S, '\n')
        
    print('Total power (integrated over full range): %g [GW]' %(ii.sum(axis=0).sum(axis=0)*dx*dy*1e6*1e-9)) 
    print('Peak power calculated using FWHM:         %g [GW]' %(imax*1e-9*1e6*2*np.pi*(calculate_fwhm_x(mwf)/2.35)*(calculate_fwhm_y(mwf)/2.35)))
    print('Max irradiance: %g [GW/mm^2]'    %(imax*1e-9)) 
    label4irradiance = 'Irradiance (W/$mm^2$)'
    
    
    pylab.figure(figsize = (21,21))

# THIS PART NEEDS TO BE NEATENED....
    [x1, x2, y1, y2] = mwf.get_limits()
    
    
    
    pylab.imshow(ii, extent=[x1*S, x2*S, y1*S, y2*S])
    
    pylab.xlim(-width, width)
    pylab.ylim(-width, width)
############################################
    
    pylab.set_cmap('bone')
    #pylab.set_cmap('hot')
    pylab.axis('tight')
    #pylab.colorbar(orientation='horizontal')
    if scale == "m":
        pylab.xlabel('x (m)')
        pylab.ylabel('y (m)')
        
    if scale == "mm":
        pylab.xlabel('x (mm)')
        pylab.ylabel('y (mm)')
        
    if scale == "um":
        pylab.xlabel('x (um)')
        pylab.ylabel('y (um)')

    if scale == "nm":
        pylab.xlabel('x (nm)')
        pylab.ylabel('y (nm)')

    pylab.title(title_fig)
    
###########################################
    
    irr_y = ii[:, np.max(np.where(xa == xc))]
    irr_x = ii[np.max(np.where(ya == yc)), :]

    
    pylab.figure()
    
    pylab.plot(ya * S, irr_y)
    pylab.xlabel(('y ({})'.format(scale)))
    
    pylab.xlim(np.min(ya[np.where(irr_y >= imax * i_y_min)])
               * S, np.max(ya[np.where(irr_y >= imax * i_y_min)]) * S)
    pylab.ylim(0,np.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Vertical cut,  xc = ' + str(int(xc * 1e6)) + scale)
    pylab.grid(True)
    
    pylab.figure()
    
    pylab.plot(xa * S, irr_x)
    pylab.xlabel('x ({})'.format(scale))
    pylab.xlim(np.min(xa[np.where(irr_x >= imax * i_x_min)])
              * S, np.max(xa[np.where(irr_x >= imax * i_x_min)]) * S)

    pylab.ylim(0,np.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Horizontal cut, yc = ' + str(int(yc * S)) + scale)
    pylab.grid(True)
    
    
    
    dd = np.zeros(shape=(nx, 2), dtype=float)
    dd[:, 0] = xa
    #for idx in range(nx): dd[idx, 1] = sum(ii[:, idx])
    dd[:, 1] = irr_x
    plt.show()
    
    
    if Fineplot == True:
        pylab.figure()
        
        pylab.plot(xa * S, irr_x)
        pylab.xlabel('x ({})'.format(scale))
        pylab.xlim(X1*S, X2*S)
    
        pylab.ylim(0,np.max(ii)*1.1)
        pylab.ylabel(label4irradiance)
        pylab.title('Horizontal cut - Fineplot')
        pylab.grid(True)
        
        
        
        dd = np.zeros(shape=(nx, 2), dtype=float)
        dd[:, 0] = xa
        #for idx in range(nx): dd[idx, 1] = sum(ii[:, idx])
        dd[:, 1] = irr_x
        plt.show()
    return dd

#%%
def cuts_plot(mwf, title_fig, i_x_min, i_y_min, width, scale, X1, X2, Fineplot = None):
    from wpg.useful_code.wfrutils import calculate_peak_pos, get_mesh, calculate_fwhm_x, calculate_fwhm_y
    import pylab
    
    J2EV = 6.24150934e18
    """
    Plot 2D wavefront (a slice).
    
    :param mwf: 2D wavefront structure , 
    :param title_fig: Figure title
    :param isHlog: if True, plot the horizontal cut in logarithmic scale
    :param isVlog: if True, plot the vertical cut in logarithmic scale
    :param i_x_min: Intensity threshold for horizontral cut, i.e. x-axis  limits are [min(where(i_x<i_x_min):max(where(i_x<i_x_min)]
    :param i_y_min: Intensity threshold for vertical cut,
    :param orient: 'x' for returning horizontal cut, 'y' for vertical cut
    :param onePlot: if True, put intensity map and plot of cuts on one plot, as  subplots
    :param bPlotPha: if True, plot the cuts of WF phase
    :return: 2-column array containing horizontal or vertical cut data in dependence of 'orient' parameter
    """
    
    [xc, yc] = calculate_peak_pos(mwf) # get centre

    if scale == "m":
        S = 1
        print('Coordinates of center, [m]:', xc * S, yc * S)
    elif scale == "mm":
        S = 1e3
        print('Coordinates of center, [mm]:', xc * S, yc * S)        
    elif scale == "um":
        S = 1e6
        print('Coordinates of center, [um]:', xc * S, yc * S)
    elif scale == "nm":
        S = 1e9
        print('Coordinates of center, [nm]:', xc * S, yc * S)
    else:
        assert("Scale should be m, mm, um or nm")


    

    ii = mwf.get_intensity(slice_number=0, polarization='horizontal')


    ii = ii*mwf.params.photonEnergy/J2EV#*1e3
    imax = np.max(ii)
    [nx, ny, xmin, xmax, ymin, ymax] = get_mesh(mwf)
    ph = mwf.get_phase(slice_number=0, polarization='horizontal')
    dx = (xmax-xmin)/(nx-1); dy = (ymax-ymin)/(ny-1)
    
    
    
    xa = np.linspace(xmin, xmax, nx); 
    ya = np.linspace(ymin, ymax, ny); 
    
    if scale == "m":
        print('stepX, stepY [m]:', dx * S, dy * S, '\n')
    elif scale == "mm":
        print('stepX, stepY [mm]:', dx * S, dy * S, '\n')    
    elif scale == "um":\
        print('stepX, stepY [um]:', dx * S, dy * S, '\n')
    elif scale == "nm":
        print('stepX, stepY [nm]:', dx * S, dy * S, '\n')
        
    print('Total power (integrated over full range): %g [GW]' %(ii.sum(axis=0).sum(axis=0)*dx*dy*1e6*1e-9)) 
    print('Peak power calculated using FWHM:         %g [GW]' %(imax*1e-9*1e6*2*np.pi*(calculate_fwhm_x(mwf)/2.35)*(calculate_fwhm_y(mwf)/2.35)))
    print('Max irradiance: %g [GW/mm^2]'    %(imax*1e-9)) 
    label4irradiance = 'Irradiance (W/$mm^2$)'
     
    irr_y = ii[:, np.max(np.where(xa == xc))]
    irr_x = ii[np.max(np.where(ya == yc)), :]

    ##################################################################
    
    fig = pylab.subplots(1,3,figsize=(60,20))
    
    pylab.subplot(131)            
    pylab.plot(ya * S, irr_y)
    pylab.xlabel(('y ({})'.format(scale)))
    
    pylab.xlim(np.min(ya[np.where(irr_y >= imax * i_y_min)])
               * S, np.max(ya[np.where(irr_y >= imax * i_y_min)]) * S)
    pylab.ylim(0,np.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Vertical cut,  xc = ' + str(int(xc * 1e6)) + scale)
    pylab.grid(True)
    
    pylab.subplot(132)
    
    pylab.plot(xa * S, irr_x)
    pylab.xlabel('x ({})'.format(scale))
    pylab.xlim(np.min(xa[np.where(irr_x >= imax * i_x_min)])
              * S, np.max(xa[np.where(irr_x >= imax * i_x_min)]) * S)

    pylab.ylim(0,np.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Horizontal cut, yc = ' + str(int(yc * S)) + scale)
    pylab.grid(True)
    
    
    

    
    
    pylab.subplot(133)
    pylab.plot(xa * S, irr_x)
    pylab.xlabel('x ({})'.format(scale))
    pylab.xlim(X1*S, X2*S)

    pylab.ylim(0,np.max(ii)*1.1)
    pylab.ylabel(label4irradiance)
    pylab.title('Horizontal cut - Fineplot')
    pylab.grid(True)
    

    
    plt.savefig("out/" + title_fig + ".png")
    plt.show()



#%%

if __name__ == "__main__":
    plot_intensity()

    np_to_sample()
    opt_from_array()
    log()
    plot_treyfront()