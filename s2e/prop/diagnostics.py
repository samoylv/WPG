import sys
#sys.path.insert(0, '../WPG') # AB
#sys.path.insert(0,'/diskmnt/a/lsamoylv/WPG')
#sys.path.insert(0,'/Users/lsamoylv/code/WPG')
sys.path.insert(0,'/data/S2E/packages/WPG/')

import os
import pylab as plt

import numpy as np
from wpg import Wavefront


def show_diagnostics(prop_out_number):
    # read prop_out_.h5
    if not prop_out_number == 'prop_out_1.h5':
        #prop_out_file = "prop_out_{}.h5".format(prop_out_number.zfill(7))
        prop_out_file = "{}.h5".format(prop_out_number.zfill(7))
    else:
        prop_out_file = prop_out_number

    if not os.path.exists(prop_out_file):
        print 'Input file {} not found.'.format(prop_out_file)
        return

    wf = Wavefront()
    wf.load_hdf5(prop_out_file)
    # show two figures window 1: image of I(x,y) integral intensity, with real
    # x and y axis and title with file name
    J2eV = 6.24150934e18
    mesh = wf.params.Mesh
    tmin = mesh.sliceMin
    tmax = mesh.sliceMax
    dt = (tmax - tmin) / (mesh.nSlices - 1)
    dx = (mesh.xMax - mesh.xMin) / (mesh.nx - 1)
    dy = (mesh.yMax - mesh.yMin) / (mesh.ny - 1)

    wf_intensity = wf.get_intensity(polarization='horizontal')
    total_intensity = wf_intensity.sum(axis=-1)
    data = total_intensity * dt
    plt.figure()
    plt.imshow(data*dx*dy*1e6*J2eV/wf.params.photonEnergy,extent=[mesh.xMin*1e6,mesh.xMax*1e6,mesh.yMin*1e6,mesh.yMax * 1e6])
    title = 'Number of photons per %.2f x %.2f $nm ^2$ pixel'  %  (dx*1e9, dx*1e9)
    plt.title(title)
    plt.colorbar()
    plt.xlabel(r'[$\mu$m]')

    # window 2: plot of 2 curves:
    #(1) history/parent/parent/temporal_struct - before propagating
    temporal_struct = wf.custom_fields['history']['parent']['parent']['misc']['temporal_struct']
    t0 = (temporal_struct[:, 0].max() + temporal_struct[:, 0].min()) / 2

    plt.figure()
    plt.plot(temporal_struct[:, 0] - t0, temporal_struct[:, 1] * 1e-9, 'b',label = 'original')
    plt.hold(True)
    #(2) integral intensity I(t) after propagating

    t = np.linspace(tmin, tmax, wf.params.Mesh.nSlices)
    pulse_energy = wf.get_intensity().sum(axis=0).sum(axis=0) #check it
    plt.plot(t * 1e15, pulse_energy*dx*dy*1e6*1e-9,'r', label = 'propag')

    title = 'The propagated pulse energy %.2f %s ' % (pulse_energy.sum(axis=0) * dx * dy * 1e6 * dt * 1e3, 'mJ')
    plt.title(title)
    plt.xlabel('time [fs]')
    plt.ylabel('Instantaneous power [GW]')
    plt.legend()
    plt.grid(True)

    sz0 = wf.custom_fields['misc']['spectrum0']
    sz1 = wf.custom_fields['misc']['spectrum1']
    plt.figure()
    plt.plot(sz0[:,0],sz0[:,1], label='before propagating')
    plt.hold(True)
    plt.plot(sz1[:,0],sz1[:,1],label='after propagating')
    plt.grid(True)
    plt.title('Spectrum (x=y=0)')
    plt.xlabel('[eV]')
    plt.ylabel('[arb. unit]')
    plt.legend()

    plt.show()

def main():
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("--input-file", dest="in_fname", default="prop_out_1.h5",help="Input wavefront file: prop_out_***.h5")
    (options, args) = parser.parse_args()

    if not options.in_fname :   # if filename is not given
        parser.error('Input filename not specified, use --input-file options')
        return

    show_diagnostics(options.in_fname)

if __name__ == "__main__":
    main()
