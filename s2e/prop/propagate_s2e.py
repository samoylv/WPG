# -*- coding: utf-8 -*-
""" Module for interface functions to wpg library. """

from glob import glob
from wpg import Wavefront, Beamline
import errno
import h5py
import multiprocessing
import numpy
import os
import wpg

MIRROR_DATA_DIR = os.path.join(os.path.abspath(os.path.dirname(__file__)), 'data_common/')

def mkdir_p(path):
    """
    Create directory tree, if not exists (mkdir -p)

    :param path: Path to be created
    """
    if path == '':
        return
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


def add_history(wf_file_name, history_file_name):
    """
    Add history from pearent file to propagated file

    :param wf_file_name: output file
    :param history_file_name: peraent file
    """
    with h5py.File(wf_file_name) as wf_h5:
        with h5py.File(history_file_name) as history_h5:
            if 'history' in wf_h5:
                del wf_h5['history']

            wf_h5.create_group('/history/parent/')
            wf_h5.create_group('/history/parent/detail')

            for k in history_h5:
                if k=='history':
                    try:
                        history_h5.copy('/history/parent', wf_h5['history']['parent'])
                    except KeyError:
                        pass

                elif not k == 'data':
                    history_h5.copy(k,wf_h5['history']['parent']['detail'])
                else:
                    wf_h5['history']['parent']['detail']['data'] = h5py.ExternalLink(history_file_name,'/data')


def calculate_fwhm(wfr):
    """
    Calculate FWHM of the beam calculating number of point bigger then max/2 throuhgt center of the image

    :param wfr:  wavefront
    :return: {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y} in [m]
    """
    intens = wfr.get_intensity(polarization='total').sum(axis=-1);


    mesh = wfr.params.Mesh
    dx = (mesh.xMax-mesh.xMin)/mesh.nx
    dy = (mesh.yMax-mesh.yMin)/mesh.ny

    x_center = intens[intens.shape[0]//2,:]
    fwhm_x = len(x_center[x_center>x_center.max()/2])*dx

    y_center = intens[:,intens.shape[1]//2]
    fwhm_y = len(y_center[y_center>y_center.max()/2])*dy
    return {'fwhm_x':fwhm_x, 'fwhm_y': fwhm_y}

def get_intensity_on_axis(wfr):
    """
    Calculate intensity (spectrum in frequency domain) along z-axis (x=y=0)

    :param wfr:  wavefront
    :return: [z,s0] in [a.u.] if frequency domain
    """

    wf_intensity = wfr.get_intensity(polarization='horizontal')
    mesh = wfr.params.Mesh;
    zmin = mesh.sliceMin;
    zmax = mesh.sliceMax;
    sz = numpy.zeros((mesh.nSlices, 2), dtype='float64')
    sz[:,0] = numpy.linspace(zmin, zmax, mesh.nSlices);
    sz[:,1] = wf_intensity[mesh.nx/2, mesh.ny/2, :] / wf_intensity.max()

    return sz

def propagate(in_fname, out_fname, get_beamline):
    """
    Propagate wavefront

    :param in_file: input wavefront file
    :param out_file: output file
    :param get_beamline: function to build beamline
    """

    print('#'*80)
    print("Setup initial wavefront.")
    wf=Wavefront()

    # Load wavefront data.
    print('Load ' + in_fname)
    wf.load_hdf5(in_fname)

    # Get beamline.
    bl0 = get_beamline()

    # Switch to frequency domain.
    wpg.srwlib.srwl.SetRepresElecField(wf._srwl_wf, 'f')

    # Save spectrum for later reference.
    sz0 = get_intensity_on_axis(wf);
    wf.custom_fields['/misc/spectrum0'] = sz0

    # Propagate.
    bl0.propagate(wf)

    # Save spectrum after propagation for later reference.
    sz1 = get_intensity_on_axis(wf);
    wf.custom_fields['/misc/spectrum1'] = sz1

    # Switch back to time domain.
    wpg.srwlib.srwl.SetRepresElecField(wf._srwl_wf, 't')

    #Resizing: decreasing Range of Horizontal and Vertical Position:
    wpg.srwlib.srwl.ResizeElecField(wf._srwl_wf, 'c', [0, 0.5, 1, 0.5,  1]);

    add_custom_data(wf, bl0)

    print('Saving propagated wavefront to ' + out_fname)
    mkdir_p(os.path.dirname(out_fname))
    wf.store_hdf5(out_fname)

    print('Saving history.')
    add_history(out_fname, in_fname)

    print('ALL DONE.')
    print('#'*80)



def stepwise(in_fname, get_beamline):
    """
    Propagate wavefront stepwise, dumping the wavefront at every step.

    :param in_file: input wavefront file
    :param get_beamline: function to build beamline
    """

    print('#'*80)
    print("Setup initial wavefront.")
    wf=Wavefront()

    # Load wavefront data.
    print('Load ' + in_fname)
    wf.load_hdf5(in_fname)

    # Get beamline.
    bl0 = get_beamline()

    beamline = bl0.propagation_options
    if len(beamline) > 1:
        raise RuntimeError("Beamline configuration not supported.")
    beamline = beamline[0]
    elements = beamline['optical_elements']
    options = beamline['propagation_parameters']
    if len(elements) != len(options):
        raise RuntimeError("Beamline configuration not supported.")

    i = 0
    for element, option in zip(elements, options):

        print('\n')
        print('#'*80)
        print('Propagation step %d.' %(i))
        print('Setting up incremental beamline.')
        beamline_step = Beamline()
        beamline_step.append(element, option) ### <== CHECKME

        # Switch to frequency domain.
        wpg.srwlib.srwl.SetRepresElecField(wf._srwl_wf, 'f')

        # Save spectrum for later reference.
        sz0 = get_intensity_on_axis(wf);
        wf.custom_fields['/misc/spectrum0'] = sz0

        # Propagate.
        beamline_step.propagate(wf)

        # Save spectrum after propagation for later reference.
        sz1 = get_intensity_on_axis(wf);
        wf.custom_fields['/misc/spectrum1'] = sz1

        # Switch back to time domain.
        wpg.srwlib.srwl.SetRepresElecField(wf._srwl_wf, 't')

        incremental_filename = "%04d.h5" % (i)
        print('Saving propagated wavefront to ' + incremental_filename)
        mkdir_p(os.path.dirname(incremental_filename))
        wf.store_hdf5(incremental_filename)

        print('Done with propagation step %d.' %(i))
        print('#'*80)

        # Increment running index.
        i += 1



def add_custom_data(wf, bl0):
    """ Utility to store additional info on the wavefront . """

    # Add more metadata.
    fwhm = calculate_fwhm(wf)


    wf.custom_fields['/misc/xFWHM'] = fwhm['fwhm_x']
    wf.custom_fields['/misc/yFWHM'] = fwhm['fwhm_y']
    wf.custom_fields['/params/beamline/printout'] = str(bl0)

    wf.custom_fields['/info/contact'] = [
        'Name: Liubov Samoylova', 'Email: liubov.samoylova@xfel.eu',
        'Name: Alexey Buzmakov', 'Email: buzmakov@gmail.com']
    wf.custom_fields['/info/data_description'] = 'This dataset contains infromation about wavefront propagated through beamline (WPG and SRW frameworks).'
    wf.custom_fields['/info/method_description'] = """WPG, WaveProperGator (http://github.com/samoylv/WPG)is an interactive simulation framework for coherent X-ray wavefront propagation.\nSRW, Synchrotron Radiation Workshop (http://github.com/ochubar/SRW),  is a physical optics computer code  for simulation of the radiation wavefront propagation through optical systems of beamlines as well as  detailed characteristics of Synchrotron Radiation (SR) generated by relativistic electrons in magnetic fields of arbitrary configuration."""
    wf.custom_fields['/info/package_version'] = '2014.1'


def propagate_wrapper(params):
    """
    Wrapper for passing parameters as a tupple from multiprocessing module
    """
    (in_fname, out_fname, get_beamline) = params
    return propagate(in_fname, out_fname, get_beamline)

def directory_process(in_dname, out_dname, get_beamline, cpu_number):
    """
    Process directory with in_dname\FELsource_out*.h5 files and store it after propagation in out_dname\prop_out*.h5 files

    :param in_dname: input directory name
    :param out_dname: ouput directory name
    :param get_beamline: function to build beamline
    :param cpu_number: NUmber of CPUs for parallel computing

    """
    input_dir = in_dname
    input_files = glob(os.path.join(input_dir, 'FELsource_out*.h5'))
    out_files = []
    for name in input_files:
        in_file_name = os.path.split(name)[-1]
        out_file_name = in_file_name.replace('FELsource_out','prop_out')
        print ('out_file_name: %s' % out_file_name)
        out_files.append(os.path.join(out_dname, out_file_name))

    print('Found {} HDF5 files in {}'.format(len(input_files), in_dname))

    batch_params = zip(input_files, out_files, [get_beamline]*len(input_files))

    p=multiprocessing.Pool(processes=cpu_number)
    p.map(propagate_wrapper, batch_params, chunksize=1)
    p.close()
    p.join()

def main():
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--input-file", dest="in_fname", help="Input wavefront file")
    parser.add_argument("--output-file", dest="out_fname", help="Output wavefront file")

    parser.add_argument("--input-directory", dest="in_dname", help="Input directory with wavefront files")
    parser.add_argument("--output-directory", dest="out_dname", help="Output directory with wavefront files")
    parser.add_argument("-n", "--cpu-number", dest="cpu_number", default=int((multiprocessing.cpu_count()+1)/2),
               help="Number of cores for batch wavefronts propagation, default value NUMBER_OF_CPU/2")
    parser.add_argument("--beamline-file", dest="beamline_file", help="Python file with beamline description")

    options = parser.parse_args()


    if not (options.in_fname or options.in_dname):   # if filename is not given
        parser.error('Input filename or directiry not specified, use --input-file or --input-directory options')
        return

    if not (options.out_fname or options.out_dname):   # if filename is not given
        parser.error('Output filename or directiry not specified, use --output-file or --output-directory options')
        return

    if not(options.beamline_file):
        parser.error('Beamline description file specified, use --beamline-file option')
        return
    else:
        import imp
        custom_beamline = imp.load_source('custom_beamline', options.beamline_file)
        get_beamline = custom_beamline.get_beamline

    if options.in_dname and options.out_dname:
        print('Input directory {}, output directory {}, number of cores {}'.format(
            options.in_dname, options.out_dname, options.cpu_number))
        print('Batch propagation started')
        directory_process(options.in_dname, options.out_dname, get_beamline, int(options.cpu_number))
        print('Batch propagation finished')

    elif options.in_fname and options.out_fname:
        print('Input file {}, output file {}'.format(options.in_fname, options.out_fname, get_beamline))
        propagate(options.in_fname, options.out_fname, get_beamline)

if __name__ == '__main__':
    main()
