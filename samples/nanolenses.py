
# coding: utf-8

# In[ ]:

get_ipython().magic(u'pylab inline')


# In[ ]:

import os
import sys
import numpy as np
import pylab as plt
import h5py
import datetime

sys.path.insert(0,'..')


# In[ ]:

import wpg
from wpg.converters.nanolenses import EdfFile


## Reading from EDF

# In[ ]:

def get_array_from_edf(edf_file_name):
    print 'Processing EDF file {}'.format(edf_file_name)
    edf_structure = EdfFile.EdfFile(edf_file_name,'r')
    edf_header = edf_structure.GetHeader(0)
    print 'Number of images in EDF file {}'.format(edf_structure.GetNumImages())
    print 'EDF header:'
    for k,v in edf_header.items():
        print '\t{:20}:\t{}'.format(k,v)
    print edf_structure.GetData(0).dtype
    image_data = edf_structure.GetData(0).astype('float32')
    image_min = np.float32(edf_header['TUD_ScalingMin'])
    image_max = np.float32(edf_header['TUD_ScalingMax'])
    
    image_data_scaled = (image_data-image_data.min())/(image_data.max()-image_data.min())
    image_data_scaled = image_data_scaled*(image_max-image_min)+image_min
    return image_data_scaled

def get_wavefront_from_edf(edf_real_file_name, edf_imag_file_name):
    wavefront_real = get_array_from_edf(edf_real_file_name)
    wavefront_imag = get_array_from_edf(edf_imag_file_name)
    res={'real': wavefront_real, 'imag': wavefront_imag}
    return res
    


# In[ ]:

nano_re_file = '../wpg/converters/nanolenses/sources/testillu_re.edf'
nano_im_file = '../wpg/converters/nanolenses/sources/testillu_im.edf'


# In[ ]:

wf = get_wavefront_from_edf(nano_re_file, nano_im_file)


### Where is th mesh and polarization?

# In[ ]:

plt.figure(figsize=(15,5))
plt.subplot(131)
plt.imshow(wf['real'])
plt.colorbar()
plt.title('Real part')
plt.subplot(132)
plt.imshow(wf['imag'])
plt.colorbar()
plt.title('Imaginary part')
plt.subplot(133)
plt.imshow(np.sqrt(wf['real']**2+wf['imag']**2))
plt.colorbar()
plt.title('Intensity')


## Writing EDF

# In[ ]:

wf = wpg.Wavefront()
wf.load_hdf5('../wpg/tests/tests_data/my_gauss.h5')


# In[ ]:

wf.params.Mesh.nSlices


# In[ ]:

real_part = wf.get_real_part(polarization='vertical').sum(axis=-1)
plt.imshow(real_part)


### What do with polarisation?

# In[ ]:

imag_part = wf.get_imag_part(polarization='vertical').sum(axis=-1)
plt.imshow(imag_part)
plt.colorbar()


# In[ ]:

def store_array_to_edf(in_array, file_name):
    array_min = in_array.min()
    array_max = in_array.max()
    out_array = (in_array - array_min)/(array_max - array_min)*np.iinfo(np.uint16).max
    out_array = out_array.astype('uint16')
    edf_headers = {
        'TUD_ScalingMin':array_min,
        'TUD_ScalingMax':array_max,
        'TUD_MatlibVersion':'tomo_version_9.0',
        'TUD_FilenameCurr':file_name,
        'TUD_FilenamePrev':file_name,
        'VersionNumber':'edf_version_10.5',
        'TUD_DateCurr':str(datetime.datetime.now()),
        'TUD_DatePrev':str(datetime.datetime.now()),
        'TUD_Step_y':'1',
        'TUD_Step_x':'1',
        'TUD_FileType':'imagedata',
        'TUD_FilenameOrig':file_name,
        'TUD_HeaderLines':'26',
        'TUD_DateOrig':str(datetime.datetime.now()),
        'TUD_HeaderSize':'2048',
        'TUD_RotAxisPosY':'-10000',
        'TUD_RotAxisPosX':'-1000'
    }
    #todo add custom headers
    edf_file = EdfFile.EdfFile(file_name)
    edf_file.WriteImage(edf_headers, out_array,0)


# In[ ]:

edf_file = store_array_to_edf(imag_part,'test.edf')


# In[ ]:

res = get_array_from_edf('test.edf')


# In[ ]:

plt.imshow(res)
plt.colorbar()

