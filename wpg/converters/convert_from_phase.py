__author__ = 'makov'

import h5py
import numpy

def convert_from_phase(in_file_name,out_file_name):

    def convert_item(name,obj):
        node_name=name.encode()
        if not isinstance(obj, h5py.Group):
            if node_name.startswith('data/'):
                if node_name=='data/arrEhor':
                    node_name='data/arrEx'
                elif node_name=='data/arrEver':
                    node_name='data/arrEy'
                else:
                    raise ValueError('wrong node name: '+ node_name)

                s=obj.shape
                h5out.create_dataset(node_name, shape=(s[1],2*s[2],s[3]), chunks=True,
                    compression='gzip', compression_opts=9)
                h5out[node_name][:,0::2,:]=obj[0]
                h5out[node_name][:,1::2,:]=obj[1]
            elif node_name=='VERSION':
                node_name='version'
                h5out[node_name]=obj.value
            else:
                node_name='header/'+node_name
                if obj.dtype==numpy.dtype('|S1'):
                    #if coloumn string
                    val=''.join(obj.value)
                    if node_name=='header/wEFieldUnit' and val=='GW/cm^2':
                        val='arbitrary'
                    if node_name== 'header/wFloatType':
                        val='float'
                    h5out[node_name]=val
                else:
                    h5out[node_name]=obj.value

    with h5py.File(in_file_name,'r') as h5in:
        with h5py.File(out_file_name,'w') as h5out:
            h5in.visititems(convert_item)

if __name__=="__main__":
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-i", "--input-file", dest="in_filename",
        help="Input file name (phase)", metavar="FILE")
    parser.add_option("-o", "--output-file", dest="out_filename",
        help="Outpu file name (srwpy)", metavar="FILE")
    (options, args) = parser.parse_args()

    convert_from_phase(options.in_filename,options.out_filename)