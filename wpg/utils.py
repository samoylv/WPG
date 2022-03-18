# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import h5py
import sys
import numpy
from collections import defaultdict
import array
import collections

__author__ = 'A. Buzmakov'


def store_dict_hdf5(hdf5_file_name, input_dict):
    """ Store dictionary in hdf5 file.
    :param hdf5_file_name:
    :param input_dict:
    """

    def store_group(group, pearent_goup):
        """
        Store group (not a value)in HDF5 file.
        :param group:
        :param pearent_goup:
        """

        for (k, v) in list(group.items()):
            if isinstance(v, dict):
                if k not in pearent_goup:
                    tmp_group = pearent_goup.create_group(k)
                else:
                    tmp_group = pearent_goup[k]
                store_group(v, tmp_group)
            else:
                store_value(k, v, pearent_goup)

    def store_value(name, value, group):
        """
        Store value (scalar, string, array, etc.) in HDF5 file
        :param name:
        :param value:
        :param group:
        """
        if value is not None:
            if name in group:
                del group[name]
            try:
                if isinstance(name, bytes):
                    name = name.decode()
                if isinstance(value, str):
                    group[name] = value
                elif isinstance(value, numpy.ndarray):
                    group.create_dataset(name, data=value, chunks=True,
                            compression='gzip', compression_opts=1)    # compression='lzf'
                elif isinstance(value, (list, tuple)):
                    # Convert to numpy array.
                    value = numpy.array(value)

                    # Handle unicode as hdf does not like unicode arrays.
                    if value.dtype == numpy.dtype('<U31'):
                        value = [v.encode('utf8') for v in value]

                    # Create the dataset.
                    group.create_dataset(name,
                                         data=value,
                                         chunks=True,
                                         compression='gzip',
                                         compression_opts=1,
                                         )    # compression='lzf'
#                     if numpy.allclose(value, 0):
#                         group.create_dataset(name, data=value, chunks=True,
#                             compression='gzip', compression_opts=1)    # compression='lzf'
#                     else:
#                         if all([isinstance(a, numpy.bytes_) for a in value]):
#                             value = numpy.array([a.decode('utf-8') for a in value])
#                         group.create_dataset(name, data=value)
#                 elif isinstance(value, str):
#                         group.create_dataset(name, data=numpy.array(value.encode('utf-8')))

#                 elif isinstance(value, list):
#                     if all([isinstance(a, str) for a in value]):
#                         value = numpy.array([a.encode('utf-8') for a in value])
#                     else:
#                         value = numpy.array(value)
#                     group.create_dataset(name, data=numpy.array(value))
                else:
                    group.create_dataset(name, data=value)
            except ValueError:  # if h5py not support compression
                group.create_dataset(name, data=value)
            except TypeError:
                group.create_dataset(name, data=value)
            except Exception:
                print("Error at name='{}' value='{}' group='{}'".format(name, value, group))
                raise

    with h5py.File(hdf5_file_name, 'w') as res_file:
        store_group(input_dict, res_file)


def load_dict_slash_hdf5(hdf5_file_name):
    """Load dictionary with fields of wavefront
    :param hdf5_file_name: Input file name
    """
    out_dict = {}

    def add_item(name, obj):
        if not isinstance(obj, h5py.Group):
            out_dict.update({name.encode(): obj[()]})
    with h5py.File(hdf5_file_name, 'r') as h5_file:
        h5_file.visititems(add_item)

    return out_dict


def update_dict_slash_string(input_dict, keys_string, value):
    """
    Update dictionary from slash separated keys_string by value
    :param input_dict: dictionary to be updated
    :param keys_string: slash separated keys_string
    :param value: value
    """
    try:
        keys_string = keys_string.decode('utf-8')
    except AttributeError:
        pass

    keys = keys_string.split('/')
    tdict = input_dict
    for k in keys[:-1]:
        if k not in tdict:
            tdict[k] = {}
        tdict = tdict[k]
    tdict[keys[-1]] = value


def tree():
    """
    Simple tree
    :return:
    """
    return defaultdict(tree)


def dicts(t):
    """
    Convert tree to dict
    :param t:
    :return:
    """
    try:
        return dict((k, dicts(t[k])) for k in t)
    except TypeError:
        return t


def set_value(dic, keys_chain, value):
    """
    Set value in dictionary by  chain of keys

    :param value: value
    :param dic: dic
    :param keys_chain: list of keys
    """
    node = dic
    for key in keys_chain[:-1]:
        if key in node:
            node = node[key]
        else:
            node[key] = {}
            node = node[key]
    node[keys_chain[-1]] = value


def get_value(dic, keys_chain):
    """
    Get node from dictionary by  chain of keys

    :param dic: dict from which value will taken
    :param keys_chain: list of keys
    :return res: return value
    """
    res = dic
    for key in keys_chain:
        res = res[key]
    return res


def set_value_attr(obj, keys_chain, value):
    """
    Return node from dictionary by  chain of keys

    :param obj:
    :param keys_chain: list of keys
    :param value:
    """

    class glossary_folder(object):

        """Glossary folder"""
        pass

    node = obj
    for key in keys_chain[:-1]:
        if key not in node.__dict__:
            node.__dict__[key] = glossary_folder()
        node = node.__dict__[key]
    node.__dict__[keys_chain[-1]] = value


def load_hdf5_tree(hdf5_file_name):
    """convert h5file to the tree"""
    out_dict = {}

    def add_item(name, obj):
        if not isinstance(obj, h5py.Group):
            tmp = {}
            if type(obj.value) == numpy.ndarray:
                tmp['value'] = 'array shape = {}'.format(obj.shape)
            else:
                tmp['value'] = obj.value
            tmp['attrs'] = dict((k.encode(), v) for k, v in list(obj.attrs.items()))
            out_dict.update({name.encode(): tmp})

    with h5py.File(hdf5_file_name, 'r') as h5_file:
        h5_file.visititems(add_item)

    return out_dict


def print_hdf5(hdf5_file_name):
    """print h5file to the tree"""
    from pprint import PrettyPrinter

    class MyPrettyPrinter(PrettyPrinter):

        def format(self, *args, **kwargs):
            repr, readable, recursive = PrettyPrinter.format(
                self, *args, **kwargs)
#             if repr:
#                 if repr[0] in ('"', "'"):
#                     repr = repr.decode('string_escape')
#                 elif repr[0:2] in ("u'", 'u"'):
#                     repr = repr.decode('unicode_escape').encode(
#                         sys.stdout.encoding)
            return repr, readable, recursive

    def pprint(obj, stream=None, indent=1, width=80, depth=None):
        printer = MyPrettyPrinter(
            stream=stream, indent=indent, width=width, depth=depth)
        printer.pprint(obj)

    def pformat(obj, stream=None, indent=1, width=80, depth=None):
        printer = MyPrettyPrinter(
            stream=stream, indent=indent, width=width, depth=depth)
        return printer.pformat(obj)

    pprint(load_hdf5_tree(hdf5_file_name))


def srw_obj2str(obj, start_str=''):
    fields = [field for field in dir(obj) if not field.startswith(
        '__') if not isinstance(getattr(obj, field), collections.Callable)]
    res = ''
    for f in fields:
        val = getattr(obj, f)
        if 'glossary_folder' in str(val):
            continue
        s = ''
        if not isinstance(val, array.array):
            if 'srwlib.' in str(val):
                s = val.__doc__ + '\n' + srw_obj2str(val, start_str + '\t')
            else:
                s = str(val)
        else:
            s = 'array of size ' + str(len(val))
        res += '{0}{1} = {2}\n'.format(start_str, f, s)
    return res
