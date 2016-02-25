# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__author__ = 'A. Buzmakov, L. Samoylova'

import numpy as np
import pylab as plt

def get_opd(transmission):
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    phase = np.array(transmission.arTr[1::2]).reshape((ny,nx)) # check, ny,nx order
    return phase

def get_absorption(transmission):
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    amplitud = np.array(transmission.arTr[::2]).reshape((ny,nx)) # check, ny,nx order
    return amplitud

def show_transmission(transmission):
    mesh = transmission.mesh
    #nx = mesh.nx
    #ny = mesh.ny
    #print nx*ny*2
    #numpy_data = np.array(transmission.arTr).reshape((ny,nx,2)) # check, ny,nx order
    plt.figure(figsize=(10,7))
    plt.subplot(121)
    #plt.imshow(numpy_data[...,0], extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart))
    plt.imshow(get_absorption(transmission), extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart))
    plt.colorbar(orientation='horizontal');
    plt.set_cmap('bone')
    plt.title('Absorption')
    plt.subplot(122)
    #plt.imshow(numpy_data[...,1], extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart))
    plt.imshow(get_opd(transmission), extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart))
    plt.colorbar(orientation='horizontal');
    plt.set_cmap('bone')
    plt.title('OPD [m]')
    
#show_transmission(opTrErM1)