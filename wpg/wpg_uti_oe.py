# -*- coding: utf-8 -*-
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import numpy as np
import pylab as plt

__author__ = 'A. Buzmakov, L. Samoylova'




def get_opd(transmission):
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    phase = np.array(transmission.arTr[1::2]).reshape(
        (ny, nx))  # check, ny,nx order
    return phase


def get_absorption(transmission):
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    amplitud = np.array(transmission.arTr[::2]).reshape(
        (ny, nx))  # check, ny,nx order
    return amplitud


def show_transmission(transmission):
    mesh = transmission.mesh

    plt.figure(figsize=(10, 7))
    plt.subplot(121)

    plt.imshow(get_absorption(transmission),
               extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart),
               cmap=plt.cm.bone)
    plt.colorbar(orientation='horizontal')
    plt.title('Absorption')

    plt.subplot(122)
    plt.imshow(get_opd(transmission),
               extent=(mesh.xStart, mesh.xFin, mesh.yFin, mesh.yStart),
               cmap=plt.cm.bone)
    plt.colorbar(orientation='horizontal')
    plt.title('OPD [m]')
