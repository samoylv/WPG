# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

__author__ = 'A. Buzmakov, L. Samoylova'

def get_opd(transmission):
    """
    Extract optical path difference of  transmission object

    :param transmission: SRWLOptT struct, see srwlib.h
    :return: optical path difference map of transmission object
    """
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    phase = np.array(transmission.arTr[1::2]).reshape(
        (ny, nx))  # check, ny,nx order
    return phase


def get_absorption(transmission):
    """
    Extract absorption coefficient of  transmission object

    :param transmission: SRWLOptT struct, see srwlib.h
    :return: absorption(x,y) map of transmission object
    """
    mesh = transmission.mesh
    nx = mesh.nx
    ny = mesh.ny
    amplitud = np.array(transmission.arTr[::2]).reshape(
        (ny, nx))  # check, ny,nx order
    return amplitud


def show_transmission(transmission):
    """
    Plot absorption and OPD maps of transmission object

    :param transmission: SRWLOptT struct, see srwlib.h
    """
    mesh = transmission.mesh

    plt.figure(figsize=(10, 7))
    plt.subplot(121)

    plt.imshow(get_absorption(transmission),
               extent=(mesh.xStart*1e3, mesh.xFin*1e3, mesh.yFin*1e3, mesh.yStart*1e3),
               cmap=plt.cm.bone, aspect='auto')
    plt.colorbar(orientation='horizontal')
    plt.title('Transmission');plt.xlabel('mm')

    plt.subplot(122)
    plt.imshow(get_opd(transmission),
               extent=(mesh.xStart*1e3, mesh.xFin*1e3, mesh.yFin*1e3, mesh.yStart*1e3),
               cmap=plt.cm.bone, aspect='auto')
    plt.colorbar(orientation='horizontal')
    plt.title('OPD [m]');plt.xlabel('mm')
    plt.show()


def propagation_parameters(sx, zx, sy, zy, mode = "fresnel"):
    """
    wrapper for propagation parameters
    
    :param sx: horizontal scaling factor 
    :param zx: horizontal zoom factor
    :param sy: vertical scaling factor
    :param zy: vertical zoom factor
    :param mode: treatment of the quadratic phase terms, SRW propagation parameter [3]
        (cpp/src/core/sroptdrf.h, ChooseLocalPropMode), case-insensitive:
        'fresnel' - 0, standard propagation through the angular representation;
        'quadratic' - 1, semi-analytical treatment of the quadratic phase;
        'fraunhofer' - 2, the second variant of the semi-analytical treatment of the quadratic
        phase in SRW (not a far-field propagator, the name is kept for compatibility);
        'diverge' - 3, propagation from waist; 'converge' - 4, propagation to waist

    :return propagation parameters:
    """
    modes = {"fresnel": 0, "quadratic": 1, "fraunhofer": 2, "diverge": 3, "converge": 4}
    if str(mode).lower() not in modes:
        raise ValueError('mode must be one of {}, given: {!r}'.format(list(modes), mode))
    m = modes[str(mode).lower()]

    return [0,0,1,m,0,sx,zx/sx,sy,zy/sy,0,0,0]

