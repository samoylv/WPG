# -*- coding: utf-8 -*-

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

__author__ = 'A. Buzmakov'


import os
import sys
sys.path.insert(0, os.path.join('..', '..'))

import numpy as np

import wpg
from wpg.generators import build_gauss_wavefront
from wpg.beamline import Beamline
from wpg.optical_elements import Drift, Use_PP
from wpg.srwlib import srwl


def main():
    d2waist = 270.
    # beam parameters:
    qnC = 0.1  # [nC] e-bunch charge
    thetaOM = 3.6e-3
    ekev = 5.0

    # calculate angular divergence:
    theta_fwhm = (17.2 - 6.4 * np.sqrt(qnC)) * 1e-6 / ekev ** 0.85
    theta_rms = theta_fwhm / 2.35
    sigX = 12.4e-10 / (ekev * 4 * np.pi * theta_rms)

    # define limits
    xmax = theta_rms * d2waist * 3.5
    xmin = - xmax
    ymin = xmin
    ymax = xmax
    nx = 600
    ny = nx
    nz = 5
    tau = 0.12e-15

    srw_wf = build_gauss_wavefront(
        nx, ny, nz, ekev, xmin, xmax, ymin, ymax, tau, sigX, sigX, d2waist)
    wf = wpg.Wavefront(srw_wf)  
    b = Beamline()
    b.append(Drift(5), Use_PP())
    b.propagate(wf)
    # srwl.ResizeElecField(srw_wf, 'c', [0, 0.25, 1, 0.25, 1])

    if not os.path.exists('tests_data'):
        os.mkdir('tests_data')

    wf_hdf5_out_file_path = os.path.join('tests_data', 'my_gauss.h5')
    wf.store_hdf5(wf_hdf5_out_file_path)

    wf_out = wpg.Wavefront()
    wf_out.load_hdf5(wf_hdf5_out_file_path)
    return wf

if __name__ == "__main__":
    main()
