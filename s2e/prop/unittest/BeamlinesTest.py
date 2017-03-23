##########################################################################
#                                                                        #
# Copyright (C) 2017 Carsten Fortmann-Grote                              #
# Contact: Carsten Fortmann-Grote <carsten.grote@xfel.eu>                #
#                                                                        #
# This file is part of WPG/s2e/prop                                      #
# and is free software: you can redistribute it and/or modify            #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# It is distributed in the hope that it will be useful,                  #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.  #
#                                                                        #
##########################################################################

""" Test module for the XFELPhotonPropagator.

    @author : CFG
    @institution : XFEL
    @creation 20170317

"""
import os, sys, shutil
import numpy
import unittest


# Import the class to test.
sys.path.insert(0,'../../..')


from s2e.prop import simple_beamline, exfel_spb_day1_beamline, exfel_spb_kb_beamline
from s2e.prop import propagate_s2e

from wpg import Beamline, Wavefront
from wpg.generators import build_gauss_wavefront


class BeamlinesTest(unittest.TestCase):
    """
    Test class for the beamline modules in prop.
    """

    @classmethod
    def setUpClass(cls):
        """ Setting up the test class. """

        setupTestWavefront()
        cls.__fel_source = "source.h5"

    @classmethod
    def tearDownClass(cls):
        """ Tearing down the test class. """
        os.remove("source.h5")
        del cls.__fel_source

    def setUp(self):
        """ Setting up a test. """
        self.__files_to_remove = []
        self.__dirs_to_remove = []

    def tearDown(self):
        """ Tearing down a test. """

        for f in self.__files_to_remove:
            if os.path.isfile(f): os.remove(f)
        for d in self.__dirs_to_remove:
            if os.path.isdir(d): shutil.rmtree(d)

    def testSimpleBeamline(self):
        """ Testing the construction of a simple beamline. """
        bl = simple_beamline.get_beamline()

        self.assertIsInstance(bl, Beamline)

    def testSimpleBeamlinePropagation(self):
        """ Testing propagation through the simple beamline. """

        # Setup output and cleanup.
        output_file = "prop_out.h5"
        self.__files_to_remove.append(output_file)

        # Propagate.
        propagate_s2e.propagate(self.__fel_source, output_file, simple_beamline.get_beamline )

        # Check that output was generated.
        self.assertIn( output_file, os.listdir(".") )

    def testEXFELSPBDay1Beamline(self):
        """ Testing the construction of the Day1 EUXFEL SPB-SFX beamline. """
        bl = exfel_spb_day1_beamline.get_beamline()

        self.assertIsInstance(bl, Beamline)

    def testEXFELSPBDay1BeamlinePropagation(self):
        """ Testing propagation through the Day1 EUXFEL SPB-SFX beamline. """

        # Setup output and cleanup.
        output_file = "prop_out.h5"
        self.__files_to_remove.append(output_file)

        # Propagate.
        propagate_s2e.propagate(self.__fel_source, output_file, exfel_spb_day1_beamline.get_beamline )

        # Check that output was generated.
        self.assertIn( output_file, os.listdir(".") )

    def testEXFELSPBKBBeamline(self):
        """ Testing the construction of the EUXFEL SPB-SFX KB beamline. """
        bl = exfel_spb_kb_beamline.get_beamline()

        self.assertIsInstance(bl, Beamline)

    def testEXFELSPBKBBeamlinePropagation(self):
        """ Testing propagation through the EUXFEL SPB-SFX KB beamline. """

        # Setup output and cleanup.
        output_file = "prop_out.h5"
        self.__files_to_remove.append(output_file)

        # Propagate.
        propagate_s2e.propagate(self.__fel_source, output_file, exfel_spb_kb_beamline.get_beamline )

        # Check that output was generated.
        self.assertIn( output_file, os.listdir(".") )

def setupTestWavefront():
    """ Utility to setup a Gaussian wavefront. Geometry corresponds to SPB-SFX Day1 configuration. """

    ### Geometry ###
    src_to_hom1 = 257.8 # Distance source to HOM 1 [m]
    src_to_hom2 = 267.8 # Distance source to HOM 2 [m]
    src_to_crl = 887.8  # Distance source to CRL [m]
    src_to_exp = 920.42 # Distance source to experiment [m]

    # Central photon energy.
    ekev = 8.4 # Energy [keV]

    # Pulse parameters.
    qnC = 0.5               # e-bunch charge, [nC]
    pulse_duration = 9.e-15 # [s]
    pulseEnergy = 1.5e-3    # total pulse energy, J

    # Coherence time
    coh_time = 0.24e-15     # [s]

    # Distance to first HOM.
    z1 = src_to_hom1

    # Angular distribution
    theta_fwhm = 2.124e-6 # Beam divergence        # From Patrick's raytrace.

    wlambda = 12.4*1e-10/ekev # wavelength [AKM]
    w0 = wlambda/(numpy.pi*theta_fwhm) # beam waist
    zR = (numpy.pi*w0**2)/wlambda #Rayleigh range
    fwhm_at_zR = theta_fwhm*zR #FWHM at Rayleigh range
    sigmaAmp = w0/(2*numpy.sqrt(numpy.log(2))) #sigma of amplitude

    # Number of points in each x and y dimension.
    np=100

    bSaved=False
    dx = 10.e-6;
    range_xy = dx*(np-1)
    #print ('range_xy = ', range_xy)
    nslices = 10;

    srwl_wf = build_gauss_wavefront(np, np, nslices, ekev, -range_xy/2, range_xy/2,
             -range_xy/2, range_xy/2 ,coh_time/numpy.sqrt(2),
             sigmaAmp, sigmaAmp, src_to_hom1,
             pulseEn=pulseEnergy, pulseRange=8.)

    wf = Wavefront(srwl_wf)

    wf.store_hdf5("source.h5")

if __name__ == '__main__':
    unittest.main()

