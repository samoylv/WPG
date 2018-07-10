""" Test module for wpg utilities

    :author : Carsten Fortmann-Grote
    :institution : European XFEL GmbH, Holzkoppel 4, 22869 Schenefeld, Germany
    :creation date: 20170329

"""

import numpy
import os
import shutil
import sys
import unittest

# Add wpg to global namespace.
sys.path.insert(0,os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)),'..','..')))


from wpg import Wavefront

from wpg.generators import build_gauss_wavefront

# Import the module to test.
from wpg import wpg_uti_wf as utilities

def setup_gauss_wavefront(sanity=True):
    """ Gaussian wavefront builder """
    if sanity:
        np = 700
        nslices = 100
    else:
        np = 10
        nslices=10

    photon_energy = 1.6e3
    theta_fwhm = 6.0e-6 ### ADJUST ME

    wlambda = 1.24*1e-6/photon_energy # wavelength [AKM]
    w0 = wlambda/(numpy.pi*theta_fwhm) # beam waist
    zR = (numpy.pi*w0**2)/wlambda #Rayleigh range
    sigmaAmp = w0/(2*numpy.sqrt(numpy.log(2))) #sigma of amplitude
    src_to_aperture_distance = 170.0
    pulse_energy = 4e-4 # [J]

    # Coherence time
    pulse_duration = 30.0e-15 # [s]
    coh_time = pulse_duration/10.0     # estimate, [s]


    # expected beam radius at HOM1 position to get the range of the wavefront
    range_xy = w0*numpy.sqrt(1+(src_to_aperture_distance/zR)**2) *5.5

    srwl_wf = build_gauss_wavefront(np, np, nslices, photon_energy/1e3, -range_xy/2, range_xy/2,
                                    -range_xy/2, range_xy/2 ,coh_time/numpy.sqrt(2),
                                    sigmaAmp, sigmaAmp, src_to_aperture_distance,
                                    pulseEn=pulse_energy, pulseRange=8.)

    return Wavefront(srwl_wf)

class WavefrontUtilitiesTest(unittest.TestCase):
    """
    Test class for the wpg_uti_wf module.
    """

    @classmethod
    def setUpClass(cls):
        """ Setting up the test class. """
        pass

    @classmethod
    def tearDownClass(cls):
        """ Tearing down the test class. """
        pass

    def setUp(self):
        """ Setting up a test. """
        self.__files_to_remove = []
        self.__dirs_to_remove = []

        self.__wavefront_sane = setup_gauss_wavefront(sanity=True)
        self.__wavefront_insane = setup_gauss_wavefront(sanity=False)


    def tearDown(self):
        """ Tearing down a test. """

        for f in self.__files_to_remove:
            if os.path.isfile(f): os.remove(f)
        for d in self.__dirs_to_remove:
            if os.path.isdir(d): shutil.rmtree(d)

        del self.__wavefront_sane
        del self.__wavefront_insane

        self.__wavefront_sane = None
        self.__wavefront_insane = None

    def testCheckSampling(self):
        """ Test the check_sampling utility. """

        report_sane = utilities.check_sampling(self.__wavefront_sane)
        report_insane = utilities.check_sampling(self.__wavefront_insane)


        self.assertIn("Horizontal Fresnel zone extension within", report_sane)
        self.assertIn("Vertical Fresnel zone extension within", report_sane)
        self.assertIn("Horizontal Fresnel zone extension NOT within", report_insane)
        self.assertIn("Vertical Fresnel zone extension NOT within", report_insane)


if __name__ == '__main__':
    unittest.main()

