""" Test module for the Converters suite.
"""
import os, sys, shutil
import numpy
import unittest


# Import the class to test.
sys.path.insert(0,'../../..')


from wpg import Wavefront
from wpg.converters import genesis_v2
from wpg.wpg_uti_wf import plot_intensity_map

class ConvertersTest(unittest.TestCase):
    """
    Test class for the converter modules.
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

    def tearDown(self):
        """ Tearing down a test. """

        for f in self.__files_to_remove:
            if os.path.isfile(f): os.remove(f)
        for d in self.__dirs_to_remove:
            if os.path.isdir(d): shutil.rmtree(d)

    def testGenesisV2(self):
        """ Testing the conversion of genesis1.3 v2 output into a wpg wavefront."""

        # Convert.
        wavefront = genesis_v2.read_genesis_file("lcls/lcls.out", "lcls/lcls.out.dfl")

        # Check type.
        self.assertIsInstance(wavefront, Wavefront)

        # Plot.
        plot_intensity_map(wavefront)

if __name__ == '__main__':
    unittest.main()

