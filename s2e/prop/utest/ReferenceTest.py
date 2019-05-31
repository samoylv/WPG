##########################################################################
#                                                                        #
# Copyright (C) 2017 Carsten Fortmann-Grote                              #
# Contact: Carsten Fortmann-Grote <carsten.grote@xfel.eu>                #
#                                                                        #
# This file is part of WPG/s2e/prop.                                     #
# See WPG/LICENSE for licensing information.                             #
#                                                                        #
##########################################################################

""" Test module for the reference tests.

    :author : Carsten Fortmann-Grote
    :institution : European XFEL GmbH, Holzkoppel 4, 22869 Schenefeld, Germany
    :creation 20170411

"""
import os, sys, shutil
import h5py
import numpy
import requests
import unittest


# Import the class to test.
sys.path.insert(0,'../../..')


from s2e.prop import exfel_spb_kb_beamline
from s2e.prop import propagate_s2e

# Make sure WPG is in $PYTHONPATH.
from wpg import Beamline, Wavefront


class ReferenceTest(unittest.TestCase):
    """
    Test class for WPG/SRW reference test.
    """

    @classmethod
    def setUpClass(cls):
        """ Setting up the test class. """

        # Get FELsource data if not present locally.
        if "FELsource_5keV_9fs_slicing12_3-29fs_104x104.h5" not in os.listdir("."):
            try:
                cls.__fel_source = wgetData(url = "https://docs.xfel.eu/alfresco/d/a/workspace/SpacesStore/343db429-a2df-4ae9-ac2e-c289a34c6fc1/FELsource_5keV_9fs_slicing12_3-29fs_104x104.h5")
            except:
                raise RuntimeError("Unable to download test data. Please try again later. If problem persists, contact support.")
                sys.exit()
        else:
            cls.__fel_source = "FELsource_5keV_9fs_slicing12_3-29fs_104x104.h5"

        # Get reference data if not present locally.
        if "reference_prop_out.h5" not in os.listdir("."):
            try:
                cls.__reference_prop_out = wgetData(url="https://docs.xfel.eu/alfresco/d/a/workspace/SpacesStore/9b263892-1cb0-42e3-b05a-7507f66c67db/reference_prop_out.h5")
            except:
                raise RuntimeError("Unable to download test data. Please try again later. If problem persists, contact support.")
                sys.exit()
        else:
            cls.__reference_prop_out = "reference_prop_out.h5"

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

    def test9fsKB(self):
        """ Reference test with a 9fs 5 keV pulse in the SPB-SFX KB nanofocus beamline."""

        # Ensure cleanup.
        self.__files_to_remove.append("prop_out.h5")

        # Get beamline interface.
        bl = exfel_spb_kb_beamline.get_beamline

        # Propagate.
        propagate_s2e.propagate(self.__fel_source, "prop_out.h5", bl)

        # Read in horizontal field component.
        Ehor = h5py.File("prop_out.h5", 'r')['data/arrEhor'].value.flatten()

        # Read in reference data.
        Ehor_reference = h5py.File(self.__reference_prop_out, 'r')['data/arrEhor'].value.flatten()

        for e, e_ref in zip(Ehor, Ehor_reference):
            self.assertAlmostEqual(e, e_ref)


def wgetData(url=None):
    """ Download source data from alfresco repository. """

    # Local filename where data will be saved.
    local_filename = url.split('/')[-1]


    # Make https request.
    print "Attempting to download %s." % (url)
    r = requests.get(url, stream=True)

    # Write to local file in chunks of 1 MB.
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=1024):
            if chunk: # filter out keep-alive new chunks
                f.write(chunk)

    # After successful write, close the https connection.
    f.close()

    # Return.
    print "Download completed and saved to %s." % (local_filename)
    return local_filename

if __name__ == '__main__':
    unittest.main()


