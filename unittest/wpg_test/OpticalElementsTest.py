""" Test module for wpg.OpticalElements.

    :author : Carsten Fortmann-Grote
    :institution : European XFEL GmbH, Holzkoppel 4, 22869 Schenefeld, Germany
    :creation date: 20170329

"""
import os, sys, shutil
import unittest


# Add sources to global namespace.
sys.path.insert(0,os.path.abspath(os.path.join(os.path.abspath(os.path.dirname(__file__)),'..','..')))

# Import the class to test.
from wpg import optical_elements


class EmptyTest(unittest.TestCase):
    """
    Test class for the screen class.
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

    def testDefaultConstruction(self):
        """ Testing the construction of an empty element."""
        empty = optical_elements.Empty()

        # Check inheritance.
        self.assertIsInstance(empty, optical_elements.Empty)
        self.assertIsInstance(empty, optical_elements.WPGOpticalElement)
        self.assertIsInstance(empty, object)

if __name__ == '__main__':
    unittest.main()

