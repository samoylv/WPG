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


class ScreenTest(unittest.TestCase):
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
        """ Testing the construction of a screen instance."""
        screen = optical_elements.Screen()

        # Check inheritance.
        self.assertIsInstance(screen, optical_elements.Empty)
        self.assertIsInstance(screen, optical_elements.WPGOpticalElement)
        self.assertIsInstance(screen, object)

        # Check default handling of parameters.
        self.assertEqual(screen._Screen__filename, os.path.abspath('screen.h5'))

    def testShapedConstruction(self):
        """ Check construction of a Screen with parameters. """
        screen = optical_elements.Screen(filename="wf_at_screen.h5")

        # Check inheritance.
        self.assertIsInstance(screen, optical_elements.Empty)
        self.assertIsInstance(screen, optical_elements.WPGOpticalElement)
        self.assertIsInstance(screen, object)

        # Check default handling of parameters.
        self.assertEqual(screen._Screen__filename, os.path.abspath('wf_at_screen.h5'))

    def testMisshapedConstruction(self):
        """ Check construction of a Screen with faulty parameters. """

        self.assertRaises(TypeError, optical_elements.Screen, filename=['not','a', 'string'])
        self.assertRaises(TypeError, optical_elements.Screen, filename=1.436)
        self.assertRaises(TypeError, optical_elements.Screen, filename=1)
        self.assertRaises(IOError, optical_elements.Screen, filename=__file__)
        self.assertRaises(IOError, optical_elements.Screen, filename='/not/an/existing/directory/screen.h5')

if __name__ == '__main__':
    unittest.main()

