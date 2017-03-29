##########################################################################
#                                                                        #
# Copyright (C) 2017 Carsten Fortmann-Grote                              #
# Contact: Carsten Fortmann-Grote <carsten.grote@xfel.eu>                #
#                                                                        #
# This file is part of wpg/s2e/prop.                                     #
# It is free software: you can redistribute it and/or modify             #
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

import unittest
import sys

# Import suites to run.
from BeamlinesTest import BeamlinesTest

# Define the encapsulating test suite.
def suite():
    suites = [
            unittest.makeSuite(BeamlinesTest, 'test'),
             ]

    return unittest.TestSuite(suites)

# Run the top level suite and return a success status code. This enables running an automated git-bisect.
if __name__=="__main__":

    result = unittest.TextTestRunner(verbosity=2).run(suite())

    if result.wasSuccessful():
        print '---> OK <---'
        sys.exit(0)

    sys.exit(1)
