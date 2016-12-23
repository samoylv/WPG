from wpg.useful_code.srwutils import AuxTransmAddSurfHeightProfileScaled as profile_tool
import numpy as np


def defineOPD(opTrErMirr,
              mdatafile,
              ncol,
              delim,
              orient,
              theta,
              scale=1.,
              stretching=1.):
    """
    Define optical path difference (OPD) from mirror profile

    :param opTrErMirr: The struct with wave front distortions from
                       mirror susrface errors
    :param mdatafile: An ascii file with mirror profile data
    :param ncol: Number of columns in the mirror profile file
    :param delim: Delimiter between numbers in an row, can be space (' '), tab '\t', etc.
    :param orient: Mirror orientation, 'x' (horizontal) or 'y' (vertical)
    :param theta: Mirror incidence angle
    :param scale: Scaling factor for the mirror profile height errors
    :param stretching: Scaling factor for the mirror profile x-axis (a hack, should be removed ASAP)
    """

    heightProfData = np.loadtxt(mdatafile).T
    heightProfData[0, :] = heightProfData[0, :] * stretching
    profile_tool(opTrErMirr, heightProfData, orient, theta, scale)
