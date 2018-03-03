# The below is to ignore numpy errors from pylint:
#     pylint: disable=E1101

"""This module holds a class which represents a wall profile.

.. module:: wallprofile
.. moduleauthor:: James Grisham
"""

import numpy as np

def _interp(x0, y0, x1, y1, x):
    """interp is a simple linear interpolation function.

    The interp function takes two points defined in the x-y plane along
    with an x-coordinate and returns the value of y at the given
    x-coordinate using a linear interpolation.  Also works for
    extrapolation, but should be used with care.  It is only meant for
    local use inside this module.

    :param x0: first x-coordinate of data to be interpolated.
    :param y0: first y-coordinate of data to be interpolated.
    :param x1: second x-coordinate of data to be interpolated.
    :param y1: second y-coordinate of data to be interpolated.
    :param x: x-coordinate at which y-value is desired.
    :returns: a float which is the interpolated value.
    """
    slope = (y1 - y0)/(x1 - x0)
    yintercept = y0 - slope*x0
    return slope*x + yintercept

class WallProfile:
    """
    This class represents a profile of pressure and skin friction
    along a solid surface.
    """

    def __init__(self, x, pressure, tau):
        """ Constructor for `WallProfile` class.

        :param x: streamwise coordinate.
        :param pressure: wall static pressure.
        :param tau: wall shear stress.
        """
        self.xcoord = x
        self.pressure = pressure
        self.shear_stress = tau
        self.sep_reattach_points = None

    def is_separated(self):
        """ Checks for boundary layer separation.

        This method checks for sign changes in cf.  If cf does
        change sign at some point, it returns true, which indicates
        that the boundary layer does separate at some point along the
        given surface.

        :returns: boolean value indicating whether or not the boundary layer separates.
        """
        zero_indices = np.where(np.diff(np.sign(self.shear_stress)))[0]
        if zero_indices.size == 0:
            return False
        else:
            return True

    def compute_zeros(self):
        """ Computes locations at which separation/reattachment of the BL occur.

        This method uses the wall shear stress profile to determine where the
        boundary layer separates/reattaches.  This is accomplished by determining
        the streamwise locations at which the wall shear stress changes sign.

        :returns: list of floats at which separation/reattachment occur.
        """

        # Finding indices of locations at which tau = 0
        zero_indices = np.where(np.diff(np.sign(self.shear_stress)))[0]
        if zero_indices.size == 0:
            return None

        # Need to use quadratic interpolation here.
        streamwise_locations = [_interp(self.tau[i], self.xcoord[i], self.tau[i+1], \
                self.xcoord[i+1], 0.0) for i in zero_indices]
        self.sep_reattach_points = streamwise_locations
        return streamwise_locations
