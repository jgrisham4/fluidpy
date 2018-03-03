"""This module contains utilities for dealing with boundary layer profiles.

..module:: Module for boundary layer profiles.
..moduleauthor:: James Grisham
..todo:: Add y+, u+, VanDriest transformation, etc.

"""

from scipy.integrate import simps
from scipy.interpolate import interp1d

class BLProfile:
    """
    This class represents a boundary layer profile.
    """

    def __init__(self, normal_coord, streamwise_velocity, edge_velocity,
                 normal_velocity=None, density=None, edge_density=None):
        """Constructor for BLProfile class.

        ..note:: This module assumes that numpy arrays are passed in for
            things like the normal coordinate, and streamwise velocity, etc.

        :param normal_coord: normal coordinate.
        :param streamwise_velocity: streamwise velocity.
        :param edge_velocity: velocity at the edge of the boundary layer.
        :param normal_velocity: normal velocity.
        :param density: optional input for a compressible profile.
        :param edge_density: density at the edge of the boundary layer.
        :type normal_coord: numpy.ndarray.
        :type streamwise_velocity: numpy.ndarray.
        :type edge_velocity: float.
        :type normal_velocity: numpy.ndarray or None.
        :type density: numpy.ndarray or None.
        :type edge_density: numpy.ndarray or None.
        """
        self.normal_coord = normal_coord
        self.streamwise_velocity = streamwise_velocity
        self.normal_velocity = normal_velocity
        self.edge_velocity = edge_velocity
        self.density = density
        self.edge_density = edge_density
        self.uque = streamwise_velocity/edge_velocity
        self.delta_0 = None
        self.delta_1 = None
        self.delta_2 = None
        self.delta_1i = None
        self.delta_2i = None
        self.shape_factor = None
        self.shape_factor_i = None

    def compute_thickness(self, uque_threshold=0.99):
        """Computes the thickness of the boundary layer.

        This method computes the thickness of the boundary layer by interpolating
        for the normal distance required for :math:`u/U_e \\approx 0.99`.

        :param uque_threshold: value used to interpolate for edge of boundary layer.
        :returns: the thickness of the boundary layer.
        """

        # Constructing fit
        fit = interp1d(self.uque, self.normal_coord)

        # Trying to interpolate, catching error if profile isn't large enough
        try:
            self.delta_0 = fit(uque_threshold)
        except ValueError:
            print("Profile isn't large enough to find boundary layer thickness.")
            raise

        return self.delta_0

    def compute_integral_thicknesses(self, integrate_to_inf=False):
        """Computes displacement and momentum thicknesses for the given profile.

        This method computes the displacement and momentum thicknesses for the
        given profile using numerical integration via SciPy's Simpson's rule
        function.

        :param integrate_to_inf: used to set upper bound on integration.

        ..note:: The integrate_to_inf argument sets the upper bound on the integration
        to be the boundary layer thickness, :math:`\\delta`, if false, otherwise the
        integral is accomplished over the entire profile.
        """

        # Finding indices to truncate arrays if necessary
        if not integrate_to_inf:
            bl_indices = [i for i, yc in enumerate(self.normal_coord) if yc <= self.delta_0]
        else:
            bl_indices = [i for i in range(0, len(self.normal_coord))]

        # Computing incompressible integral thicknesses
        normal_coord = self.normal_coord[bl_indices]
        streamwise_velocity = self.streamwise_velocity[bl_indices]
        self.delta_1i = simps(1.0 - streamwise_velocity/self.edge_velocity, normal_coord)
        self.delta_2i = simps((1.0 - streamwise_velocity/self.edge_velocity)*\
                streamwise_velocity/self.edge_velocity, normal_coord)
        self.shape_factor_i = self.delta_1i/self.delta_2i

        # Checking to see if density profile is provided
        # If so, compressible integral thicknesses are calculated.
        if self.density is not None:
            if self.edge_density is None:
                self.edge_density = self.density[-1]
            density = self.density[bl_indices]
            self.delta_1 = simps(1.0 - density/self.edge_density*\
                    streamwise_velocity/self.edge_velocity, normal_coord)
            self.delta_2 = simps(density*streamwise_velocity/\
                    (self.edge_density*self.edge_velocity)*\
                    (1.0 - streamwise_velocity/self.edge_velocity), normal_coord)
            self.shape_factor = self.delta_1/self.delta_2
