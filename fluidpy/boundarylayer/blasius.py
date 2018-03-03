"""Contains functions for computing the Blasius solution.
"""

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import minimize

def _deriv(ycoord, eta):
    """Derivative function.
    """
    return np.array([ycoord[1], ycoord[2], -ycoord[0]*ycoord[2]])

def _solve_ode(fppp0):
    """Function for solving ODE.

    :param fppp0: initial condition on third first-order system.
    """
    eta = np.linspace(0.0, 6.0, 1000)
    ycoord0 = np.array([0.0, 0.0, fppp0])
    ycoord = odeint(_deriv, ycoord0, eta)
    return (ycoord[-1, 1]-1.0)**2

def compute_blasius_solution():
    """Computes Blasius solution.

    This function computes the Blasius solution by numerically integrating
    the system of ordinary differential equations along with optimization
    to satisfy the boundary conditions.
    :returns: a tuple which contains (eta, f').
    """

    # Computing proper initial condition using optimization
    optresult = minimize(_solve_ode, 2.0)

    # Solving the system of ODEs again using the proper initial condition
    eta = np.linspace(0.0, 4.5, 1000)
    ycoord0 = np.array([0.0, 0.0, optresult.x])
    ycoord = odeint(_deriv, ycoord0, eta)
    fprime = ycoord[:, 1]

    return (eta, fprime)

def boundary_layer_thickness(streamwise_location, reynolds_x):
    """Computes the boundary layer thickness for incompressible flow over a flat plate.

    This function uses the result from Blasius' solution to compute
    the boundary layer thickness for steady flow over a flat plate.

    :param streamwise_location: corresponds to x for flow over a flat plate.
    :param reynolds_x: Reynolds number based on the given streamwise location.
    :returns: boundary layer thickness.
    """
    return 5.0*streamwise_location/np.sqrt(reynolds_x)

def displacement_thickness(streamwise_location, reynolds_x):
    """Computes the incompressible displacement thickness for flat plate.

    This function uses the result from Blasius' solution to compute
    the incompressible displacement thickness for steady flow over a
    flat plate.

    :param streamwise_location: corresponds to x for flow over a flat plate.
    :param reynolds_x: Reynolds number based on the given streamwise location.
    :returns: displacement thickness.
    """
    return 1.721*streamwise_location/np.sqrt(reynolds_x)

def momentum_thickness(streamwise_location, reynolds_x):
    """Computes the incompressible momentum thickness for flat plate.

    This function uses the result from Blasius' solution to compute
    the incompressible momentum thickness for steady flow over a
    flat plate.

    :param streamwise_location: corresponds to x for flow over a flat plate.
    :param reynolds_x: Reynolds number based on the given streamwise location.
    :returns: momentum thickness.
    """
    return 0.664*streamwise_location/np.sqrt(reynolds_x)

def skin_friction(reynolds_x):
    """Computes skin friction for incompressible flow over a flat plate.

    This function uses the result from Blasius' solution to compute the
    skin friction coefficient for incompressible, constant property flow
    over a flat plate.

    :param reynolds_x: Reynolds number based on the given streamwise location.
    :returns: skin friction coefficient.
    """
    return 0.664/np.sqrt(reynolds_x)
