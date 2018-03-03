# The below is to ignore numpy errors from pylint:
#     pylint: disable=E1101

"""Functions for computations related to Prandtl-Meyer expansions.

This module contains functions for Prandtl-Meyer expansions.  It
includes a function for computing the Prandtl-Meyer function, nu,
and another for determining the Mach number which corresponds to
a certain value of the Prandtl-Meyer function, nu.  Because the
equation is nonlinear, this is accomplished using an optimization
function from SciPy.
"""

import numpy as np
from scipy.optimize import minimize

def pmfunction(mach, gamma=1.4):
    """Prandtl-Meyer function.

    :param mach: Mach number.
    :param gamma: ratio of specific heats.
    :returns: Prandtl-Meyer function for given Mach number.
    """
    return np.sqrt((gamma+1.0)/(gamma-1.0))*\
            np.arctan(np.sqrt((gamma-1.0)/(gamma+1.0)*(mach**2-1.0))) -\
            np.arctan(np.sqrt(mach**2 - 1.0))

def compute_mach(nu_value, gamma=1.4, machguess=2.5, mach_bounds=(1.0, 50.0)):
    """Computes the Mach number which corresponds to a given value of the PM function.

    Function for computing the Mach number which corresponds to a given
    value of the Prandtl-Meyer function.

    :param nu_value: value of the Prandtl-Meyer function.
    :param gamma: ratio of specific heats.
    :param machguess: initial guess for Mach number.
    :param mach_bounds: bounds used in iterative procedure.
    :returns: Mach number which corresponds to the given Prandtl-Meyer function.
    """
    res = minimize(lambda m: (nu_value - pmfunction(m, gamma=gamma))**2/nu_value**2,
                   machguess, bounds=[mach_bounds])
    return res.x[0]
