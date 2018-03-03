"""
This module contains functions for computing the pressure, density, temperature
and total pressure ratios across an oblique shock.

.. todo:: Need to add total pressure ratio.
"""

import math

def pressure_ratio(mach, theta, gamma=1.4):
    """Computes the pressure ratio across a oblique shock.

    :param mach: Mach number.
    :param theta: deflection angle in radians.
    :param gamma: ratio of specific heats.
    :returns: pressure ratio across the oblique shock.
    """
    beta_tmp = wave_angle(mach, theta, gamma=gamma)
    return 1.0 + 2.0*gamma/(gamma+1.0)*(mach**2*math.sin(beta_tmp)**2 - 1.0)

def density_ratio(mach, theta, gamma=1.4):
    """Computes the density ratio across a oblique shock.

    :param mach: Mach number.
    :param theta: deflection angle in radians.
    :param gamma: ratio of specific heats.
    :returns: density ratio across the oblique shock.
    """
    beta_tmp = wave_angle(mach, theta, gamma=gamma)
    return (gamma+1.0)*mach**2*math.sin(beta_tmp)**2/(
        (gamma-1.0)*mach**2*math.sin(beta_tmp)**2 + 2.0)

def temperature_ratio(mach, theta, gamma=1.4):
    """Computes the temperature ratio across a oblique shock.

    :param mach: Mach number.
    :param theta: deflection angle in radians.
    :param gamma: ratio of specific heats.
    :returns: temperature ratio across the oblique shock.
    """
    beta_tmp = wave_angle(mach, theta, gamma=gamma)
    return 1.0 + 2.0*(gamma-1.0)/(gamma+1.0)**2*(mach**2*math.sin(beta_tmp)**2 - 1.0)/(
        mach**2*math.sin(beta_tmp)**2)*(gamma*mach**2*math.sin(beta_tmp)**2 + 1.0)

def wave_angle(mach, theta, gamma=1.4):
    """Computes the wave angle.

    This function computes the wave angle when provided with the Mach number
    and the deflection angle.

    :param mach: Mach number.
    :param theta: deflection angle in radians.
    :param gamma: ratio of specific heats.
    :returns: returns the wave angle in radians.
    """

    # Assuming theta is in radians
    delta = 1.0  # weak soln, delta = 0 is strong solution

    # Equations
    llambda = ((mach**2 - 1.0)**2 - 3.0*(1.0 + (gamma - 1.0)/2.0*mach**2)*
               (1.0 + (gamma+1.0)/2.0*mach**2)*(math.tan(theta))**2)**(1.0/2.0)

    chi = ((mach**2 - 1)**3 - 9*(1 + (gamma-1)/2.0*mach**2) \
            *(1 + (gamma-1)/2.0*mach**2 +
              (gamma+1)/4*mach**4)*(math.tan(theta))**2)/llambda**3
    numerator = mach**2 - 1 + 2.0*llambda*math.cos((4*math.pi*delta + math.acos(chi))/3.0)
    denominator = 3.0*(1 + (gamma-1)/2.0*mach**2)*math.tan(theta)

    return math.atan2(numerator, denominator)
