"""
This module contains functions for finding the pressure, temperature,
density and total pressure ratios across a normal shock for a calorically
perfect gas.

.. todo:: Need to add total pressure ratio.
"""

def pressure_ratio(mach, gamma=1.4):
    """Computes the pressure ratio across a normal shock.

    :param mach: Mach number upstream of the shock.
    :param gamma: ratio of specific heats.
    :returns: pressure ratio across a normal shock.
    """
    return 1.0 + 2.0*gamma/(gamma+1.0)*(mach**2 - 1.0)

def density_ratio(mach, gamma=1.4):
    """Computes the density ratio across a normal shock.

    :param mach: Mach number upstream of the shock.
    :param gamma: ratio of specific heats.
    :returns: density ratio across a normal shock.
    """
    return (gamma+1.0)*mach**2/((gamma-1.0)*mach**2+2.0)

def temperature_ratio(mach, gamma=1.4):
    """Computes the pressure ratio across a normal shock.

    :param mach: Mach number upstream of the shock.
    :param gamma: ratio of specific heats.
    :returns: temperature ratio across a normal shock.
    """
    return 1.0 + 2.0*(gamma-1.0)/(gamma+1.0)**2*(gamma*mach**2 + 1.0)/mach**2*(mach**2-1.0)

def mach_post_shock(mach, gamma=1.4):
    """Computes the Mach number after the shock.

    :param mach: Mach number upstream of the shock.
    :param gamma: ratio of specific heats.
    :returns: Mach number downstream of the shock.
    """
    return ((1.0 + (gamma-1.0)/2.0*mach**2)/(gamma*mach**2 - (gamma-1.0)/2.0))**(0.5)
