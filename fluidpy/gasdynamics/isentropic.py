"""This module holds isentropic relations for temperature
and pressure.
"""

def t0qt(mach, gamma=1.4):
    """Ratio of stagnation to static temperature.

    :param mach: Mach number.
    :param gamma: ratio of specific heats.
    :returns: ratio of total to static temperature.
    :rtype: float
    """
    return 1.0 + (gamma-1.0)/2.0*mach**2

def p0qp(mach, gamma=1.4):
    """Ratio of stagnation to static pressure.

    :param mach: Mach number.
    :param gamma: ratio of specific heats.
    :returns: ratio of total to static pressure.
    :rtype: float
    """
    return t0qt(mach, gamma=gamma)**(gamma/(gamma-1.0))
