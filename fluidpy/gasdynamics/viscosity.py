"""This module contains a function for computing the viscosity of air.
"""

def viscosity(temperature):
    """Computes viscosity for air as a function of temperature.

    This function computes viscosity as a function of temperature
    using Sutherland's law.

    :param temperature: temperature in Kelvin.
    :returns: viscosity in kg/(m.s)
    """
    tref = 273.15
    muref = 1.716e-5
    sconst = 110.4
    return muref*(temperature/tref)**(3/2) * (tref + sconst)/(temperature + sconst)
