"""
This module has a function for determining the thermodynamic
state at a given altitude using the standard atmosphere.
"""

import os
import sys
import numpy as np

def standard_atmosphere(altitude):
    """
    standard_atmosphere(altitude) computes the pressure and temperature
    at the given altitude in meters above sea level using the international standard atmosphere.
    @param altitude either a float or a numpy.ndarray which contains the altitude in meters.
    @return a tuple which contains the pressure and temperature.
    """

    # This needs to be re-written

