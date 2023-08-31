#!/usr/bin/env python

"""generate.py: module is dedicated to generate timeserise data for analysis"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
import pandas as pd

def create_synthetic_data(Amp, Phi, T, t, cad, stn="mit"):
    """
    This function is responsible for creating a
    synthetic magnetic field with following sets
    of parameters.
    
    Parameter:
    ----------
    Amp (list) - Magntitue at different freuency components (m)
    Phi (list) - Phase at different freuency components (m)
    T (list) - Periods of different freuency components (m)
    t (float) - Total time (in seconds)
    """
    data = {}
    t = np.linspace(0, t, t, endpoint=False)
    f = np.zeros(len(t))
    for A, Phi, T in zip(Amp, Phi, T):
        f += A * np.sin(2 * np.pi * t / T + np.deg2rad(Phi))
    data[stn] = pd.DataFrame()
    data[stn]["dTEC"], data[stn]["time"], data[stn]["dt"] = f, t, cad
    return data