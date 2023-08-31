#!/usr/bin/env python

"""wave.py: module is dedicated to run wavelete analysis runs"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import numpy as np
import pywt

def get_families(short=False):
    """
    Get all wavelete families in PyWavelets
    """
    return pywt.families(short)

def get_wavelist(kind="all", family="morl"):
    """
    Get all wavelists
    """
    return pywt.wavelist(family, kind)

def cwt_wavelete_transform(time, signal, waveletname="morl", scales=np.arange(1, 300), dt=60.0):
    """
    Conduct wavelete transform
    
    Parameters:
    -----------
    
    signal: Y-signal (array)
    waveletname: Name of the wavelete (str)
    scales: Scales of the CWT (array)
    dt: Timeseries time cadence in sec
    """
    [coefficients, frequencies] = pywt.cwt(signal, scales, waveletname, dt)
    power = (abs(coefficients)) ** 2
    period = 1. / (frequencies*dt)
    fdata = {
        "time": time,
        "waveletname": waveletname,
        "scales": scales,
        "dt": dt,
        "coefficients": coefficients,
        "frequencies": frequencies,
        "power": power,
        "period": period
    }
    return fdata


import matplotlib.pyplot as plt
def plot_wavelet(
    fdata, waveletname = 'morl', 
    cmap = plt.cm.seismic, 
    title = 'Wavelet Transform (Power Spectrum) of signal', 
    ylabel = 'Period (Minutes)', 
    xlabel = 'Time'
):
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel(xlabel, fontsize=18)
    im = ax.pcolor(fdata["time"], fdata["period"], fdata["power"], cmap=cmap)
    cbar_ax = fig.add_axes([0.95, 0.5, 0.03, 0.25])
    cbar = fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    cbar.set_label("Power")
    return