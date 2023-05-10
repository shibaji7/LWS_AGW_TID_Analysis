#!/usr/bin/env python

"""tidFilters.py: filter module."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import copy
import datetime as dt
import music
import numpy as np
import pandas as pd


def filterBasedOnLimits(
    fd,
    dates,
    gateLimits,
    timeLimit,
    timeRes=120,
    numtaps=101,
    cutoffs_min=[15, 60],
    pfilt_thresh=1.0,
):
    """Interpolates the data in a musicArray object along the beams of the radar.  This method will ensure that no
    rangegates are missing data.  Ranges outside of metadata['gateLimits'] will be set to 0.
    The result is stored as a new musicDataObj in the given musicArray object.
    """
    cutoff_high, cutoff_low = (1.0 / (60 * cutoffs_min[0]), 1.0 / (60 * cutoffs_min[1]))
    dataObj = music.musicArray(fd, sTime=dates[0], eTime=dates[1], fovModel="GS")
    dataObj.get_data_sets()
    music.defineLimits(dataObj, gateLimits=gateLimits)
    new_times = music.filterTimes(
        timeLimit[0], timeLimit[1], timeRes=timeRes, numTaps=numtaps
    )
    music.defineLimits(dataObj, timeLimits=new_times)
    dataObj.active.applyLimits()
    music.beamInterpolation(dataObj)
    music.timeInterpolation(dataObj, timeRes=timeRes)
    music.determineRelativePosition(dataObj)
    filt = music.filter(
        dataObj, numtaps=numtaps, cutoff_low=cutoff_low, cutoff_high=cutoff_high
    )
    dataset = music.getDataSet(dataObj, "active")
    time, bmnum, slist, lat, lon, srange = [], [], [], [], [], []
    o_data = copy.copy(dataset.data)
    #o_data[np.abs(o_data) < pfilt_thresh] = np.nan
    data = []
    for i, t in enumerate(dataset.time):
        for j, b in enumerate(dataset.fov.beams):
            L = len(dataset.data[i, j, :])
            bmnum.extend(L * [b])
            time.extend(L * [t])
            lat.extend(dataset.fov.latCenter[b, :].tolist())
            lon.extend(dataset.fov.lonCenter[b, :].tolist())
            srange.extend(dataset.fov.slantRCenter[b, :].tolist())
            slist.extend(np.arange(gateLimits[0], gateLimits[1] + 1).tolist())
            data.extend(o_data[i, j, :].tolist())
    dat = pd.DataFrame()
    dat["time"], dat["bmnum"], dat["p_l"] = time, bmnum, data
    dat["glat"], dat["glon"], dat["srange"] = lat, lon, srange
    dat["slist"] = slist
    return dat


def fftBasedOnLimits(
    fd,
    beam,
    srange,
    timeLimit,
    timeRes=60,
    flim=[1e-6, 1e-2],
):
    """
    FFT based on limits
    """
    frame = fd.frame.copy()
    frame = frame[
        (frame.time>=timeLimit[0]-dt.timedelta(seconds=timeRes/2)) &
        (frame.time<timeLimit[1]+dt.timedelta(seconds=timeRes/2)) &
        (frame.bmnum==beam) &
        (frame.srange==srange)
    ]
    
    N = int((timeLimit[1]-timeLimit[0]).total_seconds()/timeRes)
    t = timeRes*np.arange(N)
    if N>len(frame):
        t = t[:len(frame)]
    if N<len(frame):
        t = timeRes*np.arange(len(frame))
    
    from astropy.timeseries import LombScargle
    frequency = np.linspace(
        flim[0], 
        flim[1], 
        int(10**(np.log10(flim[1])-np.log10(flim[0])))
    )
    ls = LombScargle(t, np.asarray(frame.p_l))
    o = pd.DataFrame()
    o["freq"], o["power"] = ls.autopower() #frequency, ls.power(frequency)
    return o
