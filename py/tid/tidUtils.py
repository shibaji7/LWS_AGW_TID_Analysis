#!/usr/bin/env python

"""tidUtils.py: utility module to support other functions."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import configparser
import os
import h5py

import numpy as np


def get_config(key, section="LWS"):
    config = configparser.ConfigParser()
    config.read("config/conf.ini")
    val = config[section][key]
    return val


def get_gridded_parameters(
    q, xparam="beam", yparam="slist", zparam="v", r=0, rounding=True
):
    """
    Method converts scans to "beam" and "slist" or gate
    """
    plotParamDF = q[[xparam, yparam, zparam]]
    if rounding:
        plotParamDF.loc[:, xparam] = np.round(plotParamDF[xparam].tolist(), r)
        plotParamDF.loc[:, yparam] = np.round(plotParamDF[yparam].tolist(), r)
    plotParamDF = plotParamDF.groupby([xparam, yparam]).mean().reset_index()
    plotParamDF = plotParamDF[[xparam, yparam, zparam]].pivot(xparam, yparam)
    x = plotParamDF.index.values
    y = plotParamDF.columns.levels[1].values
    X, Y = np.meshgrid(x, y)
    # Mask the nan values! pcolormesh can't handle them well!
    Z = np.ma.masked_where(
        np.isnan(plotParamDF[zparam].values), plotParamDF[zparam].values
    )
    return X, Y, Z


def get_folder(date):
    """
    Get folder by date
    """
    fold = os.path.join(get_config("FOLDER"), date.strftime("%Y-%m-%d"))
    os.makedirs(fold, exist_ok=True)
    return fold

def to_date(ts):
    """
    Convert to date
    """
    ts = dt.datetime.fromordinal(int(ts)) + dt.timedelta(days=ts%1)-dt.timedelta(days=366)
    t = ts.replace(second=0, microsecond=0) + dt.timedelta(seconds=30*np.rint(ts.second/30))
    return t

def read_tec_mat_files(fname):
    """
    Read TEC files
    """
    tec = h5py.File(fname)
    times = np.concatenate(tec["UTT"],axis=0)
    times = [to_date(ts) for ts in times]
    return tec, times

def fetch_tec_by_datetime(ts, mat, times):
    t = ts.replace(second=0, microsecond=0) + dt.timedelta(seconds=30*np.rint(ts.second/30))
    idx = times.index(t)
    dset = mat[mat["fulltimedata"][idx][0]]
    ipplat, ipplon, dtec = dset[2], dset[3], dset[5]
    return ipplat, ipplon, dtec