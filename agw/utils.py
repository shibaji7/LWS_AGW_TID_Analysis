#!/usr/bin/env python

"""utils.py: utility module for data parsing and helping analsysis or plots."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

sys.path.extend(["agw/", "agw/tid/"])

import configparser
import datetime as dt
import glob
from types import SimpleNamespace

import numpy as np
import pandas as pd
from loguru import logger


class RecursiveNamespace(SimpleNamespace):
    """ """

    @staticmethod
    def map_entry(entry):
        if isinstance(entry, dict):
            return RecursiveNamespace(**entry)
        return entry

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        for key, val in kwargs.items():
            if type(val) == dict:
                setattr(self, key, RecursiveNamespace(**val))
            elif type(val) == list:
                setattr(self, key, list(map(self.map_entry, val)))


def get_config(key, section="LWS"):
    config = configparser.ConfigParser()
    config.read("config/conf.ini")
    val = config[section][key]
    return val


def read_logs(_filestr="config/record.logs/*.log"):
    """
    Read all log files and store them each records
    to a list of dictionary
    """
    files = sorted(glob.glob(_filestr))
    _modes = pd.DataFrame()
    for f in files:
        logger.info(f"Log file {f}")
        _modes = pd.concat([_modes, pd.read_csv(f, parse_dates=["date"])])
    _modes.ut_start = [
        r["date"] + dt.timedelta(hours=r["ut_start"]) for i, r in _modes.iterrows()
    ]
    _modes.ut_end = [
        r["date"] + dt.timedelta(hours=r["ut_end"]) for i, r in _modes.iterrows()
    ]
    logger.info(f"Total _mode log entry {len(_modes)}")
    return _modes


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
    else:
        plotParamDF[xparam] = plotParamDF[xparam].tolist()
        plotParamDF[yparam] = plotParamDF[yparam].tolist()
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
