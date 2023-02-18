#!/usr/bin/env python

"""dataUtils.py: utility module to for all level data"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt

import numpy as np
from loguru import logger


class Gate(object):
    """Class object to hold each range cell value"""

    def __init__(self, bm, i, params=["v", "w_l", "gflg", "p_l", "v_e"]):
        """
        initialize the parameters which will be stored
        bm: beam object
        i: index to store
        params: parameters to store
        """
        for p in params:
            if len(getattr(bm, p)) > i:
                setattr(self, p, getattr(bm, p)[i])
            else:
                setattr(self, p, np.nan)
        return


class Beam(object):
    """Class to hold one beam object"""

    def __init__(self, nrange_scatter):
        """initialize the instance"""
        self.nrange_scatter = nrange_scatter
        return

    def set(
        self,
        time,
        d,
        s_params=["bmnum", "noise.sky", "tfreq", "scan", "nrang"],
        v_params=["v", "w_l", "gflg", "p_l", "slist", "v_e"],
        k=None,
    ):
        """
        Set all parameters
        time: datetime of beam
        d: data dict for other parameters
        s_param: other scalar params
        v_params: other list params
        """
        for p in s_params:
            if p in d.keys():
                if p == "scan" and d[p] != 0:
                    setattr(self, p, 1)
                else:
                    setattr(self, p, d[p]) if k is None else setattr(self, p, d[p][k])
            else:
                setattr(self, p, None)
        for p in v_params:
            if p in d.keys():
                setattr(self, p, np.asarray(d[p]))
            else:
                setattr(self, p, np.asarray([]))
        if hasattr(self, "frang") and hasattr(self, "slist") and hasattr(self, "rsep"):
            setattr(
                self, 
                "srange", 
                np.asarray(self.frang + self.slist*self.rsep)
            )
        if hasattr(self, "intt.sc") and hasattr(self, "intt.us"):
            setattr(
                self,
                "intt",
                getattr(self, "intt.sc") + 1.0e-6 * getattr(self, "intt.us")
            )
        if self.nrange_scatter is not None:
            sl, gf = np.asarray(self.slist), np.asarray(self.gflg)
            gf[sl <= self.nrange_scatter] = -1
            self.gflg = np.asarray(gf)
        self.time = time
        return

    def copy(self, bm):
        """Copy all parameters"""
        for p in bm.__dict__.keys():
            setattr(self, p, getattr(bm, p))
        return


class Scan(object):
    """Class to hold one scan (multiple beams)"""

    def __init__(self, stime=None, etime=None):
        """
        initialize the parameters which will be stored
        stime: start time of scan
        etime: end time of scan
        """
        self.stime = stime
        self.etime = etime
        self.beams = []
        return

    def update_time(self):
        """
        Update stime and etime of the scan.
        up: Update average parameters if True
        """
        self.stime = min([b.time for b in self.beams])
        self.etime = max([b.time for b in self.beams])
        self.scan_time = 60*np.rint((self.etime-self.stime).total_seconds()/60.)
        for b in self.beams:
            setattr(b, "scan_time", self.scan_time)
        return