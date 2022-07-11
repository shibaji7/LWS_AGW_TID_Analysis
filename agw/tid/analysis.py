#!/usr/bin/env python

"""
    analysis.py: module to analyze the data
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import os
import pickle
from loguru import logger

import numpy as np
from astropy.timeseries import LombScargle
from tqdm import tqdm


class Periodograms(object):
    """
    This class takes filtered or detranded dataset
    and estimate the wave frequencies and wavelengths.
    """

    def __init__(self, rad, dates, conf, df, file):
        """
        First estimate avalibale frequencies in each
        2-hour window. Then calculate possible wavelengths.

        Parameters
        ----------
        rad - 3 char radar code
        dates - [sdate, edate]
        conf - Configuration obj
        df - Dataframe
        file - Name of the output file
        """
        self.rad = rad
        self.dates = dates
        self.conf = conf
        self.df = df
        self.pname = self.conf.periodograms.pname
        self.results = {}
        self.file = file
        if not os.path.exists(file):
            self._compute_freq_()
            self._compute_wv_()
            self._save_()
        else:
            self._load_()
        return

    def _compute_freq_(self):
        """
        Compute frequencies for range cell.
        """
        conf = self.conf.periodograms.freq
        slist, beams = self.df.slist.unique(), self.df.bmnum.unique()
        result = {}
        for bm in tqdm(beams, desc="On beams"):
            if bm not in result.keys():
                result[bm] = {}
            for gate in tqdm(slist, desc="On gate"):                
                u = self.df[(self.df.bmnum == bm) & (self.df.slist == gate)]
                if self.__check_eligibility__(u, conf, ty="f"):
                    logger.info(f"Beam, Gate : {bm},{gate}")
                    ts = u["time"].tolist()[0]
                    t, y = np.array(
                        [(x - ts).total_seconds() for x in u["time"]]
                    ), np.array(u[self.pname])
                    LS = LombScargle(
                        t,
                        y,
                        fit_mean=conf.ls.fit_mean,
                        center_data=conf.ls.center_data,
                        nterms=conf.ls.nterms,
                        normalization=conf.ls.normalization,
                    )
                    f = np.linspace(1e-6, 1e-2, 10000)
                    power = LS.power(f)
                    result[bm][gate] = {"f": f, "pow": power}
        self.results["freq"] = result
        return

    def _compute_wv_(self):
        """
        Compute wavelengths for time instances cell.
        """
        conf = self.conf.periodograms.wv
        beams = self.df.bmnum.unique()
        itrs = int(
            (self.dates[-1] - self.dates[0]).total_seconds() / (60 * conf.dtau_mins)
        )

        result = {}
        for bm in tqdm(beams, desc="On beams"):
            if bm not in result.keys():
                result[bm] = {}
            for i in tqdm(range(itrs - 1), desc="On time"):
                logger.info(f"Beam, Time : {bm},{i}")
                t_start = self.dates[0] + dt.timedelta(minutes=i * conf.dtau_mins)
                t_end = self.dates[0] + dt.timedelta(minutes=(i + 1) * conf.dtau_mins)
                u = self.df[
                    (self.df.bmnum == bm)
                    & (self.df.time >= t_start)
                    & (self.df.time <= t_end)
                ]
                if self.__check_eligibility__(u, conf, ty="w"):
                    l, y = np.array(u["srange"])*1e3, np.array(u[self.pname])
                    LS = LombScargle(
                        l,
                        y,
                        fit_mean=conf.ls.fit_mean,
                        center_data=conf.ls.center_data,
                        nterms=conf.ls.nterms,
                        normalization=conf.ls.normalization,
                    )
                    w = np.linspace(1e-10, 1e-4, 10000)
                    power = LS.power(w)
                    result[bm][i] = {"wvn": w, "pow": power}
        self.results["wv"] = result
        return

    def __check_eligibility__(self, u, conf, ty="f"):
        """
        Check the DF is eligible for LS Periodograms
        """
        e = False
        if ty == "f":
            min_no_echoes = conf.min_no_echoes
            min_window_hour = conf.min_window_hour
            if len(u) >= min_no_echoes:
                dtau = (
                    np.abs((u.time.tolist()[-1] - u.time.tolist()[0]).total_seconds())
                    / 3600
                )
                if dtau >= min_window_hour:
                    e = True
        if ty == "w":
            min_no_echoes = conf.min_no_echoes
            if len(u) >= min_no_echoes:
                e = True
        return e

    def _save_(self):
        """
        Save all dataset to local
        """
        with open(self.file, "wb") as h:
            pickle.dump(self.results, h, protocol=pickle.HIGHEST_PROTOCOL)
        return
    
    def _load_(self):
        """
        Load all dataset from local
        """
        with open(self.file, "rb") as h:
            self.results = pickle.load(h)
        return
