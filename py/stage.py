#!/usr/bin/env python

"""
    stage.py: module to stagging data into raw format
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import sys

import datetime as dt
import json
import multiprocessing as mp
import os
import time
import traceback
from functools import partial
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pydarn
import utils
from boxcar import BxFilter
from fetch_fit_data import FetchData
from loguru import logger
from plots import RTI, IntervalPlots, histogram_plots
from scipy.fft import rfft, rfftfreq
from scipy.interpolate import interp1d
from scipy.signal import get_window
from scipy.stats import t as T
from tqdm import tqdm


class StagingUnit(object):
    """
    Store all the Raw SuperDARN data to local or
    ARC (remote super computer) for further processing.
    """

    def __init__(self, rad, dates, conf):
        """
        Parameters
        ----------
        rad - 3 char radar code
        dates - [sdate, edate]
        """
        tqdm.pandas()
        self.proc_start_time = time.time()
        self.rad = rad
        self.dates = dates
        self.conf = conf
        self.modify_config(dates)
        self.data_exists = False
        self.types = ["raw", "filt", "dtrnd"]
        hdw_data = pydarn.read_hdw_file(self.rad)
        self.lats, self.lons = pydarn.Coords.GEOGRAPHIC(hdw_data.stid)
        self._fetch()
        return

    def modify_config(self, dates):
        """
        Modify Config Parameters
        """
        self.base = self.conf.files.base.format(
            date=self.dates[0].strftime("%Y-%m-%d"), run_id=self.conf.run_id
        )
        if not os.path.exists(self.base):
            os.makedirs(self.base, exist_ok=True)
        return

    def fetch_file(self, p):
        """
        Return filename
        """
        stime, etime = self.dates[0].strftime("%H%M"), self.dates[-1].strftime("%H%M")
        file = self.base + self.conf.files.csv.format(
            rad=self.rad, stime=stime, etime=etime, kind=p
        )
        return file

    def __get_location__(self, row):
        """
        Get AACGMV2 magnetic location
        """
        lat, lon, dtime = (
            self.lats[row["slist"], row["bmnum"]],
            self.lons[row["slist"], row["bmnum"]],
            row["time"],
        )
        row["glat"], row["glon"] = lat, lon
        return row

    def _create_sequence_(self, scans):
        self.fscans = []
        fscans = [scans[i : i + 3] for i in range(len(scans) - 3)]
        p0 = mp.Pool(self.conf.cores)
        self.filter, j = BxFilter(thresh=self.conf.filter.thresh), 0
        partial_filter = partial(self.filter.doFilter)
        for s in p0.map(partial_filter, fscans):
            if np.mod(j, 10) == 0:
                logger.info(f"Running median filter for scan at - {j}")
            self.fscans.append(s)
            j += 1
        return

    def _fetch(self):
        """
        Fetch radar data from repository
        """
        logger.info(" Fetching data...")
        dates = (
            [
                self.dates[0] - dt.timedelta(minutes=self.conf.filter.w_mins / 2),
                self.dates[1] + dt.timedelta(minutes=self.conf.filter.w_mins / 2),
            ]
            if self.conf.filter.w_mins is not None
            else self.dates
        )
        logger.info(
            f" Read radar file {self.rad} for {[d.strftime('%Y.%m.%dT%H.%M') for d in self.dates]}"
        )
        self.fd = FetchData(self.rad, dates, ftype=self.conf.ftype)
        if not os.path.exists(self.fetch_file("raw")):

            _, scans, self.data_exists = self.fd.fetch_data(by="scan")
            if self.data_exists:
                # self._create_sequence_(scans)
                scan_time = (scans[0].etime - scans[0].stime).total_seconds()
                self.frame = self.fd.scans_to_pandas(scans)
                if len(self.frame) > 0:
                    self.frame["srange"] = self.frame.frang + (
                        self.frame.slist * self.frame.rsep
                    )
                    self.frame["intt"] = (
                        self.frame["intt.sc"] + 1.0e-6 * self.frame["intt.us"]
                    )
                    self.frame = self.frame.progress_apply(
                        self.__get_location__, axis=1
                    )
                    self.frame["scan_time"] = scan_time
                    self._save(self.frame, "raw")
                else:
                    logger.info(f" Radar file does not exist!")
            else:
                self.data_exists = False
        else:
            self.frame = pd.read_csv(self.fetch_file("raw"), parse_dates=["time"])
            self.data_exists = True
        return

    def _save(self, df, p):
        """
        Print data into parquet/gzip file for later processing.
        """
        df.to_csv(self.fetch_file(p), index=False, header=True)
        return

    def fetch_file(self, p):
        """
        Return filename
        """
        stime, etime = self.dates[0].strftime("%H%M"), self.dates[-1].strftime("%H%M")
        file = self.base + self.conf.files.csv.format(
            rad=self.rad, stime=stime, etime=etime, kind=p
        )
        return file


class Filter(StagingUnit):
    """
    For each 1-hour interval, select only ionospheric backscatter
    and reject ground scatter by requiring backscatter to satisfy
    one of the following conditions:
        a) Based on Doppler velocity|VLOS| <= 30 m/s or
           Backscatter is flagged as ionospheric scatter by traditional method.
        b) Power is greater than 3 dB and Errors in VLOS & W < 100 m/s
        c) Remove near-range data (slant range <= 765 km)
    """

    def __init__(self, rad, dates, conf, filters=["a", "b", "c"]):
        """
        Parameters
        ----------
        rad - 3 char radar code
        dates - [sdate, edate]
        conf - Configuration obj
        filters - combinations of three filtering criteria
        """
        self.proc_start_time = time.time()
        super().__init__(rad, dates, conf)
        self.filters = filters
        if self.data_exists:
            self._detrnd()
            # self._analysis()
        self._compile_sample_output()
        return

    def _detrnd(self):
        """
        Detreand the dataset with w_min window
        """
        if not os.path.exists(self.fetch_file("dtrnd")):
            logger.info(f" Start detranding data {len(self.frame)}")
            self.dtnd_frame = self.frame.progress_apply(self.__trnd_support__, axis=1)
            self._save(self.dtnd_frame, "dtrnd")
        else:
            self.dtnd_frame = pd.read_csv(
                self.fetch_file("dtrnd"), parse_dates=["time"]
            )
        logger.info(f" Done Dtrnd {len(self.dtnd_frame)}")
        return

    def __trnd_support__(self, row):
        """
        Detrend the velocity data - Support Function.
        """
        w_mins = self.conf.filter.w_mins
        pval = np.nan
        b, g, t, p_raw = row["bmnum"], row["slist"], row["time"], row["p_l"]
        if (t >= self.dates[0]) & (t <= self.dates[-1]):
            t_start, t_end = t - dt.timedelta(minutes=w_mins / 2), t + dt.timedelta(
                minutes=w_mins / 2
            )
            x = self.frame[
                (self.frame.bmnum == b)
                & (self.frame.slist == g)
                & (self.frame.time >= t_start)
                & (self.frame.time >= t_end)
            ]
            pval = p_raw - np.nanmedian(x["p_l"])
        row["p_l"] = pval
        return row

    def _analysis(self):
        """
        Analyze the data
        1. Run resample and FFT algo
        2. Run MUSIC
        """
        self.dtnd_frame.dropna(inplace=True)
        stime, etime = self.dates[0].strftime("%H%M"), self.dates[-1].strftime("%H%M")
        file = self.base + self.conf.files.pkl.format(
            rad=self.rad, stime=stime, etime=etime
        )
        self.pg = Periodograms(self.rad, self.dates, self.conf, self.dtnd_frame, file)
        return

    def _compile_sample_output(self, beam=13, fq=[(10, 35)]):
        """
        Create figures:
        1. RTI plots.
        2. Series plot.
        """
        rti = RTI(100, self.dates, fig_title="", num_subplots=3)
        stime, etime = self.dates[0].strftime("%H%M"), self.dates[-1].strftime("%H%M")
        rti_file = self.base + self.conf.files.png.format(
            rad=self.rad, stime=stime, etime=etime, kind="rti"
        )
        rti.addParamPlot(
            self.frame,
            beam,
            self.rad.upper()
            + ":"
            + str(beam)
            + ", "
            + self.dates[0].strftime("%Y-%m-%d"),
        )
        #         rti.addParamPlot(
        #             self.dtnd_frame,
        #             beam,
        #             self.rad.upper()
        #             + ":"
        #             + str(beam)
        #             + ", "
        #             + self.dates[0].strftime("%Y-%m-%d"),
        #             p_max=21, p_min=-20
        #         )
        # rti.add_series(self.dtnd_frame, fq[0][0], fq[0][1])
        rti.save(rti_file)
        rti.close()
        #         ip = IntervalPlots(self.pg.results, fq)
        #         series_file = self.base + self.conf.files.png.format(
        #             rad=self.rad, stime=stime, etime=etime, kind="series"
        #         )
        #         ip.save(series_file)
        #         ip.close()
        #         file = self.base + "histograms.png"
        #         histogram_plots(self.pg.frequencies, file)
        return


class StagingHopper(object):
    """
    For each RBSP log entry in the log files
    this class forks a process to fetch fitacf
    data from the SD repository and store to
    local or ARC.
    """

    def __init__(self, _filestr="config/record.logs/*.log", cores=4, run_first=None):
        """
        Params
        ------
        _filestr - Regular expression to search files
        cores - Mutiprocessing cores
        run_first - Firt N modes to run
        """
        tqdm.pandas()
        self._filestr = _filestr
        self.cores = cores
        self.run_first = run_first
        self._logs = utils.read_logs(_filestr)
        self.load_params()
        self._run()
        return

    def load_params(self):
        """
        Load parameters from params.json
        """
        with open("config/params.json") as f:
            self.conf = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
        return

    def _proc(self, o):
        """
        Process method to invoke stagging unit
        """
        try:
            rad, dates = o["rad"], [o["ut_start"], o["ut_end"]]
            logger.info(
                f"Filtering radar {rad} for {[d.strftime('%Y.%m.%dT%H.%M') for d in dates]}"
            )
            Filter(rad, dates, self.conf)
        except:
            err = traceback.format_exc()
            logger.error(
                f"Error in filtering radar {rad} for {[d.strftime('%Y.%m.%dT%H.%M') for d in dates]} \n {err}"
            )
        return

    def _run(self):
        """
        Run parallel threads to access and save data
        """
        if self.run_first:
            logger.info(f"Start parallel procs for first {self.run_first} entries")
        else:
            logger.info(f"Start parallel procs for first {len(self._logs)} entries")
        rlist = self._logs if self.run_first is None else self._logs[: self.run_first]
        for f in rlist.to_dict(orient="records"):
            self._proc(f)
        logger.info(f"Job complete!")
        return


if __name__ == "__main__":
    "__main__ function"
    start = time.time()
    StagingHopper(run_first=None)
    end = time.time()
    logger.info(f" Interval time {'%.2f'%(end - start)} sec.")
