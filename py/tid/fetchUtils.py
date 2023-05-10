#!/usr/bin/env python

"""fetchUtils.py: utility module to fetch fitacf<v> level data."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import bz2
import datetime as dt
import glob
import os

import numpy as np
import pandas as pd
import pydarn
import tidUtils
from dataUtils import Beam, Gate, Scan
from fanUtils import Fan
from filterUtils import Boxcar
from loguru import logger
from rtiUtils import RTI, plot_SDTEC_TS
from tqdm import tqdm


class FetchData(object):
    """Class to fetch data from fitacf files for one radar for atleast a day"""

    def __init__(
        self,
        rad,
        date_range,
        ftype="fitacf",
        files=None,
        verbose=True,
        nrange_scatter=None,
    ):
        """
        initialize the vars
        rad = radar code
        date_range = [ start_date, end_date ]
        files = List of files to load the data from
        e.x :   rad = "sas"
                date_range = [
                    datetime.datetime(2017,3,17),
                    datetime.datetime(2017,3,18),
                ]
        """
        self.rad = rad
        self.date_range = date_range
        self.files = files
        self.verbose = verbose
        self.regex = tidUtils.get_config("SD_DATA")
        self.ftype = ftype
        if (rad is not None) and (date_range is not None) and (len(date_range) == 2):
            self._create_files()
        self.s_params = [
            "bmnum",
            "noise.sky",
            "tfreq",
            "scan",
            "nrang",
            "intt.sc",
            "intt.us",
            "mppul",
            "rsep",
            "cp",
            "frang",
            "smsep",
            "lagfr",
            "channel",
            "mplgs",
            "nave",
            "noise.search",
            "mplgexs",
            "xcf",
            "noise.mean",
            "ifmode",
            "bmazm",
            "rxrise",
            "mpinc",
        ]
        self.nrange_scatter = nrange_scatter
        self.v_params = ["v", "w_l", "gflg", "p_l", "slist"]
        self.hdw_data = pydarn.read_hdw_file(self.rad)
        self.lats, self.lons = pydarn.Coords.GEOGRAPHIC(self.hdw_data.stid)
        return

    def _create_files(self):
        """
        Create file names from date and radar code
        """
        if self.files is None:
            self.files = []
        reg_ex = self.regex
        days = (self.date_range[1] - self.date_range[0]).days + 2
        ent = -1
        for d in range(-1, days):
            e = self.date_range[0] + dt.timedelta(days=d)
            fnames = sorted(
                glob.glob(
                    reg_ex.format(
                        year=e.year,
                        rad=self.rad,
                        ftype=self.ftype,
                        date=e.strftime("%Y%m%d"),
                    )
                )
            )
            for fname in fnames:
                tm = fname.split(".")[1]
                sc = fname.split(".")[2]
                dus = dt.datetime.strptime(
                    fname.split(".")[0].split("/")[-1] + tm + sc, "%Y%m%d%H%M%S"
                )
                due = dus + dt.timedelta(hours=2)
                if (ent == -1) and (dus <= self.date_range[0] <= due):
                    ent = 0
                if ent == 0:
                    self.files.append(fname)
                if (ent == 0) and (dus <= self.date_range[1] <= due):
                    ent = -1
        return

    def _parse_data(self, data, s_params, v_params, by):
        """
        Parse data by data type
        data: list of data dict
        params: parameter list to fetch
        by: sort data by beam or scan
        """
        _b, _s = [], []
        if self.verbose:
            logger.info("Started converting to beam data.")
        for d in data:
            time = dt.datetime(
                d["time.yr"],
                d["time.mo"],
                d["time.dy"],
                d["time.hr"],
                d["time.mt"],
                d["time.sc"],
                d["time.us"],
            )
            if time >= self.date_range[0] and time <= self.date_range[1]:
                bm = Beam(self.nrange_scatter)
                bm.set(time, d, s_params, v_params)
                _b.append(bm)
        if self.verbose:
            logger.info("Converted to beam data.")
        if by == "scan":
            if self.verbose:
                logger.info("Started converting to scan data.")
            scan, sc = 0, Scan(None, None)
            sc.beams.append(_b[0])
            for _ix, d in enumerate(_b[1:]):
                if d.scan == 1 and d.time != _b[_ix].time:
                    sc.update_time()
                    _s.append(sc)
                    sc = Scan(None, None)
                    sc.beams.append(d)
                else:
                    sc.beams.append(d)
            sc.update_time()
            _s.append(sc)
            if self.verbose:
                logger.info("Converted to scan data.")
        return _b, _s, True

    def convert_to_pandas(
        self,
        beams,
    ):
        """
        Convert the beam data into dataframe
        """
        if "time" not in self.s_params:
            self.s_params.append("time")
        _o = dict(
            zip(
                self.s_params + self.v_params,
                ([] for _ in self.s_params + self.v_params),
            )
        )
        for b in beams:
            l = len(getattr(b, "slist"))
            for p in self.v_params:
                _o[p].extend(getattr(b, p))
            for p in self.s_params:
                _o[p].extend([getattr(b, p)] * l)
        L = len(_o["slist"])
        for p in self.s_params + self.v_params:
            if len(_o[p]) < L:
                l = len(_o[p])
                _o[p].extend([np.nan] * (L - l))
        return pd.DataFrame.from_records(_o)

    def scans_to_pandas(
        self,
        scans,
        start_scnum=0,
        flt=False,
    ):
        """
        Convert the scan data into dataframe
        """
        if "time" not in self.s_params:
            self.s_params.append("time")
        if "srange" not in self.v_params:
            self.v_params.append("srange")
        if "intt" not in self.s_params:
            self.s_params.append("intt")
        _o = dict(
            zip(
                self.s_params + self.v_params + ["scnum", "scan_time"],
                ([] for _ in self.s_params + self.v_params + ["scnum", "scan_time"]),
            )
        )
        for idn, s in enumerate(scans):
            for b in s.beams:
                l = len(getattr(b, "slist"))
                for p in self.v_params:
                    _o[p].extend(getattr(b, p))
                for p in self.s_params:
                    _o[p].extend([getattr(b, p)] * l)
                _o["scnum"].extend([idn + start_scnum] * l)
                _o["scan_time"].extend([getattr(b, "scan_time")] * l)
            L = len(_o["slist"])
            for p in self.s_params + self.v_params:
                if len(_o[p]) < L:
                    l = len(_o[p])
                    _o[p].extend([np.nan] * (L - l))
        return pd.DataFrame.from_records(_o)

    def __get_location__(self, row):
        """
        Get locations
        """
        lat, lon, dtime = (
            self.lats[row["slist"], row["bmnum"]],
            self.lons[row["slist"], row["bmnum"]],
            row["time"],
        )
        row["glat"], row["glon"] = lat, lon
        return row

    def pandas_to_beams(
        self,
        df,
    ):
        """
        Convert the dataframe to beam
        """
        if "time" not in self.s_params:
            self.s_params.append("time")
        beams = []
        for bm in np.unique(df.bmnum):
            o = df[df.bmnum == bm]
            d = o.to_dict(orient="list")
            for p in self.s_params:
                d[p] = d[p][0]
            b = Beam(self.nrange_scatter)
            b.set(o.time.tolist()[0], d, self.s_params, self.v_params)
            beams.append(b)
        return beams

    def pandas_to_scans(
        self,
        df,
    ):
        """
        Convert the dataframe to scans
        """
        if "time" not in self.s_params:
            self.s_params.append("time")
        scans = []
        for sn in np.unique(df.scnum):
            o = df[df.scnum == sn]
            beams = self.pandas_to_beams(o)
            sc = Scan(None, None)
            sc.beams.extend(beams)
            sc.update_time()
            scans.append(sc)
        return scans

    def fetch_data(
        self,
        by="beam",
    ):
        """
        Fetch data from file list and return the dataset
        params: parameter list to fetch
        by: sort data by beam or scan
        """
        data = []
        for f in self.files:
            with bz2.open(f) as fp:
                fs = fp.read()
            if self.verbose:
                logger.info(f"File:{f}")
            reader = pydarn.SuperDARNRead(fs, True)
            records = reader.read_fitacf()
            data += records
        if (by is not None) and (len(data) > 0):
            data = self._parse_data(data, self.s_params, self.v_params, by)
            return data
        else:
            return (None, None, False)

    def plot_RTI(
        self,
        beams=[],
        nGates=100,
        date_range=None,
        angle_th=100.0,
        vhm=None,
        tec_mat_file=None,
        tec_param="cdvTECgrid2",
        power_vlim=[20, 50],
        tec_vlim=[-0.05, 0.05],
    ):
        """
        Plot RTI plots by beams
        """
        logger.info("Inside RTI plots")
        date_range = date_range if date_range else self.date_range
        beams = beams if beams and len(beams) > 0 else self.frame.bmnum.unique()

        # Reload dataset
        file = os.path.join(tidUtils.get_folder(self.date_range[0]), f"{self.rad}.csv")
        mfile = os.path.join(
            tidUtils.get_folder(self.date_range[0]), f"{self.rad}_med.csv"
        )
        self.frame = pd.read_csv(file, parse_dates=["time"])
        self.medframe = pd.read_csv(mfile, parse_dates=["time"])

        if tec_mat_file:
            tec, tec_times = tidUtils.read_tec_mat_files(tec_mat_file)

        for b in beams:
            logger.info(f" Beam:{b}")
            file = (
                tidUtils.get_folder(self.date_range[0]) + f"/{self.rad}-{'%02d'%b}.png"
            )
            rt = RTI(
                100,
                date_range,
                (self.lats, self.lons),
                [date_range[0] + dt.timedelta(hours=14), date_range[1]],
                f"{self.date_range[0].strftime('%Y-%m-%d')}/{self.rad}/{b}",
                num_subplots=3,
                angle_th=angle_th,
                vhm=vhm,
            )
            rt.addParamPlot(
                self.frame,
                b,
                "",
                xlabel="",
                cbar=True,
                plot_fov=False,
                vlim=power_vlim,
            )
            if tec_mat_file:
                rt.ovearlay_TEC(
                    tec,
                    tec_times,
                    b,
                    plot_fov=False,
                    xlabel="",
                    tec_param=tec_param,
                    vlim=tec_vlim,
                    cmap="magma",
                )

                axis = rt.ovearlay_TEC(
                    tec,
                    tec_times,
                    b,
                    cbar_xOff=0.15,
                    tec_param=tec_param,
                    plot_fov=False,
                    vlim=tec_vlim,
                    cmap="magma",
                )
                rt.addParamPlot(
                    self.frame,
                    b,
                    "",
                    xlabel="",
                    plot_fov=False,
                    ax=axis,
                    vlim=power_vlim,
                )
            rt.save(file)
            rt.close()
        return

    def TS1D(
        self,
        bm,
        rg,
        date_range=None,
        tec_mat_file=None,
        tec_param="cdvTECgrid2",
        power_vlim=[20, 50],
        tec_vlim=[-0.05, 0.05],
    ):
        """
        Conduct 1D TS analysis
        """
        logger.info("Inside TS1D plots")
        date_range = date_range if date_range else self.date_range

        # Reload dataset
        file = os.path.join(tidUtils.get_folder(self.date_range[0]), f"{self.rad}.csv")
        mfile = os.path.join(
            tidUtils.get_folder(self.date_range[0]), f"{self.rad}_med.csv"
        )
        self.frame = pd.read_csv(file, parse_dates=["time"])
        self.medframe = pd.read_csv(mfile, parse_dates=["time"])

        if tec_mat_file:
            tec, tec_times = tidUtils.read_tec_mat_files(tec_mat_file)

        logger.info(f" Beam:{bm}, Gate:{rg}")
        srange = rg * self.frame.rsep.tolist()[0] + self.frame.frang.tolist()[0]
        sd_frame = self.frame[(self.frame.bmnum == bm) & (self.frame.srange == srange)]
        tec_frame = tidUtils.fetch_tec_by_beam(tec, bm, param=tec_param)[:, rg]
        fname = os.path.join(
            tidUtils.get_folder(self.date_range[0]), f"{self.rad}-{bm}-{rg}.png"
        )
        plot_SDTEC_TS(
            sd_frame.time,
            # tidUtils.smooth(sd_frame["p_l"], replace_nans=True, byvalue=None),
            sd_frame["p_l"],
            tec_times[: len(tec_frame)],
            # tidUtils.smooth(tec_frame),
            tec_frame,
            fname,
            [date_range[0] + dt.timedelta(hours=14), date_range[1]],
            txt=f"{date_range[0].strftime('%Y-%m-%d')} / {self.rad} / Beam:{bm} / Range:{srange}",
            sd_ylim=power_vlim,
            tec_ylim=tec_vlim,
        )
        return

    def plot_FoV(self, scan_num=None, date=None, tec_mat_file=None):
        """
        Plot FoV plots by scan_num/date
        """
        if scan_num:
            scan = self.scans[scan_num]
            file = (
                tidUtils.get_folder(scan.stime)
                + f"/{self.rad}-Fan-{scan.stime.strftime('%H.%M')}.png"
            )
            fov = Fan([self.rad], scan.stime)
            fov.generate_fov(self.rad, self.frame, tec_mat_file=tec_mat_file)
            fov.save(file)
            fov.close()
        elif date:
            file = (
                tidUtils.get_folder(scan.stime)
                + f"/{self.rad}-Fan-{date.strftime('%H.%M')}.png"
            )
            fov = Fan([self.rad], date)
            fov.generate_fov(self.rad, self.frame, tec_mat_file=tec_mat_file)
            fov.save(file)
            fov.close()
        else:
            tec, tec_times = tidUtils.read_tec_mat_files(tec_mat_file)
            for scan in self.scans:
                file = (
                    tidUtils.get_folder(scan.stime)
                    + f"/{self.rad},{scan.stime.strftime('%H-%M')}.png"
                )
                fov = Fan([self.rad], scan.stime, tec=tec, tec_times=tec_times)
                fov.generate_fov(self.rad, self.frame, laytec=True)
                fov.save(file)
                fov.close()
        return

    def scanGenerator(self, df=None, scans=None):
        """
        Generate yield generators for scans
        """
        if df is not None:
            scans = self.pandas_to_scans(df)
        for sc in scans:
            yield sc

    def to_netcdf(self, scans, fname, params=["p_l", "gflg"], near_range_scatter=7):
        """
        Convert to NetCDF files
        """
        print(self.lats.shape)
        def extract_3D_data(px="p_l"):
            dat = (
                np.zeros((len(scans), self.lats.shape[1], self.lats.shape[0])) * np.nan
            )
            for i, s in enumerate(scans):
                for j, b in enumerate(s.beams):
                    if len(getattr(b, "slist")) > 0:
                        dat[i, j, getattr(b, "slist")] = getattr(b, px)
            dat = np.ma.masked_invalid(dat)
            return dat

        from netCDF4 import Dataset

        ds = Dataset(fname, "w")
        # Create all dataset dinemtions
        beams = ds.createDimension("beams", self.lats.shape[1])
        gates = ds.createDimension("gates", self.lats.shape[0])
        lscans = ds.createDimension("scans", len(scans))
        # Create all dataset
        beams = ds.createVariable("beams", "i2", ("beams",))
        gates = ds.createVariable("gates", "i2", ("gates",))
        lscans = ds.createVariable("scans", "i2", ("scans",))
        beams[:] = np.arange(self.lats.shape[1])
        gates[:] = np.arange(self.lats.shape[0])
        lscans[:] = np.arange(len(scans))
        lat_grid = ds.createVariable("gdlat", "f4", ("beams", "gates"))
        lon_grid = ds.createVariable("glong", "f4", ("beams", "gates"))
        lat_grid[:], lon_grid[:] = self.lats.T, self.lons.T
        lat_grid.units, lon_grid.units = "deg", "deg"
        stimes = ds.createVariable("start_time", "i4", ("scans",))
        etimes = ds.createVariable("end_time", "i4", ("scans",))
        stimes[:], etimes[:] = (
            [(s.stime - self.date_range[0]).total_seconds() for s in scans],
            [(s.etime - self.date_range[0]).total_seconds() for s in scans],
        )
        stimes.units, etimes.units = (
            f"Sec since {self.date_range[0].strftime('%Y-%m-%dT%H:%M:%S')}",
            f"Sec since {self.date_range[0].strftime('%Y-%m-%dT%H:%M:%S')}",
        )
        for param in params:
            p = ds.createVariable(param, "f4", ("scans", "beams", "gates"))
            p[:] = extract_3D_data(param)
        ds.close()
        return

    def med_filter(self, param, fname):
        """
        Run median filtering
        """
        bx = Boxcar(
            thresh=param["thresh"],
            w=None,
            nrange_scatter=self.nrange_scatter,
        )
        scan_stacks = [self.scans[i - 1 : i + 2] for i in range(1, len(self.scans) - 1)]
        self.filtered_scans = [bx.do_filter(scan_stack) for scan_stack in scan_stacks]
        self.medframe = self.scans_to_pandas(self.filtered_scans)
        self.medframe = self.medframe.progress_apply(self.__get_location__, axis=1)
        logger.info(f"m-Data length {self.rad}: {len(self.medframe)}")
        self.medframe.to_csv(fname, index=False, header=True, float_format="%g")
        return

    def to_geom(self):
        """
        Save lat/lons. csv
        """
        file = os.path.join(
            tidUtils.get_folder(self.date_range[0]), f"{self.rad}_geom.csv"
        )
        o = pd.DataFrame()
        rsep, frang = self.frame.rsep.iloc[0], self.frame.frang.iloc[0]
        for b in range(self.lats.shape[1]):
            o["beam%02d.lat" % (b + 1)] = self.lats[:, b]
            o["beam%02d.lon" % (b + 1)] = self.lons[:, b]
            o["beam%02d.slist" % (b + 1)] = range(len(self.lats[:, b]))
            o["beam%02d.srange" % (b + 1)] = np.arange(len(self.lats[:, b]))*rsep + frang
        o.to_csv(file, index=False, header=True, float_format="%g")
        return

    @staticmethod
    def fetch(
        rad,
        date_range,
        ftype="fitacf",
        files=None,
        verbose=True,
        med_filter=None,
        nrange_scatter=7,
        to_scan=False,
        to_med_fetch=False,
    ):
        """
        Static method to fetch datasets
        """
        tqdm.pandas()
        file = os.path.join(tidUtils.get_folder(date_range[0]), f"{rad}.csv")
        mfile = os.path.join(tidUtils.get_folder(date_range[0]), f"{rad}_med.csv")
        nfile = os.path.join(tidUtils.get_folder(date_range[0]), f"{rad}.nc")
        fd = FetchData(rad, date_range, ftype, files, verbose, nrange_scatter)
        
        if os.path.exists(file):
            fd.frame = pd.read_csv(file, parse_dates=["time"])
            logger.info(f"Data length {rad}: {len(fd.frame)}")
            if to_scan:
                fd.scans = fd.pandas_to_scans(fd.frame)
                logger.info(f"# Scans {rad}: {len(fd.scans)}")
            if to_med_fetch and os.path.exists(mfile):
                fd.medframe = pd.read_csv(mfile, parse_dates=["time"])
                logger.info(f"m-Data length {rad}: {len(fd.medframe)}")
                fd.filtered_scans = fd.pandas_to_scans(fd.medframe)
                logger.info(f"# Scans {rad}: {len(fd.filtered_scans)}")
        else:
            #if med_filter:
            _, scans, data_exists = fd.fetch_data(by="scan")
            if data_exists:
                fd.frame = fd.scans_to_pandas(scans)
                fd.scans = scans
                logger.info(f"Data length {rad}: {len(fd.frame)}")
                if len(fd.frame) > 0:
                    fd.frame = fd.frame[
                        (fd.frame.slist<101)
                        & (fd.frame.bmnum<np.max(fd.frame.bmnum))
                    ]
                    fd.frame = fd.frame.progress_apply(fd.__get_location__, axis=1)
                    fd.frame.to_csv(file, index=False, header=True, float_format="%g")
                    # fd.to_netcdf(fd.scans, nfile)
                    # fd.med_filter(med_filter, mfile)
                    # fd.to_netcdf(fd.filtered_scans, nfile.replace(".nc", "_med.nc"))
                else:
                    logger.info(f"Data does not exists: {rad}!")
        return fd
