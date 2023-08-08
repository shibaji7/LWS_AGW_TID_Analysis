#!/usr/bin/env python3

"""gitm.py: fetch GEMINI3D data"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"


from loguru import logger
import datetime as dt
import glob
import h5py
import numpy as np
import utils
import pandas as pd
from scipy.io import savemat


class GEMINI2D(object):
    
    def __init__(self, date_str, cfg, folder="dataset/GEMINI3D/", param="nsall", grid_file="grid.mat"):
        self.date_str = date_str
        self.cfg = cfg
        self.folder = folder
        self.param = param
        self.grid_file = grid_file
        self.load_grid()
        self.load_dataset()
        return
    
    def load_grid(self):
        """
        Load grid file
        """
        self.grid = pd.DataFrame()
        with h5py.File(self.folder + self.grid_file, "r") as fkey:
            self.grid["glat"] = np.array(fkey.get("glat")[:]).tolist()
            self.grid["glon"] = np.mod( (fkey.get("glon")[:] + 180), 360 ) - 180
            self.grid["alt"] = fkey.get("alt")[:]
            self.grid["glat"] = fkey.get("glat")[:]
        return
    
    def load_dataset(self):
        """
        Load all dataset available
        """
        self.files = glob.glob(self.folder + self.date_str + "*")
        self.files.sort()
        self.dataset = dict()
        for fname in self.files:
            with h5py.File(fname, "r") as fkey:
                day = dt.datetime.strptime(fname.split("_")[0].split("/")[-1], "%Y%m%d")
                seconds = int(fname.split("_")[1].split(".")[0])
                day = day + dt.timedelta(seconds=seconds)
                o = self.grid.copy()
                o[self.param] = fkey.get(self.param)[:]
                self.dataset[day] = o
                logger.info(f"Load dataset: {day}")
        return
    
    def fetch_dataset_by_locations(self, time, lats, lons, alts, dh=500, dlat=0.2, dlon=0.2):
        """
        Fetch data by lat/lon limits
        """
        logger.info(f"Time: {time}")
        df = self.dataset[time]
        df = df[df.alt>=0]
        out, ix = np.zeros((len(alts), len(lons))) * np.nan, 0
        for lat, lon in zip(lats, lons):
            uf = df[
                (df.glat>=lat-dlat)
                & (df.glat<=lat+dlat)
                & (df.glon>=lon-dlon)
                & (df.glon<=lon+dlon)
            ]
            uf = uf.groupby("alt").mean().reset_index()
            if len(uf) > 0:
                out[:, ix] = utils.interpolate_by_altitude(
                        np.array(uf.alt)/1e3, alts, utils.smooth(np.array(uf[self.param]), 51),
                        self.cfg.scale, self.cfg.kind, method="extp"
                    ) * 1e-6
            ix += 1
        self.param_val = out
        return out
    
    @staticmethod
    def create_object(
        date_str,
        cfg,
        folder="dataset/GEMINI3D/",
        param="nsall",
        grid_file="grid.mat",
        time=dt.datetime(2016,7,8,0, 13,12),
        lats=[],
        lons=[],
        alts=[],
        to_file=None,
    ):
        mobj = {}
        gem = GEMINI2D(date_str, cfg, folder, param, grid_file)
        mobj["ne"] = gem.fetch_dataset_by_locations(time, lats, lons, alts)
        if to_file:
            savemat(to_file, mobj)
        return gem

    
if __name__ == "__main__":
    GEMINI2D.create_object("20160708", None)