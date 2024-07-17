#!/usr/bin/env python3

"""density_model.py: fetch GEMINI2D/3D IRI data"""

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


class IRI2D(object):
    
    def __init__(
        self, 
        time, 
        lats, 
        lons, 
        alts, 
        to_file=None, 
        cfg=None, 
        tid_prop = dict(),
        iri_version=20,
    ):
        """
        Setup all parameters
        """
        self.time = time
        self.lats = lats
        self.lons = lons
        self.alts = alts
        self.to_file = to_file
        self.tid_prop = tid_prop
        self.cfg = cfg
        self.alt_range = [alts[0], alts[-1], alts[1]-alts[0]]
        self.iri_version = iri_version
        self.load_data()
        self.create_TIDs()
        self.save_to_file()
        return
    
    def load_data(self, iri_version=None):
        """
        Load all IRI dataset
        """
        import iricore
        iri_version = iri_version if iri_version else self.iri_version
        self.edens = np.zeros((len(self.alts), len(self.lats)))
        for i in range(len(self.lats)):
            iriout = iricore.iri(
                self.time, 
                self.alt_range,
                self.lats[i],
                self.lons[i],
                iri_version,
            )
            self.edens[:,i] = iriout.edens * 1e-6
        return
    
    def create_TIDs(self):
        """
        Create TIDs in the system
        """
        dist = np.linspace(
            0, self.cfg.max_ground_range_km, self.cfg.number_of_ground_step_km
        )
        A = np.sin(2*np.pi*(dist - self.tid_prop["v_x"]*self.tid_prop["t"])/self.tid_prop["lamb_x"]) 
        # Rescale to 0-1
        A = (A - A.min())/(A.max() - A.min())
        # Resacle to A limits
        a_lim = self.tid_prop["a_lim"]
        A = (A * (a_lim[1]-a_lim[0])) + a_lim[0]
        for i in range(len(self.alts)):
            self.edens[i, :] = self.edens[i,:] * A
        return
    
    def save_to_file(self, to_file=None):
        """
        Save to file
        """
        to_file = to_file if to_file else self.to_file
        logger.info(f"IRI save to file: {to_file}")
        if to_file:
            mobj = {}
            mobj["ne"] = self.edens
            savemat(to_file, mobj)
        return


class GEMINI2D(object):
    
    def __init__(self, date_str, cfg, folder="dataset/GEMINI3D/", param="nsall", grid_file="coordinates.mat"):
        self.date_str = date_str
        self.cfg = cfg
        self.folder = f"{folder}{date_str}/"
        self.param = param
        self.grid_file = grid_file
        self.load_grid()
        self.load_datasets()
        return
    
    def load_grid(self):
        """
        Load grid file
        """
        self.grid = pd.DataFrame()
        with h5py.File(self.folder + self.grid_file, "r") as fkey:
            self.grid["glat"] = np.array(fkey.get("glat")[0]).tolist()
            self.grid["glon"] = np.mod( (fkey.get("glon")[0] + 180), 360 ) - 180
            self.grid["alt"] = fkey.get("galt")[0,:]
            self.grid["glat"] = fkey.get("glat")[0,:]
        return
    
    def load_dataset(self, day, fname=None, reset=True):
        """
        Load dataset
        """
        if reset:
            del self.dataset
            self.dataset = dict()
        if day not in list(self.dataset.keys()):
            fname = fname if fname else self.files[self.dates.index(day)]
            logger.info(f"Loading matlab file: {fname}")
            with h5py.File(fname, "r") as fkey:
                o = self.grid.copy()
                o[self.param] = fkey.get(self.param)[0,:]
                self.dataset[day] = o
                logger.info(f"Load dataset: {day}")
        return
    
    def load_datasets(self):
        """
        Load all dataset pointers
        """
        self.files = glob.glob(self.folder + self.date_str + "*")
        self.files.sort()
        self.dates = []
        self.dataset = dict()
        for fname in self.files:
            day = dt.datetime.strptime(fname.split("_")[0].split("/")[-1], "%Y%m%d")
            seconds = int(fname.split("_")[1].split(".")[0])
            day = day + dt.timedelta(seconds=seconds)
            self.dates.append(day)
        return
    
    def fetch_dataset_by_locations(self, time, lats, lons, alts, dlat=0.2, dlon=0.2, to_file=None):
        """
        Fetch data by lat/lon limits
        """
        self.load_dataset(time)
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
                method = "intp" if uf.alt.max()/1e3 > max(alts) else "extp"
                #print(np.array(uf.alt).max()/1e3, np.array(uf.alt).min()/1e3, max(alts), min(alts), method)
                out[:, ix] = utils.interpolate_by_altitude(
                        np.array(uf.alt)/1e3, alts, np.array(uf[self.param]),
                        self.cfg.scale, self.cfg.kind, method=method
                    ) * 1e-6
            ix += 1
        self.param_val = out
        if to_file:
            mobj = {}
            mobj["ne"] = out
            savemat(to_file, mobj)
        return out
    
    @staticmethod
    def create_object(
        date_str,
        cfg,
        folder="dataset/GEMINI3D/",
        param="nsall",
        grid_file="coordinates.mat",
        time=dt.datetime(2016,7,8,0,13,12),
        lats=[],
        lons=[],
        alts=[],
        dlat=0.2, 
        dlon=0.2,
        to_file=None,
    ):
        gem = GEMINI2D(date_str, cfg, folder, param, grid_file)
        gem.fetch_dataset_by_locations(time, lats, lons, alts, dlat, dlon, to_file)
        return gem
    
    @staticmethod
    def get_time_keys(
        date_str,
        folder="dataset/GEMINI3D/",
    ):
        files = glob.glob(folder + date_str + f"/{date_str}*")
        files.sort()
        dates = []
        for fname in files:
            day = dt.datetime.strptime(fname.split("_")[0].split("/")[-1], "%Y%m%d")
            seconds = int(fname.split("_")[1].split(".")[0])
            day = day + dt.timedelta(seconds=seconds)
            dates.append(day)
        return dates

    
if __name__ == "__main__":
    GEMINI2D.create_object("20160708", None)