#!/usr/bin/env python3

"""density_model_3D.py: fetch GEMINI3D and IRI3D data"""

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
import os
from scipy.io import savemat


from multiprocessing import Pool

class Iono3D(object):
    
    def __init__(
        self, 
        time, 
        bearing_obj, 
        to_file=None, 
        cfg=None, 
        tid_prop = dict(
            lamb_x = 1000,
            v_x = 0.1,
            a_lim = [0.1, 1],
            t = 0,
            source="center",
        ),
        iri_version=20,
    ):
        """
        Setup all parameters
        """
        self.time = time
        self.bearing_obj = bearing_obj
        self.to_file = to_file
        self.tid_prop = tid_prop
        self.cfg = cfg
        self.iri_version = iri_version
        self.iono = {}
        self.create_grids()
        if not os.path.exists(self.to_file):
            self.load_iri_data()
            self.save_to_file()
        self.create_TIDs()
        return
    
    def create_grids(self):
        lats = np.arange(
            self.bearing_obj["lat_start"], 
            self.bearing_obj["lat_start"] + self.bearing_obj["lat_inc"]*int(self.bearing_obj["num_lat"]), 
            self.bearing_obj["lat_inc"]
        )
        self.lats = lats[:int(self.bearing_obj["num_lat"])]
        lons = np.arange(
            self.bearing_obj["lon_start"], 
            self.bearing_obj["lon_start"] + self.bearing_obj["lon_inc"]*self.bearing_obj["num_lon"], 
            self.bearing_obj["lon_inc"]
        )
        self.lons = lons[:int(self.bearing_obj["num_lon"])]
        alts = np.arange(
            self.bearing_obj["ht_start"], 
            self.bearing_obj["ht_start"] + self.bearing_obj["ht_inc"]*self.bearing_obj["num_ht"], 
            self.bearing_obj["ht_inc"]
        )
        self.alts = alts[:int(self.bearing_obj["num_ht"])]
        self.edens = (
            np.zeros((len(self.lats), len(self.lons), len(self.alts)))*np.nan
        )
        if self.tid_prop["source"] == "center":
            self.tid_source = (np.mean(lats), np.mean(lons))
        return
    
    def load_iri_data(self, iri_version=None):
        """
        Load all IRI dataset
        """
        logger.info(f"Started IRI calculations of shape: {self.edens.shape}")
        alt_range = [self.alts[0], self.alts[-1], self.alts[1]-self.alts[0]]
        import iricore
        iri_version = iri_version if iri_version else self.iri_version
        for i, lat in enumerate(self.lats):
            for j, lon in enumerate(self.lons):
                iriout = iricore.iri(
                    self.time, 
                    alt_range,
                    lat,
                    lon,
                    iri_version,
                )
                self.edens[i,j,:] = iriout.edens * 1e-6
        logger.info(f"Stopped populating IRI")
        (self.iono["iono_en_grid"], self.iono["iono_en_grid_5"]) = (
            self.edens, self.edens
        )
        return
    
    def create_TIDs(self):
        """
        Create TIDs in the system
        """
        logger.info(f"Create TIDs: {self.tid_prop['t']}")
        import geopy
        dist = np.zeros((len(self.lats), len(self.lons)))*np.nan
        for i, lat in enumerate(self.lats):
            for j, lon in enumerate(self.lons):
                dist[i,j] = geopy.distance.geodesic(self.tid_source, (lat, lon)).km
        A = np.sin(2*np.pi*(dist - self.tid_prop["v_x"]*self.tid_prop["t"])/self.tid_prop["lamb_x"]) 
        # Rescale to 0-1
        A = (A - A.min())/(A.max() - A.min())
        # Resacle to A limits
        a_lim = self.tid_prop["a_lim"]
        A = (A * (a_lim[1]-a_lim[0])) + a_lim[0]
        for i in range(len(self.alts)):
            self.edens[:,:,i] = self.edens[:,:,i] * A
        return
    
    def save_to_file(self, to_file=None):
        """
        Save to file
        """
        to_file = to_file if to_file else self.to_file
        logger.info(f"Ionosphere save to file: {to_file}")
        if to_file:
            savemat(to_file, self.iono)
        return
