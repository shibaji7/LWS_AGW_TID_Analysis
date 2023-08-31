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
from scipy.io import savemat

class Grid(object):
    """
    Generate 3D ionospheric grid
    """
    
    def __init__(
        self, 
        time, 
        lats, 
        lons, 
        alts,
        iri_version=20,
    ):
        return


class IRI3D(object):
    
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
