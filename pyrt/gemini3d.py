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

class GEMINI2D(object):
    
    def __init__(self, cfg, year=2017, folder="", param=""):
        self.cfg = cfg
        self.year = year
        self.folder = folder
        self.param = param
        self.load_dataset()
        return
    
    def load_dataset(self):
        """
        Load all dataset available
        """
        return
    
    def fetch_dataset(self, time, latlim, lonlim):
        """
        Fetch data by lat/lon limits
        """
        return
    
     def fetch_dataset_by_locations(self, time, lats, lons, alts):
        """
        Fetch data by lat/lon limits
        """
        return
    
    @staticmethod
    def create_object(
        cfg,
        year=2017,
        folder="",
        param="",
        time=dt.datetime(2017, 8, 21, 17, 30),
        lats=[],
        lons=[],
        alts=[],
        to_file=None,
    ):
        mobj = {}
        gem = GEMINI3D(cfg, year, folder, param)
        mobj["ne"], _ = gem.fetch_dataset_by_locations(time, lats, lons, alts)
        if to_file:
            savemat(to_file, mobj)
        return gem
