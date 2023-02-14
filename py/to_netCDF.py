#!/usr/bin/env python

"""to_netCDF.py: utility module to stored as a netCDF file."""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import sys

sys.path.extend(["py/txUtils/", "py/tid/", "py/davitPy/"])
from fetchUtils import FetchData
import tidUtils

rads = ["fhe", "fhw", "bks"]
fdMap = {}
dates = [dt.datetime(2022, 12, 19), dt.datetime(2022, 12, 20), dt.datetime(2022, 12, 21), dt.datetime(2022, 12, 22), dt.datetime(2022, 12, 23), dt.datetime(2022, 12, 24)]
for d in dates:
    fdMap[d] = FetchData.fetch(rads, [d, d + dt.timedelta(1)], to_netcdf=True)