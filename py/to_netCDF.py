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
import tidUtils
from fetchUtils import FetchData

rads = [
    "fhe",
    "fhw",
    "cve",
    "cvw"
]
fdMap = {}
dates = [
    dt.datetime(2017, 9, 7),
    dt.datetime(2017, 9, 8),
]
for d in dates:
    for rad in rads:
        fd = FetchData.fetch(
            rad, [d, d + dt.timedelta(1)], med_filter={"cpu": 4, "thresh": 0.7}
        )
        fd.to_geom()
