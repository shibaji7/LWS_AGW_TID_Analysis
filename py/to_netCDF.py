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

rads = [
    "bks",
    "fhe",
    "fhw",
]
fdMap = {}
dates = [
#     dt.datetime(2022, 12, 18), 
#     dt.datetime(2022, 12, 19),
#     dt.datetime(2022, 12, 20),
    dt.datetime(2022, 12, 21), 
#     dt.datetime(2022, 12, 22),
#     dt.datetime(2022, 12, 23),
#     dt.datetime(2022, 12, 24)
]
for d in dates:
    fdMap[d] = {}
    for rad in rads:
        fd = FetchData.fetch(
            rad,
            [d, d + dt.timedelta(1)],
            med_filter={"cpu": 4, "thresh":.7}
        )
        fd.plot_RTI(angle_th=110.)
        fdMap[d][r] = fd
            