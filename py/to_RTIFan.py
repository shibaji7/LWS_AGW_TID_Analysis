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

# Initializations
CASES = [
    "RTI",
    "1DTS",
    "Fan",
]
case = "1DTS"
rads = ["fhw"]
dates = [
    dt.datetime(2022, 12, 20),
]
bm_gate_cells = [(7, 25)]
fdMap = {}


# Conduct RTI analysis
if case == "RTI":
    for d in dates:
        for rad in rads:
            fd = FetchData.fetch(
                rad,
                [d, d + dt.timedelta(1)],
            )
            fd.plot_RTI(
                tec_mat_file=f"data/{d.strftime('%Y-%m-%d')}/{rad}geom.mat",
                tec_param="cdvTECgrid2",
                power_vlim=[20, 50],
                tec_vlim=[-0.3, 0.3],
            )

# Conduct 1D-Timeseries analysis
if case == "1DTS":
    for d in dates:
        for rad in rads:
            fd = FetchData.fetch(
                rad,
                [d, d + dt.timedelta(1)],
            )
            for bm_gate in bm_gate_cells:
                fd.TS1D(
                    bm_gate[0],
                    bm_gate[1],
                    tec_mat_file=f"data/{d.strftime('%Y-%m-%d')}/{rad}geom.mat",
                    tec_param="cdvTECgrid2",
                    power_vlim=[0, 40],
                    tec_vlim=[-0.3, 0.3],
                )
