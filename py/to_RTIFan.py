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
from rtiUtils import RTI
from model_vheight import chisham_vhm, standard_vhm

# Initializations
CASES = [
    "RTI",
    "1DTS",
    "Fan",
    "vhm",
]
case = "RTI"
rads = ["fhw", "fhe"]
dates = [
    dt.datetime(2017, 5, 27),
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
                tec_mat_file=None,
                tec_param=None,
                power_vlim=[3, 30],
                tec_vlim=None,
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

if case == "vhm":
    for d in dates:
        for rad in rads:
            fd = FetchData.fetch(
                rad,
                [d, d + dt.timedelta(1)],
            )
            fd.frame["vheight"] = [standard_vhm(s, hop=1.0) for s in fd.frame.srange]
            for b in fd.frame.bmnum.unique():
                rt = RTI(
                    100,
                    [d, d + dt.timedelta(1)],
                    None,
                    [d + dt.timedelta(hours=14), d + dt.timedelta(1)],
                    f"{d.strftime('%Y-%m-%d')}/{rad}/{b}",
                    num_subplots=1
                )
                rt.addParamPlot(
                    fd.frame,
                    b,
                    "",
                    vlim=[0, 600],
                        cbar=True,
                    plot_fov=False,
                    zparam="vheight",
                    label=r"$H_{virtual}$ (km)"
                )
                file = (
                    tidUtils.get_folder(d) + f"/{rad}-{'%02d'%b}.png"
                )
                rt.save(file)
                rt.close()