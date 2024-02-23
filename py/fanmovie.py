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
from fanUtils import Fan, create_movie
import numpy as np
from plots import RTI

import os

foldername = "figures/FH"
rads = ["fhe", "fhw"]
dates = [
    dt.datetime(2017, 4, 21),
]

os.makedirs(f"{foldername}/fan/", exist_ok=True)
fds = dict()
colors, i = ["r", "b"], 0
for d in dates:
    for rad in rads:
        fd = FetchData.fetch(
            rad,
            [d, d + dt.timedelta(1)],
        )
        setattr(fd, "color", colors[i])
        fds[rad] = fd
        i += 1

scan_time = fds["fhe"].frame.scan_time.iloc[0]
for d in range(int(1440/(scan_time/60))):
    date = dates[0] + dt.timedelta(minutes=d*(scan_time/60))
    fname = f"figures/FH/fan/{'%04d'%d}.png"
    if not os.path.exists(fname):
        fan = Fan(rads, date)
        ax = fan.generate_fovs(fds)
        lons = np.arange(-110, -75, 5)
        lats = np.arange(44, 46, 0.1)
        fan.save(fname)
        fan.close()

cmd = f"ffmpeg -framerate 3 -pattern_type glob -i 'figures/FH/fan/*.png' -c:v libx264 figures/FH/fan/{dates[0].strftime('%Y-%m-%d')}.mp4"
os.system(cmd)

### Create RTI Plots
os.makedirs(f"{foldername}/rti/", exist_ok=True)
for rad, col in zip(rads, colors):
    beams = np.unique(fds[rad].frame.bmnum)
    for b in beams:
        fname = f"{foldername}/rti/{rad}_{b}_fov.png"
        if not os.path.exists(fname):
            fan = Fan(rads, dates[0])
            fan.overlay_fovs(rad, beams=[b])
            fan.save(fname)
            fan.close()
            rti = RTI(100, [dates[0], dates[0]+dt.timedelta(1)], f"{dates[0].strftime('%Y-%m-%d')} / {rad.upper()} / {b}")
            rti.addParamPlot(fds[rad].frame, b, "")
            rti.save(fname.replace("_fov", ""))
            rti.close()
        cmd = f"convert {fname} {fname.replace('_fov', '')} {fname.replace('_fov', '').replace('.png', '.pdf')}"
        os.system(cmd)
        os.remove(fname)
        os.remove(fname.replace("_fov", ""))
