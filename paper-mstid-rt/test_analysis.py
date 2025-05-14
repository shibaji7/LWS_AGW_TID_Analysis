import datetime as dt
import sys
import os

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import tidUtils
from fetchUtils import FetchData
from fanUtils import Fan
from rtiUtils import RTI
import cartopy.crs as ccrs
import numpy as np
import rays
from solar import SolarDataset

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18]
# figures = [3]
rads = ["fhw", "fhe"]
dates = [
    dt.datetime(2017, 5, 27),
]

os.makedirs("paper-mstid-rt/figures", exist_ok=True)
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

if 14 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [18, 30]),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 17), "fhe", 11, [18, 30]),
    ]
    rp = rays.PlotChannels(rtos[0], nrows=2, ncols=1, xlim=[18, 24], ylim=[500, 2000])
    rp.lay_rays(add_cbar=False)
    # rp.lay_rays(text="(A)", xlabel="", zoomed_in=[[500, 1200], [150, 250]], add_cbar=False)
    rp.lay_rays(rto=rtos[1], text="(B)", add_tag=False)
    rp.save(f"paper-mstid-rt/figures/Figure14.png")
    rp.close()