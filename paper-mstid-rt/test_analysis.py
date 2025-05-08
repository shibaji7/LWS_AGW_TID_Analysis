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

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15]
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

if 15 in figures:
    sd = SolarDataset([dt.datetime(2017, 5, 27), dt.datetime(2017, 5, 28)])
    sd.create_stackplots("paper-mstid-rt/figures/Figure15.png")