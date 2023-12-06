import datetime as dt
import sys
import os

sys.path.extend(["../py/txUtils/", "../py/tid/", "../py/davitPy/"])
import tidUtils
from fetchUtils import FetchData
from fanUtils import Fan
import cartopy.crs as ccrs
import numpy as np

rads = ["fhw", "fhe"]
dates = [
    dt.datetime(2023, 6, 29),
]

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

for d in range(1440):
    date = dates[0] + dt.timedelta(minutes=d)
    fname = f"figures/{date.strftime('%Y-%m-%d-%H-%M')}.png"
    if not os.path.exists(fname):
        fan = Fan(
            rads,
            date,
        )
        ax = fan.generate_fovs(fds)
        lons = np.arange(-110, -75, 5)
        lats = np.arange(44, 46, 0.1)
        fan.save(fname)
        fan.close()