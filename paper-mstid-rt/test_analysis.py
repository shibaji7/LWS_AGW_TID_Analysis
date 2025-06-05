import datetime as dt
import os
import sys

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import rays
from fetchUtils import FetchData

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
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

if 20 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [18, 30]),
    ]
    rp = rays.PlotRays(rtos[0], nrows=1, ncols=1, arc=True)
    rp.lay_rays(
        text="",
        tag_distance=1400,
        ped_angles=[22.7, 23, 21.3, 21.7],
    )
    rp.save(f"paper-mstid-rt/figures/Figure20.png")
    rp.close()
