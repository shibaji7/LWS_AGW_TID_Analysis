import datetime as dt
import os
import sys

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import rays
from fetchUtils import FetchData

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]
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

if 13 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19, 30), "fhe", 11, [21.6, 21.8]),
        rays.RayTraceObject(
            dt.datetime(2017, 5, 27, 19, 30),
            "fhe",
            11,
            [21.6, 21.8],
            run_name="May2017_gemini_control_cosmic",
            model_name="gemini",
        ),
    ]
    rp = rays.PlotRays(
        rtos[0], nrows=2, ncols=1, lw=0.5, arc=True, figsize=(8, 4), ylim=[-300, 600]
    )
    rp.lay_rays(
        text="(A) TID",
        zoomed_in=[[500, 1200], [150, 250]],
        add_cbar=False,
        ped_angles=[21.7],
    )
    rp.lay_rays(
        rto=rtos[1],
        text="(B) Control",
        ylabel="",
        # zoomed_in=[[500, 1200], [150, 250]],
        add_tag=False,
    )
    # rp.fig.subplots_adjust(hspace=1.2)
    rp.save(f"paper-mstid-rt/figures/Figure13.png")
    rp.save(f"paper-mstid-rt/pub_figures/Figure09.png")
    rp.close()
