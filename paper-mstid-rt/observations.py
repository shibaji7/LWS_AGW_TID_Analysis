import datetime as dt
import sys
import os

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/"])
import tidUtils
from fetchUtils import FetchData
from fanUtils import Fan
from rtiUtils import RTI
import cartopy.crs as ccrs
import numpy as np

rads = ["fhw", "fhe"]
dates = [
    dt.datetime(2017, 5, 27),
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

rti = RTI(
    100,
    [dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 28)],
    fov=None,
    xlim=[dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 28)],
    ylim=[180, 3000],
    fig_title=None,
    num_subplots=2
)
ax, _ = rti.addParamPlot(
    fds["fhe"].frame, 11, 
    f"Rad: fhe / Beam: 11 / $f_0$= {(fds['fhe'].frame.tfreq.mean()/1e3).round(1)} MHz / 27 May 2017", 
    vlim=[10, 25], zparam="p_l",
    label="Power, dB", cmap="plasma", xlabel=""
)
ax.text(
    0.05, 0.95, f"(A)", 
    ha="left", va="center", 
    transform=ax.transAxes, fontdict=dict(size="small")
)
ax, _ = rti.addParamPlot(
    fds["fhe"].frame, 11, "",
    vlim=[-30, 30],
    cmap="Spectral",
)
ax.text(
    0.05, 0.95, f"(B)", 
    ha="left", va="center", 
    transform=ax.transAxes, fontdict=dict(size="small")
)
rti.save(f"figures/Figure5.png")
rti.close()


fan = Fan(
    rads, dates[0],
    txt_coord=True,
    cbar=False,
)
ax = fan.overlay_fovs("fhe", beams=[3, 11], col="b")
fan.overlay_fovs("fhw", ax=ax, col="r")
ax.overlay_station("alp", 45.0617, -83.4328)
fan.save("figures/fov.png")
fan.close()

fan = Fan(
    rads,
    dates[0],
    nrows=3, ncols=2,
    add_text=True,
)
for i, d in enumerate(range(60*18,60*21,30)):
    date = dates[0] + dt.timedelta(minutes=d)
    fan.date, fan.txt_coord, fan.cbar = date, (d==60*18), (d==60*20)
    ax = fan.generate_fovs(fds, beams={"fhe":[11,3], "fhw":[]})
    ax.text(
        0.05, 0.95, f"({chr(65+i)})", 
        ha="left", va="center", 
        transform=ax.transAxes, fontdict=dict(size="xx-small")
    )
fan.save(f"figures/Figure4.png")
fan.close()

