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
figures = [3]
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

if 3 in figures:
    rti = RTI(
        100,
        [dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 28)],
        fov=None,
        xlim=[dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 28)],
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    frame = fds["fhe"].frame.copy()
    tnums = frame[
        (frame.bmnum==11) 
        & (frame.time>=dt.datetime(2017, 5, 27, 16))
        & (frame.time<=dt.datetime(2017, 5, 27, 23))
    ].time.unique()
    dranges = [
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700, 1700,
        1700, 1700, 1700, 1700, 1700, 
    ]
    ax, _ = rti.addParamGsPlot(
        frame,
        11,
        f"Rad: fhe / $f_0$= {(int((frame.tfreq.mean()/1e3)/0.5)+1)*0.5} MHz / 27 May 2017",
        vlim=[12.5, 17.5],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        dsranges=(dranges, tnums)
    )
    ax.text(
        0.05,
        0.95,
        f"(A) Beam: 11",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ax, _ = rti.addParamGsPlot(
        frame,
        3,
        "",
        vlim=[12.5, 17.5],
        zparam="p_l",
        label="Power, dB",
        cmap="plasma",
        cbar=False,
        dsranges=(dranges, tnums)
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Beam: 03",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti.save(f"paper-mstid-rt/figures/Figure03.png")
    rti.close()

if 17 in figures:
    date = dates[0] + dt.timedelta(minutes=int(60 * 18))
    import matplotlib.pyplot as plt

    plt.rcParams.update({'font.size': 8})
    fan = Fan(
        ["fhe"],
        date,
        nrows=1,
        ncols=2,
        add_text=True,
        extent=[-105, -70, 35, 60],
    )

    fan.date, fan.txt_coord, fan.cbar = date, True, True
    ax = fan.generate_fovs(dict(fhe=fds["fhe"]), beams={"fhe": []}, discreat={"fhe": 0})
    fan.date, fan.txt_coord, fan.cbar = dates[0] + dt.timedelta(minutes=int(60 * 20)), False, False
    ax = fan.generate_fovs(dict(fhe=fds["fhe"]), beams={"fhe": []}, discreat={"fhe": 0})
    fan.save("paper-mstid-rt/figures/Figure_Analysis_FoV_18.png")
    fan.close()

if 18 in figures:
    rti = RTI(
        100,
        [dt.datetime(2017, 5, 27), dt.datetime(2017, 5, 28)],
        fov=None,
        xlim=[dt.datetime(2017, 5, 27), dt.datetime(2017, 5, 28)],
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    frame = fds["fhe"].frame.copy()
    ax, _ = rti.addParamPlot(
        frame,
        11,
        f"Rad: fhe / $f_0$= {(frame.tfreq.mean()/1e3).round(1)} MHz / 27 May 2017",
        vlim=[10, 25],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
    )
    ax.text(
        0.05,
        0.95,
        f"(A) Beam: 11",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ax, _ = rti.addParamPlot(
        frame,
        11,
        "",
        vlim=[-50, 50],
        zparam="v",
        label="Velocity, m/s",
        cmap="Spectral",
        cbar=True,
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Beam: 11",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti.save("paper-mstid-rt/figures/Figure_Analysis_RTI.png")
    rti.close()
