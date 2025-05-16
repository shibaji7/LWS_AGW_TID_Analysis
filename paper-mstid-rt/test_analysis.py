import datetime as dt
import os
import sys

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import rays
from fetchUtils import FetchData
from rtiUtils import RTI

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 14, 15, 16, 17, 18, 19]
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

if 19 in figures:
    import matplotlib.dates as mdates
    DS = rays.get_datasets_by_beams(
        "fhe", [11],
        start_time = dt.datetime(2017, 5, 27, 16), 
        end_time = dt.datetime(2017, 5, 27, 21),
        limit_elvs=[18, 30],
    )
    print(DS.head())
    frame = fds["fhe"].frame.copy()
    rti = RTI(
        100,
        [dt.datetime(2017, 5, 27, 18), dt.datetime(2017, 5, 27, 21)],
        fov=None,
        xlim=[dt.datetime(2017, 5, 27, 18), dt.datetime(2017, 5, 27, 21)],
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    ax, _ = rti.addParamGsPlot(
        frame,
        11,
        f"Rad: fhe-11/ $f_0$= {(int((frame.tfreq.mean()/1e3)/0.5)+1)*0.5} MHz / 27 May 2017",
        vlim=[12.5, 17.5],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        major_loc=mdates.MinuteLocator(byminute=range(0, 60, 30)),
        minor_loc=mdates.MinuteLocator(byminute=range(0, 60, 15)),
    )
    ax.text(
        0.05,
        0.95,
        f"(A) Observations",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ax, _ = rti.addParamPlot(
        DS,
        11,
        "",
        vlim=[-90, -85],
        zparam="p_l",
        xlabel="Time, UT",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        major_loc=mdates.MinuteLocator(byminute=range(0, 60, 30)),
        minor_loc=mdates.MinuteLocator(byminute=range(0, 60, 15)),
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Simulated",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti.save(f"paper-mstid-rt/figures/Figure19.png")
    rti.save(f"paper-mstid-rt/pub_figures/Figure08.png")
    rti.close()
