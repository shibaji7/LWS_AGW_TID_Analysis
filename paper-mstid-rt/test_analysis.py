import datetime as dt
import os
import sys

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import rays
from fetchUtils import FetchData
from fanUtils import Fan
import numpy as np
import pandas as pd
from rtiUtils import RTI

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
figures = [19]
rads = ["fhe", "fhw", "bks"]
dates = [
    dt.datetime(2017, 5, 27),
]

os.makedirs("paper-mstid-rt/figures", exist_ok=True)
os.makedirs("paper-mstid-rt/pub_figures", exist_ok=True)
fds = dict()
colors, i = ["r", "b", "m"], 0
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
    import matplotlib

    xlims = [dt.datetime(2017, 5, 27, 18), dt.datetime(2017, 5, 27, 21)]
    def get_curve(dfs, window_length=31, polyorder=3):
        dfs = dfs[
            (dfs.time >= xlims[0]) & (dfs.time <= xlims[1])
        ]
        dfs = dfs.dropna(subset=["p_l"])
        utime, lgate = (
            dfs.time.unique(), 
            dfs.groupby("time")["srange"].agg(np.nanmin)
        )
        from scipy.signal import savgol_filter, detrend
        smoothed = savgol_filter(lgate, window_length=window_length, polyorder=polyorder)
        dsmth = detrend(smoothed)
        return utime, lgate, smoothed, dsmth

    beam = 11
    DS = rays.get_datasets_by_beams(
        "fhe",
        [beam],
        start_time=dt.datetime(2017, 5, 27, 16),
        end_time=dt.datetime(2017, 5, 27, 21),
        limit_elvs=[18, 30],
        run_name="May2017_gemini_tid_cosmic2",
    )
    frame = fds["fhe"].frame.copy()
    frame = frame[frame.bmnum==beam]
    
    utime, model_gate, model_smoothed, model_dsmth  = get_curve(DS.copy(), 21, 3)

    dfs = frame.copy()
    dfs = dfs[dfs.srange>=1500]
    ugtime, obs_gate, obs_smoothed, obs_dsmth  = get_curve(dfs, 7, 1)
    model_tnum = mdates.date2num(utime)
    obs_tnum = mdates.date2num(ugtime)
    model_on_obs = np.interp(obs_tnum, model_tnum, model_dsmth)
    residuals = model_on_obs - obs_dsmth
    
    rti = RTI(
        100,
        xlims,
        fov=None,
        xlim=xlims,
        ylim=[180, 3000],
        fig_title="",
        num_subplots=3,
    )
    ax, _ = rti.addParamGsPlot(
        frame,
        beam,
        f"Rad: fhe-{beam}/ $f_0$= {(int((frame.tfreq.mean()/1e3)/0.5)+1)*0.5} MHz / 27 May 2017",
        vlim=[8, 16],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        major_loc=mdates.MinuteLocator(byminute=range(0, 60, 30)),
        minor_loc=mdates.MinuteLocator(byminute=range(0, 60, 15)),
    )
    ax.plot(ugtime, obs_smoothed, ls="-", lw=3, color="g")
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
        beam,
        "",
        vlim=[-90, -85],
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
        f"(B) Simulated",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ax.plot(utime, model_smoothed, ls="-", lw=3, color="m")

    ax = rti._add_axis()
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
    ax.xaxis.set_major_locator(mdates.MinuteLocator(byminute=range(0, 60, 30)))
    ax.xaxis.set_minor_locator(mdates.MinuteLocator(byminute=range(0, 60, 15)))
    ax.set_xlabel("Time UT", fontdict={"size": 12, "fontweight": "bold"})
    ax.set_xlim(
        [mdates.date2num(rti.drange[0]), mdates.date2num(rti.drange[1])]
    )
    ax.set_ylabel("Detrened ($\delta$) Slant Range, (km)")
    ax.set_ylim([-200, 200])
    ax.plot(utime, model_dsmth, ls="-", lw=1.0, color="m", label="Model")
    ax.plot(ugtime, obs_dsmth, ls="-", lw=1.0, color="g", label="Observations")
    ax.legend(loc="lower right", fontsize=8, framealpha=0.5)
    ax.text(
        0.05,
        0.95,
        f"(C)",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    pos = ax.get_position()
    inset_size = min(pos.width * 0.3, pos.height * 0.9)
    inset_left = min(pos.x1 + 0.15, 1.15 - inset_size)
    inset_bottom = pos.y0 + 0.5 * (pos.height - inset_size)
    inset_ax = ax.figure.add_axes([inset_left, inset_bottom, inset_size, inset_size])
    inset_ax.hist(residuals, color="0.6", edgecolor="k", linewidth=0.4)
    inset_ax.axvline(np.nanmean(residuals), color="r", ls="--", lw=0.8)
    inset_ax.set_xlabel(f"e($\delta$), km", fontsize=9, labelpad=1)
    inset_ax.set_ylabel("Counts", fontsize=9, labelpad=1)
    inset_ax.tick_params(axis="both", labelsize=6, length=2, pad=1)
    inset_ax.text(
        0.05,
        0.95,
        f"(D)",
        ha="left",
        va="center",
        transform=inset_ax.transAxes,
        fontdict=dict(size=9),
    )
    rti.save(f"paper-mstid-rt/figures/Figure19.png")
    # rti.save(f"paper-mstid-rt/pub_figures/Figure07.png")
    rti.save(f"paper-mstid-rt/pub_figures/FigureSp06.png")
    rti.close()

if 5 in figures:
    
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
        (frame.bmnum == 11)
        & (frame.time >= dt.datetime(2017, 5, 27, 16))
        & (frame.time <= dt.datetime(2017, 5, 27, 23))
    ].time.unique()
    dranges = [
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
        1700,
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
        # dsranges=(dranges, tnums),
    )
    rti.save("paper-mstid-rt/pub_figures/Figure05.png")
    rti.close()
