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

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 25]
figures = [25]
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

if 25 in figures:
    data = fds["fhe"].frame.copy()
    print(data.srange.head())
    print(data.columns)
    from model_vheight import chisham_vhm
    vheights = [chisham_vhm(s) for s in data.srange]
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(4,4))
    ax.hist(vheights, bins=30, color='skyblue', histtype="step")
    ax.set_title("Virtual Heights, Chisham[2008] model")
    ax.set_xlabel("Virtual Heights, km")
    ax.set_ylabel("Counts")
    ax.set_yscale("log")
    fig.savefig("paper-mstid-rt/figures/Figure25.png")
    # print(vheights)

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

    beam = 19
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
    ax.plot(ugtime, obs_smoothed, ls="-", lw=3, color="k")
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
        label=r"$P_r$, dB",
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
    ax.plot(utime, model_smoothed, ls="-", lw=3, color="g")

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
    ax.plot(utime, model_dsmth, ls="-", lw=1.0, color="g", label="Model")
    ax.plot(ugtime, obs_dsmth, ls="-", lw=1.0, color="k", label="Observations")
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
    # rti.save(f"paper-mstid-rt/figures/Figure19.png")
    # rti.save(f"paper-mstid-rt/pub_figures/Figure07.png")
    rti.save(f"paper-mstid-rt/pub_figures/FigureSp06.png")
    rti.close()

if 5 in figures:
    import matplotlib.dates as mdates
    from rad_fov import CalcFov
    import pydarn
    from geopy.distance import great_circle
    
    rti = RTI(
        100,
        [dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 27, 22)],
        fov=None,
        xlim=[dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 27, 22)],
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    frame = fds["fhe"].frame.copy()
    hdw = pydarn.read_hdw_file("fhe", dt.datetime(2017, 5, 27))
    fov = CalcFov(hdw=hdw, model="GS", fov_dir="front")
    station = (hdw.geographic.lat, hdw.geographic.lon)
    f11 = frame[frame.bmnum == 11].copy()
    lat, lon = fov.latCenter[10, :], fov.lonCenter[10, :]
    f11["lat"] = f11.slist.apply(lambda x: lat[int(x)-1])
    f11["lon"] = f11.slist.apply(lambda x: lon[int(x)-1])
    # Use great_circle distance as shared Y metric for both radar and TEC
    f11["gsMapped"] = f11.apply(
        lambda x: great_circle((x.lat, x.lon), station).kilometers
        if not np.isnan(x.lat) else np.nan, axis=1
    )

    ax, _ = rti.addParamGsPlot(
        f11,
        11,
        f"Rad: fhe / $f_0$= {(int((frame.tfreq.mean()/1e3)/0.5)+1)*0.5} MHz / 27 May 2017",
        vlim=[12.5, 17.5],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        yscale="gsMapped",
    )
    ax.set_ylim([100, 1500])
    ax.set_yticks([100, 500, 1000, 1500])
    ax.set_ylabel(r"GS Mapped, $(\theta,\phi)$", fontdict={"size": 12, "fontweight": "bold"})
    f11_valid = f11[(f11.gsMapped >= 0) & (f11.gsMapped <= 1500)].dropna(subset=["lat"])
    ytick_labels, yticks = [], []
    for label in ax.get_yticklabels():
        label_value = int(label.get_text())
        ilab = (f11_valid["gsMapped"] - label_value).abs().idxmin()
        ytick_labels.append(f"${f11_valid.loc[ilab, 'lat']:.1f}^\circ$\n${f11_valid.loc[ilab, 'lon']:.1f}^\circ$")
        yticks.append(f11_valid.loc[ilab, "gsMapped"])
    ax.set_yticks(yticks)
    ax.tick_params(axis="y", labelrotation=90)
    ax.set_yticklabels(ytick_labels, fontdict={"size": 8})

    import sys
    sys.path.append("paper-mstid-rt/")
    from nexrad_utils import load_fulltimedata
    df = load_fulltimedata(
        '2017-05-27', beam="9to13",
        mat_dir='/home/chakras4/Research/Individual_Studies/LWS_AGW_TID_Analysis/data',
        start_time='16:00', end_time='22:00',
        elevation_cutoff=20.0,
    )
    df["mappedGS"] = df.apply(
        lambda x: great_circle(
            (x.ipp_lat, x.ipp_lon), (hdw.geographic.lat, hdw.geographic.lon)
        ).kilometers, 
        axis=1
    )
    import tidUtils
    import matplotlib.dates as mdates
    df["mappedGS"] = df["mappedGS"].round(-1)   # 10-km bins -> ~141 unique Y values
    X, Y, Z = tidUtils.get_gridded_parameters(
        df, xparam="time", yparam="mappedGS", zparam="vtec_bp5_40", rounding=False
    )
    # X stays as datetime64 — ax has a DateConverter from addParamGsPlot's scatter.
    # Passing mdates floats would double-convert them → off-screen. Use datetime64 directly.
    tec_vals = Z.compressed() if np.ma.is_masked(Z) else Z[np.isfinite(Z)]
    tec_lim = max(abs(np.nanpercentile(tec_vals, 2)), abs(np.nanpercentile(tec_vals, 98)))
    pm = ax.pcolormesh(
        X, Y, Z.T,
        # cmap="RdBu_r",
        vmin=-0.03, vmax=0.03,
        shading="auto",
        alpha=0.7,
        lw=0.01,
        edgecolors="None",
        cmap="Blues_r",
        snap=True,
    )
    ax.text(
        0.95,
        1.05,
        f"(A) Beam: 11",
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )


    f19 = frame[frame.bmnum == 19].copy()
    lat, lon = fov.latCenter[18, :], fov.lonCenter[18, :]
    f19["lat"] = f19.slist.apply(lambda x: lat[int(x)-1])
    f19["lon"] = f19.slist.apply(lambda x: lon[int(x)-1])
    f19["gsMapped"] = f19.apply(
        lambda x: great_circle((x.lat, x.lon), station).kilometers
        if not np.isnan(x.lat) else np.nan, axis=1
    )
    ax, _ = rti.addParamGsPlot(
        f19,
        19,
        "",
        vlim=[12.5, 17.5],
        zparam="p_l",
        xlabel="Time, UT",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=False,
        yscale="gsMapped",
    )
    ax.set_ylim([100, 1500])
    ax.set_yticks([100, 500, 1000, 1500])
    ax.set_ylabel(r"GS Mapped, $(\theta,\phi)$", fontdict={"size": 12, "fontweight": "bold"})
    f19_valid = f19[(f19.gsMapped >= 0) & (f19.gsMapped <= 1500)].dropna(subset=["lat"])
    ytick_labels, yticks = [], []
    for label in ax.get_yticklabels():
        label_value = int(label.get_text())
        ilab = (f19_valid["gsMapped"] - label_value).abs().idxmin()
        ytick_labels.append(f"${f19_valid.loc[ilab, 'lat']:.1f}^\circ$\n${f19_valid.loc[ilab, 'lon']:.1f}^\circ$")
        yticks.append(f19_valid.loc[ilab, "gsMapped"])
    ax.set_yticks(yticks)
    ax.tick_params(axis="y", labelrotation=90)
    ax.set_yticklabels(ytick_labels, fontdict={"size": 8})

    df = load_fulltimedata(
        '2017-05-27', beam="17to21",
        mat_dir='/home/chakras4/Research/Individual_Studies/LWS_AGW_TID_Analysis/data',
        start_time='16:00', end_time='22:00',
        elevation_cutoff=20.0,
    )
    df["mappedGS"] = df.apply(
        lambda x: great_circle(
            (x.ipp_lat, x.ipp_lon), (hdw.geographic.lat, hdw.geographic.lon)
        ).kilometers, 
        axis=1
    )
    df["mappedGS"] = df["mappedGS"].round(-1) 
    X, Y, Z = tidUtils.get_gridded_parameters(
        df, xparam="time", yparam="mappedGS", zparam="vtec_bp5_40", rounding=False
    )
    pm = ax.pcolormesh(
        X, Y, Z.T,
        # cmap="RdBu_r",
        vmin=-0.03, vmax=0.03,
        shading="auto",
        alpha=0.7,
        lw=0.01,
        edgecolors="None",
        cmap="Blues_r",
        snap=True,
    )
    ax.text(
        0.95,
        1.05,
        f"(B) Beam: 19",
        ha="right",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti._add_colorbar(rti.fig, ax, pm, label=r"dvTEC$_{5-40}$, TECU", xOff=0.0)

    rti.save("paper-mstid-rt/pub_figures/Figure05.png")
    rti.close()
