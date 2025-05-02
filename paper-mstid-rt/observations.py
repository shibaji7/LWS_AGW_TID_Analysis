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

figures = [1, 3, 4, 7, 8]
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


if 1 in figures:
    fan = Fan(
        rads,
        dates[0],
        txt_coord=True,
        cbar=False,
    )
    ax = fan.overlay_fovs("fhe", beams=[3, 11], col="b")
    fan.overlay_fovs("fhw", ax=ax, col="r")
    # ax.overlay_station("alp", 45.0617, -83.4328)
    fan.save("paper-mstid-rt/figures/Figure01.png")
    fan.close()

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
        3,
        "",
        vlim=[10, 25],
        zparam="p_l",
        label="Power, dB",
        cmap="plasma",
        cbar=False,
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


# Similar figure for Simulation for comparison
if 4 in figures:
    fan = Fan(
        rads,
        dates[0],
        nrows=3,
        ncols=2,
        add_text=True,
    )
    gates = [38, 35, 29, 33, 37, 35]
    for i, d in enumerate(
        (np.array([17.25, 17.5, 19, 19.25, 20.0, 21.0]) * 60).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
        ax = fan.generate_fovs(
            fds, beams={"fhe": [11, 3], "fhw": []}, discreat={"fhe": 27, "fhw": 27}
        )
        fan.add_arc_fov("fhe", ax=ax, maxGate=27, lineColor="k", lineWidth=0.6)
        fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="m", lineWidth=1.2)
        ax.text(
            0.05,
            0.95,
            f"({chr(65+i)})",
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict=dict(size="xx-small"),
        )
    fan.save(f"paper-mstid-rt/figures/Figure04.png")
    fan.close()

if 5 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 17, 30), "fhe", 11),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 18, 15), "fhe", 11),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19, 30), "fhe", 11),
    ]
    rp = rays.PlotRays(rtos[0], nrows=2, ncols=2)
    rp.lay_rays(xlabel="", add_cbar=False, text="(A)")
    rp.lay_rays(
        ylabel="", xlabel="", rto=rtos[1], add_cbar=False, add_tag=False, text="(B)"
    )
    rp.lay_rays(rto=rtos[2], add_cbar=False, add_tag=False, text="(C)")
    rp.lay_rays(ylabel="", rto=rtos[3], add_tag=False, text="(D)")
    rp.save(f"paper-mstid-rt/figures/Figure05.png")
    rp.close()

if 7 in figures:
    ddates = [dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 27, 21)]
    rti = RTI(
        100,
        ddates,
        fov=None,
        xlim=ddates,
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    ax, _ = rti.addParamPlot(
        fds["fhe"].frame,
        11,
        f"Rad: fhe / Beam: 11 / $f_0$= {(fds['fhe'].frame.tfreq.mean()/1e3).round(1)} MHz / 27 May 2017",
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
        f"(A) Observation",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ds = rays.get_datasets_by_beams("fhe", [11], ddates[0], ddates[1])
    ax, _ = rti.addParamPlot(
        ds,
        11,
        "",
        vlim=[-90, -75],
        zparam="p_l",
        label=r"$P_r$, dB",
        cmap="plasma",
        cbar=True,
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Simulation",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti.save(f"paper-mstid-rt/figures/Figure07.png")
    rti.close()

if 8 in figures:
    fan = Fan(
        ["fhe"],
        dates[0],
        nrows=1,
        ncols=2,
        add_text=True,
        extent=[-105, -70, 35, 60],
    )
    date = dates[0] + dt.timedelta(minutes=int(60 * 19.5))
    ds = rays.get_datasets_by_beams("fhe", None, date, date + dt.timedelta(minutes=1))
    fan.date, fan.txt_coord, fan.cbar = date, True, True
    ax = fan.generate_fovs(dict(fhe=fds["fhe"]), beams={"fhe": []})
    ax.text(
        0.05,
        0.95,
        f"(A) Observations",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="xx-small"),
    )
    fan.date, fan.txt_coord, fan.cbar = date, False, True
    ax = fan.overlay_simulation_fovs("fhe", ds, beams={}, gate_filter=[])
    ax.text(
        0.05,
        0.95,
        f"(B) Simulations",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="xx-small"),
    )
    fan.save(f"paper-mstid-rt/figures/Figure08.png")
    fan.close()

if 9 in figures:
    fan = Fan(
        rads,
        dates[0],
        nrows=3,
        ncols=2,
        add_text=True,
    )
    gates = [35, 35, 29, 33, 37, 35]
    for i, d in enumerate(
        (
            np.array(
                [
                    17.5,
                    18,
                    19,
                    20.0,
                ]
            )
            * 60
        ).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
        ds = rays.get_datasets_by_beams(
            "fhe", None, date, date + dt.timedelta(minutes=1)
        )
        ax = fan.overlay_simulation_fovs("fhe", ds, beams={}, gate_filter=[])
        fan.add_arc_fov("fhe", ax=ax, maxGate=25, lineColor="k", lineWidth=0.6)
        fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="k", lineWidth=0.6)
        ax.text(
            0.05,
            0.95,
            f"({chr(65+i)})",
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict=dict(size="xx-small"),
        )
    fan.save(f"paper-mstid-rt/figures/Figure09.png")
    fan.close()


if 10 in figures:
    ddates = [dt.datetime(2017, 5, 27, 16), dt.datetime(2017, 5, 27, 21)]
    rti = RTI(
        100,
        ddates,
        fov=None,
        xlim=ddates,
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    ds = rays.get_datasets_by_beams("fhe", [11], ddates[0], ddates[1])
    ax, _ = rti.addParamPlot(
        ds,
        11,
        f"Rad: fhe / Beam: 11 / $f_0$= {(fds['fhe'].frame.tfreq.mean()/1e3).round(1)} MHz / 27 May 2017",
        vlim=[-90, -75],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
    )
    ax.text(
        0.05,
        0.95,
        f"(A) Simulation / R12=12",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    ds = rays.get_datasets_by_beams(
        "fhe", [11], ddates[0], ddates[1], run_name="May2017_gemini_tid_cosmic2_R12.100"
    )
    ax, _ = rti.addParamPlot(
        ds,
        11,
        "",
        vlim=[-90, -75],
        zparam="p_l",
        label=r"$P_r$, dB",
        cmap="plasma",
        cbar=True,
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Simulation / R12=100",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size="small"),
    )
    rti.save(f"paper-mstid-rt/figures/Figure10.png")
    rti.close()

if 11 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 17, 30), "fhe", 11, [18, 30]),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 18, 15), "fhe", 11, [18, 30]),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [18, 30]),
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19, 30), "fhe", 11, [18, 30]),
    ]
    rp = rays.PlotRays(rtos[0], nrows=2, ncols=2)
    rp.lay_rays(xlabel="", add_cbar=False, text="(A)", tag_distance=1400)
    rp.lay_rays(
        ylabel="",
        xlabel="",
        rto=rtos[1],
        add_cbar=False,
        add_tag=False,
        text="(B)",
        tag_distance=1400,
    )
    rp.lay_rays(
        rto=rtos[2], add_cbar=False, add_tag=False, text="(C)", tag_distance=1400
    )
    rp.lay_rays(ylabel="", rto=rtos[3], add_tag=False, text="(D)", tag_distance=1400)
    rp.save(f"paper-mstid-rt/figures/Figure11.png")
    rp.close()

if 12 in figures:
    fan = Fan(
        rads,
        dates[0],
        nrows=3,
        ncols=2,
        add_text=True,
    )
    gates = [35, 35, 29, 33, 37, 35]
    for i, d in enumerate(
        (
            np.array(
                [
                    17.5,
                    18,
                    19,
                    20.0,
                ]
            )
            * 60
        ).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
        ds = rays.get_datasets_by_beams(
            "fhe", None, date, date + dt.timedelta(minutes=1), [18, 30]
        )
        ax = fan.overlay_simulation_fovs("fhe", ds, beams={}, gate_filter=[])
        fan.add_arc_fov("fhe", ax=ax, maxGate=25, lineColor="k", lineWidth=0.6)
        fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="k", lineWidth=0.6)
        ax.text(
            0.05,
            0.95,
            f"({chr(65+i)})",
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict=dict(size="xx-small"),
        )
    fan.save(f"paper-mstid-rt/figures/Figure12.png")
    fan.close()
