import datetime as dt
import os
import sys
from pathlib import Path

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import numpy as np
import rays
from fanUtils import Fan
from fetchUtils import FetchData
from rtiUtils import RTI
from solar import SolarDataset

figures = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13]
figures = [2]
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
    ax = fan.overlay_fovs("fhe", beams=[3, 11], col="r")
    fan.overlay_fovs("fhw", ax=ax, col="b")
    # ax.overlay_station("alp", 45.0617, -83.4328)
    fan.save("paper-mstid-rt/figures/Figure01.png")
    fan.close()

if 2 in figures:
    date = dates[0] + dt.timedelta(minutes=int(60 * 19.5))
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
    fan.save("paper-mstid-rt/figures/Figure02.png")
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
        dsranges=(dranges, tnums),
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
        dsranges=(dranges, tnums),
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
    gates = [32, 32, 29, 33, 34, 35]
    for i, d in enumerate(
        # (np.array([17.5, 18, 19, 19.5, 20.0, 20.5]) * 60).astype(int)
        (np.array([18+(4/6),18+(5/6), 19, 19+(1/6), 19+(2/6), 19+(3/6)]) * 60).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        # fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
        fan.date, fan.txt_coord, fan.cbar = date, i==0, i==4
        ax = fan.generate_fovs(
            fds, beams={"fhe": [11, 3]}, discreat={"fhe": 27}
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
    fan.save(f"paper-mstid-rt/figures/Figure03.png")
    fan.save(f"paper-mstid-rt/pub_figures/Figure03.png")
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
    ax = fan.generate_fovs(dict(fhe=fds["fhe"]), beams={"fhe": []}, discreat={"fhe": 0})
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
    rp = rays.PlotRays(rtos[0], nrows=4, ncols=1, arc=True)
    rp.lay_rays(
        xlabel="",
        add_cbar=False,
        text="(A)",
        tag_distance=1400,
        ped_angles=[
            24.3,
        ],
    )
    rp.lay_rays(
        rto=rtos[1],
        add_cbar=False,
        add_tag=False,
        text="(B)",
        tag_distance=1400,
        ped_angles=[22.6],
    )
    rp.lay_rays(
        rto=rtos[2],
        add_cbar=False,
        add_tag=False,
        text="(C)",
        tag_distance=1400,
        ped_angles=[22.7, 23, 21.3, 21.7],
    )
    rp.lay_rays(
        rto=rtos[3],
        add_tag=False,
        text="(D)",
        tag_distance=1400,
        ped_angles=[21.7, 23.4],
    )
    rp.save(f"paper-mstid-rt/figures/Figure11.png")
    rp.save(f"paper-mstid-rt/pub_figures/Figure05.png")
    rp.close()

if 12 in figures:
    fan = Fan(
        rads,
        dates[0],
        nrows=3,
        ncols=2,
        add_text=True,
    )
    # gates = [38, 35, 29, 33, 37, 35]
    gates = [32, 32, 29, 33, 34, 35]
    mgates = [31, 31, 28, 28, 26, 29]
    for i, d in enumerate(
        # (np.array([17.5, 18, 19, 19.50, 20.0, 20.5]) * 60).astype(int)
        (np.array([18+(4/6),18+(5/6), 19, 19+(1/6), 19+(2/6), 19+(3/6)]) * 60).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        # fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.5), (d == 60 * 20)
        fan.date, fan.txt_coord, fan.cbar = date, i==0, i==4
        ds = rays.get_datasets_by_beams(
            "fhe", None, date, date + dt.timedelta(minutes=1), [18, 30]
        )
        ax = fan.overlay_simulation_fovs("fhe", ds, beams={}, gate_filter=[])
        fan.add_arc_fov("fhe", ax=ax, maxGate=25, lineColor="k", lineWidth=0.6)
        fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="m", lineWidth=0.6)
        fan.add_arc_fov(
            "fhe", ax=ax, maxGate=mgates[i], lineColor="lightgreen", lineWidth=0.6
        )
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
    fan.save(f"paper-mstid-rt/pub_figures/Figure07.png")
    fan.close()


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
    rp.save(f"paper-mstid-rt/pub_figures/Figure08.png")
    rp.close()

if 14 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [18, 30]),
        rays.RayTraceObject(
            dt.datetime(2017, 5, 27, 19),
            "fhe",
            11,
            [18, 30],
            run_name="May2017_gemini_control_cosmic",
            model_name="gemini",
        ),
    ]
    rp = rays.PlotChannels(rtos[0], nrows=3, ncols=1, xlim=[18, 24], ylim=[500, 2500])
    ax = rp.lay_rays(
        text="(A) TID", add_cbar=False, xlabel="", ylabel="Slant Range, km"
    )
    # rp.lay_rays(text="(A)", xlabel="", zoomed_in=[[500, 1200], [150, 250]], add_cbar=False)
    ax.axvline(22.7, ls="--", color="k", lw=0.5)
    ax = rp.lay_rays(
        rto=rtos[1],
        text="(B) Control",
        add_tag=False,
        xlabel="Elevation Angle, deg ($^\circ$)",
        ylabel="Slant Range, km",
    )
    ax.axvline(22.7, ls="--", color="k", lw=0.5)
    ax = rp.create_figure_pane("", "")
    ax.set_xlabel("Slant Range, km", fontdict={"size": 12, "fontweight": "bold"})
    ax.set_ylabel(
        r"Refractive Index, $\eta$", fontdict={"size": 12, "fontweight": "bold"}
    )
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [22.7, 22.7]),
        rays.RayTraceObject(
            dt.datetime(2017, 5, 27, 19),
            "fhe",
            11,
            [22.7, 22.7],
            run_name="May2017_gemini_control_cosmic",
            model_name="gemini",
        ),
    ]
    df = rtos[0].compile()
    ax.scatter(
        df.geometric_distance,
        df.refractive_index,
        marker="s",
        color="r",
        s=2,
        label="TID",
    )
    df = rtos[1].compile()
    ax.scatter(
        df.geometric_distance,
        df.refractive_index,
        marker="s",
        color="b",
        s=2,
        label="Control",
    )
    ax.set_ylim(0.85, 1)
    ax.set_xlim(0, 2500)
    ax.legend(loc=3)
    rp.fig.subplots_adjust(hspace=0.3)
    rp.save(f"paper-mstid-rt/figures/Figure14.png")
    rp.save(f"paper-mstid-rt/pub_figures/Figure09.png")
    rp.close()


if 15 in figures:
    sd = SolarDataset([dt.datetime(2017, 5, 27), dt.datetime(2017, 5, 28)])
    sd.create_stackplots("paper-mstid-rt/figures/Figure15.png")


if 17 in figures:
    date = dates[0] + dt.timedelta(minutes=int(60 * 18))
    import matplotlib.pyplot as plt

    plt.rcParams.update({"font.size": 8})
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
    fan.date, fan.txt_coord, fan.cbar = (
        dates[0] + dt.timedelta(minutes=int(60 * 20)),
        False,
        False,
    )
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

if 19 in figures:
    import matplotlib.dates as mdates

    DS = rays.get_datasets_by_beams(
        "fhe",
        [11],
        start_time=dt.datetime(2017, 5, 27, 16),
        end_time=dt.datetime(2017, 5, 27, 21),
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
    rti.save(f"paper-mstid-rt/pub_figures/Figure06.png")
    rti.close()


if 21 in figures:
    from solar import SolarDataset
    sol = SolarDataset([
        dt.datetime(2017, 5, 27),
        dt.datetime(2017, 5, 28), 
    ])
    sol._load_omni_()
    sol.create_stackplots_omni("paper-mstid-rt/figures/Figure_Analysis_Solar1.png")
    sol = SolarDataset([
        dt.datetime(2012, 12, 21),
        dt.datetime(2012, 12, 22) 
    ])
    sol._load_omni_()
    sol.create_stackplots_omni("paper-mstid-rt/figures/Figure_Analysis_Solar2.png")

if 22 in figures:
    nexrad = load_nexrad_data(
        date=dt.datetime(2017, 5, 27),
        mat_dir=Path(__file__).resolve().parent.parent / "data",
        downsample_step=4,
        start_time=dt.datetime(2017, 5, 27, 18),
        end_time=dt.datetime(2017, 5, 27, 21),
    )
    print(
        f"Loaded NEXRAD from {nexrad['file'].name}: "
        f"{len(nexrad['time'])} timestamps, "
        f"lat={nexrad['lat'].shape}, lon={nexrad['lon'].shape}"
    )


if 20 in figures:
    beam = 3
    rti = RTI(
        100,
        [dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 27, 22)],
        fov=None,
        xlim=[dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 27, 22)],
        ylim=[180, 3000],
        fig_title="",
        num_subplots=2,
    )
    frame = fds["fhw"].frame.copy()
    ax, _ = rti.addParamGsPlot(
        frame,
        beam,
        f"27 May 2017 / $f_0$= {(frame.tfreq.mean()/1e3).round(1)} MHz / Rad: fhw",
        vlim=[10, 20],
        zparam="p_l",
        xlabel="",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=True,
        times=[dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 27, 22)],
        sranges=[1400, 3000],
    )
    ax.text(
        0.05,
        0.95,
        f"(A) Beam: {beam}",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size=15),
    )
    beam = 11
    frame = fds["fhw"].frame.copy()
    ax, _ = rti.addParamGsPlot(
        frame,
        beam,
        "",
        vlim=[10, 20],
        zparam="p_l",
        label=r"$P_l$, dB",
        cmap="plasma",
        cbar=False,
        times=[dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 27, 22)],
        sranges=[1400, 3000],
    )
    ax.text(
        0.05,
        0.95,
        f"(B) Beam: {beam}",
        ha="left",
        va="center",
        transform=ax.transAxes,
        fontdict=dict(size=15),
    )
    rti.save("paper-mstid-rt/figures/Figure_Analysis_RTI.png")
    rti.close()



if 22 in figures:
    from nexrad_utils import get_nexrad_data_by_date
    from tec import _load_mat
    fan = Fan(
        rads,
        dates[0],
        nrows=3,
        ncols=2,
        add_text=True,
        # extent=[-110, -75, 28, 60],
        extent=[-103, -75, 28, 48],
        figsize=(2.5, 2)
    )
    gates = [35, 29, 34, 37, 34, 35]
    for i, d in enumerate(
        (np.array([18, 19, 19.5, 20.0, 20.5, 21.0]) * 60).astype(int)
        # (np.array([18+(4/6),18+(5/6), 19, 19+(1/6), 19+(2/6), 19+(3/6)]) * 60).astype(int)
    ):
        date = dates[0] + dt.timedelta(minutes=int(d))
        nexrad_data, tec_data = (
            get_nexrad_data_by_date(
                date,
                mat_dir="data",
                downsample_step=4,
            ),
            _load_mat(f"data/dvTEC/vTECdata_{date.strftime('%H%M')}.mat")
        )
        # fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
        fan.date, fan.txt_coord, fan.cbar = date, i==0, i==4
        ax = fan.add_axes()
        fan.add_arc_fov("fhe", ax=ax, maxGate=27, lineColor="k", lineWidth=0.3)
        fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="m", lineWidth=0.5)
        ax.text(
            0.05,
            0.95,
            f"({chr(65+i)})",
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict=dict(size="xx-small"),
        )
        im = ax.pcolormesh(
            nexrad_data["lon"],
            nexrad_data["lat"],
            nexrad_data["precip"][:, :, 0],
            shading="nearest",
            cmap="hot_r",
            vmin=0, vmax=7,
            transform=fan.geo,
            alpha=0.8,
            zorder=2,
        )
        if i == 1:
            cpos = [1.04, 0.1, 0.025, 0.8]
            cax = ax.inset_axes(cpos, transform=ax.transAxes)
            cb = fan.fig.colorbar(im, ax=ax, cax=cax)
            cb.ax.tick_params(labelsize="x-small")
            cb.set_label("Precipitation (mm/10 min)", fontsize="x-small")
        im = ax.pcolormesh(
            tec_data["lon"],
            tec_data["lat"],
            tec_data["tec"].T,
            shading="nearest",
            cmap="Blues",
            vmin=-0.1, vmax=0.1,
            transform=fan.geo,
            alpha=0.8,
            zorder=1
        )
        if i == 3:
            cpos = [1.04, 0.1, 0.025, 0.8]
            cax = ax.inset_axes(cpos, transform=ax.transAxes)
            cb = fan.fig.colorbar(im, ax=ax, cax=cax)
            cb.ax.tick_params(labelsize="x-small")
            cb.set_label("dTEC (TECU)", fontsize="x-small")
        fan.generate_fovs(
            fds, beams={"fhe": []}, discreat={"fhe": 27}, ax=ax,
            plot_discreat={"fhe": False},
        )
        # if i==3: break
    # fan.save(f"paper-mstid-rt/figures/Figure03.png")
    fan.fig.subplots_adjust(wspace=0.05, hspace=0.05)
    fan.fig.savefig(f"paper-mstid-rt/pub_figures/FigureSp02.png", bbox_inches="tight", dpi=300)
    fan.close()

if 23 in figures:
    scaled = pd.read_csv("data/scaled.csv")
    scaled["datetime"] = pd.to_datetime(
        scaled["datetime"].str.replace(r"\s+-\d+\s+", " ", regex=True)
    )
    plot_dates = [dt.datetime(2017, 5, 27, 12), dt.datetime(2017, 5, 28, 0)]
    scaled = scaled[
        (scaled["datetime"] >= plot_dates[0]) & (scaled["datetime"] <= plot_dates[1])
    ]
    sp = rays.StackPlots(2, 1, datetime=True, figsize=(7, 3), dpi=300)

    _, ax = sp.plot_stack_plots(
        np.array(scaled["datetime"]),
        scaled["foF2"],
        ylabel="Frequency (MHz)",
        color="#4C78A8",
        xlim=plot_dates,
        xlabel="",
        text="(A)",
        title="Digisonde @Boulder BC840 / 27 May 2017",
    )
    ax.lines[-1].remove()
    ax.scatter(
        scaled["datetime"], scaled["foF2"], s=18, color="#4C78A8", alpha=0.85, label="foF2"
    )
    ax.scatter(
        scaled["datetime"], scaled["foEs"], s=18, color="#F58518", alpha=0.85, label="foEs"
    )
    ax.set_ylim(0, 8)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=13, scatterpoints=3)
    ax.tick_params(axis="both", labelsize=13)
    ax.set_ylabel("Frequency (MHz)", fontsize=13)
    ax.set_title(
        "Digisonde @Boulder BC840 / 27 May 2017",
        fontsize=14,
    )
    ax.texts[-1].set_fontsize(13)

    _, ax = sp.plot_stack_plots(
        np.array(scaled["datetime"]),
        scaled["hmF2"],
        ylabel="Height (km)",
        color="#4C78A8",
        xlim=plot_dates,
        xlabel="Time, UT",
        text="(B)",
    )
    ax.lines[-1].remove()
    ax.scatter(
        scaled["datetime"], scaled["hmF2"], marker="D", s=18, color="#4C78A8", alpha=0.85, label="hmF2"
    )
    ax.scatter(
        scaled["datetime"], scaled["h`Es"], marker="D", s=18, color="#F58518", alpha=0.85, label="h'Es"
    )
    ax.set_ylim(90, 350)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper left", fontsize=13, scatterpoints=3)
    ax.tick_params(axis="both", labelsize=13)
    ax.set_ylabel("Height (km)", fontsize=13)
    ax.set_xlabel("Time, UT", fontsize=13)
    ax.texts[-1].set_fontsize(13)

    sp.save_fig("paper-mstid-rt/figures/Figure23.png")
    sp.save_fig("paper-mstid-rt/pub_figures/FigureSp03.png")
    sp.close()
