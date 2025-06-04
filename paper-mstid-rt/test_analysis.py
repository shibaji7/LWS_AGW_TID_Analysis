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
    rp.save(f"paper-mstid-rt/pub_figures/Figure10.png")
    rp.close()
