import datetime as dt
import os
import sys

sys.path.extend(["py/", "py/txUtils/", "py/tid/", "py/davitPy/", "paper-mstid-rt/"])
import rays
from fetchUtils import FetchData
from fanUtils import Fan
import numpy as np
import pandas as pd

figures = [1, 3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
# figures = [3]
rads = ["fhe"]
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

# if 13 in figures:
#     rtos = [
#         rays.RayTraceObject(dt.datetime(2017, 5, 27, 19, 30), "fhe", 11, [21.6, 21.8]),
#         rays.RayTraceObject(
#             dt.datetime(2017, 5, 27, 19, 30),
#             "fhe",
#             11,
#             [21.6, 21.8],
#             run_name="May2017_gemini_control_cosmic",
#             model_name="gemini",
#         ),
#     ]
#     rp = rays.PlotRays(
#         rtos[0], nrows=2, ncols=1, lw=0.5, arc=True, figsize=(8, 4), ylim=[-300, 600]
#     )
#     rp.lay_rays(
#         text="(A) TID",
#         zoomed_in=[[500, 1200], [150, 250]],
#         add_cbar=False,
#         ped_angles=[21.7],
#     )
#     rp.lay_rays(
#         rto=rtos[1],
#         text="(B) Control",
#         ylabel="",
#         # zoomed_in=[[500, 1200], [150, 250]],
#         add_tag=False,
#     )
#     # rp.fig.subplots_adjust(hspace=1.2)
#     rp.save(f"paper-mstid-rt/figures/Figure13.png")
#     rp.save(f"paper-mstid-rt/pub_figures/Figure09.png")
#     rp.close()


# if 11 in figures:
#     rtos = [
#         rays.RayTraceObject(dt.datetime(2017, 5, 27, 17, 30), "fhe", 11, [18, 30]),
#         rays.RayTraceObject(dt.datetime(2017, 5, 27, 18, 15), "fhe", 11, [18, 30]),
#         rays.RayTraceObject(dt.datetime(2017, 5, 27, 19), "fhe", 11, [18, 30]),
#         rays.RayTraceObject(dt.datetime(2017, 5, 27, 19, 30), "fhe", 11, [18, 30]),
#     ]
#     rp = rays.PlotRays(rtos[0], nrows=4, ncols=1, arc=True)
#     rp.lay_rays(
#         xlabel="",
#         add_cbar=False,
#         text="(A)",
#         tag_distance=1400,
#         ped_angles=[
#             24.3,
#         ],
#     )
#     rp.lay_rays(
#         rto=rtos[1],
#         add_cbar=False,
#         add_tag=False,
#         text="(B)",
#         tag_distance=1400,
#         ped_angles=[22.6],
#     )
#     rp.lay_rays(
#         rto=rtos[2],
#         add_cbar=False,
#         add_tag=False,
#         text="(C)",
#         tag_distance=1400,
#         ped_angles=[22.7, 23, 21.3, 21.7],
#     )
#     rp.lay_rays(
#         rto=rtos[3],
#         add_tag=False,
#         text="(D)",
#         tag_distance=1400,
#         ped_angles=[21.7, 23.4],
#     )
#     rp.save(f"paper-mstid-rt/figures/Figure11.png")
#     rp.save(f"paper-mstid-rt/pub_figures/Figure06.png")
#     rp.close()


if 11 in figures:
    rtos = [
        rays.RayTraceObject(dt.datetime(2017, 5, 27, 16) + dt.timedelta(minutes=i), "fhe", 11, [18, 30])
        for i in range(300)
    ]
    for i, rt in enumerate(rtos):
        rp = rays.PlotRays(rt, nrows=1, ncols=1, arc=True)
        rp.lay_rays(
            text="",
            tag_distance=1400,
        )
        rp.save("figures/movies/Figure%04d.png"%i)
        rp.close()


# if 12 in figures:
#     fan = Fan(
#         rads,
#         dates[0],
#         nrows=3,
#         ncols=2,
#         add_text=True,
#     )
#     # gates = [38, 35, 29, 33, 37, 35]
#     gates = [32, 32, 29, 33, 34, 35]
#     mgates = [31, 31, 28, 30, 30, 29]
#     for i, d in enumerate(
#         # (np.array([17.5, 18, 19, 19.50, 20.0, 20.5]) * 60).astype(int)
#         (np.array([18+(4/6),18+(5/6), 19, 19+(1/6), 19+(2/6), 19+(3/6)]) * 60).astype(int)
#     ):
#         date = dates[0] + dt.timedelta(minutes=int(d))
#         # fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.5), (d == 60 * 20)
#         fan.date, fan.txt_coord, fan.cbar = date, i==0, i==4
#         ds = rays.get_datasets_by_beams(
#             "fhe", None, date, date + dt.timedelta(minutes=1), [18, 30]
#         )
#         ax = fan.overlay_simulation_fovs("fhe", ds, beams={}, gate_filter=[])
#         fan.add_arc_fov("fhe", ax=ax, maxGate=25, lineColor="k", lineWidth=0.6)
#         fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="m", lineWidth=0.6)
#         fan.add_arc_fov(
#             "fhe", ax=ax, maxGate=mgates[i], lineColor="lightgreen", lineWidth=0.6
#         )
#         ax.text(
#             0.05,
#             0.95,
#             f"({chr(65+i)})",
#             ha="left",
#             va="center",
#             transform=ax.transAxes,
#             fontdict=dict(size="xx-small"),
#         )
#     fan.save(f"paper-mstid-rt/figures/Figure12.png")
#     fan.save(f"paper-mstid-rt/pub_figures/Figure07.png")
#     fan.close()

# if 4 in figures:
#     fan = Fan(
#         rads,
#         dates[0],
#         nrows=3,
#         ncols=2,
#         add_text=True,
#     )
#     gates = [32, 32, 29, 33, 34, 35]
#     for i, d in enumerate(
#         # (np.array([17.5, 18, 19, 19.5, 20.0, 20.5]) * 60).astype(int)
#         (np.array([18+(4/6),18+(5/6), 19, 19+(1/6), 19+(2/6), 19+(3/6)]) * 60).astype(int)
#     ):
#         date = dates[0] + dt.timedelta(minutes=int(d))
#         # fan.date, fan.txt_coord, fan.cbar = date, (d == 60 * 17.25), (d == 60 * 20)
#         fan.date, fan.txt_coord, fan.cbar = date, i==0, i==4
#         ax = fan.generate_fovs(
#             fds, beams={"fhe": [11, 3]}, discreat={"fhe": 27}
#         )
#         fan.add_arc_fov("fhe", ax=ax, maxGate=27, lineColor="k", lineWidth=0.6)
#         fan.add_arc_fov("fhe", ax=ax, maxGate=gates[i], lineColor="m", lineWidth=1.2)
#         ax.text(
#             0.05,
#             0.95,
#             f"({chr(65+i)})",
#             ha="left",
#             va="center",
#             transform=ax.transAxes,
#             fontdict=dict(size="xx-small"),
#         )
#     fan.save(f"paper-mstid-rt/figures/Figure03.png")
#     fan.save(f"paper-mstid-rt/pub_figures/Figure03.png")
#     fan.close()