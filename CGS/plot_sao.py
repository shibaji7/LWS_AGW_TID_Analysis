import datetime as dt
from pathlib import Path

import matplotlib.dates as mdates
import pandas as pd

from pynasonde.digisonde.digi_plots import SaoSummaryPlots
from pynasonde.digisonde.digi_utils import setsize
from pynasonde.digisonde.parsers.sao import SaoExtractor

# ----------------------------
# User-facing configuration
# ----------------------------
date = dt.datetime(2017, 5, 27)
file = "/media/chakras4/Ionosonde/BC840_20170527(147)_SAO.XML"
n_procs = 12
font_size = 16
day_mode = "auto"
selected_record_index = 10
height_prange = [0, 3]

fig_dir = Path(".")
fig_dir.mkdir(parents=True, exist_ok=True)

setsize(font_size)


def _set_day_axis(ax, base_date):
    """Constrain x-axis to one day with 6-hour ticks."""
    # ax.set_xlim([base_date, base_date + dt.timedelta(1)])
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))


def _dedupe_legend(plot_obj):
    """Combine and de-duplicate legends across all figure axes."""
    main_ax = plot_obj.axes
    handles, labels = [], []
    for axis in plot_obj.fig.get_axes():
        h, l = axis.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)
    seen = set()
    u_handles, u_labels = [], []
    for h, l in zip(handles, labels):
        if l and l not in seen:
            seen.add(l)
            u_handles.append(h)
            u_labels.append(l)
    if u_handles:
        main_ax.legend(u_handles, u_labels, loc=2)


# ---------------------------------
# 1) Full-day height-profile panel
# ---------------------------------
sao = SaoExtractor(
    file,
    extract_time_from_name=True,
    extract_stn_from_name=True
)
sao.extract_xml()
df = sao.get_scaled_datasets_xml()


sao_plot = SaoSummaryPlots(
    figsize=(8, 4),
    fig_title="BC840 / F2 (scaled, day-file multi-record), 27 May 2017",
    draw_local_time=False,
    font_size=font_size,
)
df.datetime = [x.split(" ")[0] + " " + x.split(" ")[2] for x in df.datetime]
print(df.datetime[0])
df.datetime = pd.to_datetime(df.datetime)
sao_plot.plot_TS(
    df,
    right_yparams=["hmF2"],
    left_yparams=["foF2"],
    right_ylim=[100, 450],
    left_ylim=[1, 15],
)
_set_day_axis(sao_plot.axes, date)
_dedupe_legend(sao_plot)
sao_plot.save(str(fig_dir / "stack_sao_multi_F2.png"))