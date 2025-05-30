#!/usr/bin/env python

"""
fanUtils.py: module to plot Fan plots with various transformation
"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import datetime as dt
import os

import matplotlib.pyplot as plt
import numpy as np
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
plt.rcParams["text.usetex"] = False
import glob

import cartopy
import cv2
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import tidUtils
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER


def get_color_by_number(number, cmap_name="viridis"):
    """This function takes a number and a colormap name, and returns a color."""
    cmap = cm.get_cmap(cmap_name)
    norm = plt.Normalize(0, 26)  # Normalize the number to the colormap range
    return cmap(norm(number))


class Fan(object):
    """
    This class holds plots for all radars FoVs
    """

    def __init__(
        self,
        rads,
        date,
        fig_title=None,
        nrows=1,
        ncols=1,
        coord="geo",
        cs=False,
        tec=None,
        tec_times=None,
        txt_coord=False,
        cbar=False,
        add_text=False,
        extent=[-110, -75, 35, 60],
    ):
        # if cs:
        #     plt.style.use(["science", "ieee"])
        self.cs = cs
        self.rads = rads
        self.date = date
        self.nrows, self.ncols = nrows, ncols
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(2 * ncols, 2 * nrows), dpi=1000)
        self.coord = coord
        self.tec, self.tec_times = tec, tec_times
        self.fig_title = fig_title
        self.txt_coord = txt_coord
        self.cbar = cbar
        self.add_text = add_text
        self.extent = extent
        return

    def add_axes(self):
        """
        Instatitate figure and axes labels
        """
        self._num_subplots_created += 1
        proj = cartopy.crs.Stereographic(central_longitude=-90.0, central_latitude=45.0)
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="SDCarto",
            map_projection=proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        ax.set_extent(self.extent, crs=cartopy.crs.PlateCarree())
        plt_lons = np.linspace(self.extent[0], self.extent[1], 3)
        mark_lons = np.linspace(self.extent[0], self.extent[1], 3)
        plt_lats = np.arange(40, 90, 10)
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        ax.mark_latitudes(plt_lats, fontsize="xx-small", color="darkred")
        ax.mark_longitudes(mark_lons, fontsize="xx-small", color="darkblue")
        self.proj = proj
        self.geo = cartopy.crs.PlateCarree()
        if self.txt_coord:
            ax.text(
                -0.02,
                0.99,
                "Coord: Geo",
                ha="center",
                va="top",
                transform=ax.transAxes,
                fontsize="x-small",
                rotation=90,
            )
        if self._num_subplots_created == 1 or self.add_text:
            ax.text(
                0.05,
                1.05,
                (
                    f"{self.date_string()} / {self.fig_title}"
                    if self.fig_title
                    else f"{self.date_string()}"
                ),
                ha="left",
                va="center",
                fontweight="bold",
                fontsize="xx-small",
                transform=ax.transAxes,
            )
        return ax

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d %b %Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime = self.date
        date_str = "{:{dd} {tt}} UT".format(stime, dd=dfmt, tt=tfmt)
        return date_str

    def generate_fov(
        self,
        rad,
        frame,
        beams=[],
        discreat=0,
        ax=None,
        laytec=False,
        maxGate=45,
        col="k",
        param="p_l",
        label=r"$P_l$, dB",
        p_max=18,
        p_min=12,
    ):
        """
        Generate plot with dataset overlaid
        """
        ax = ax if ax else self.add_axes()
        if laytec:
            ipplat, ipplon, dtec = tidUtils.fetch_tec_by_datetime(
                self.date, self.tec, self.tec_times
            )
            ax.overlay_tec(ipplat, ipplon, dtec, self.proj)
        ax.overlay_radar(rad, font_color=col)
        ax.overlay_fov(rad, lineColor=col)
        if len(frame) > 0:
            ax.overlay_data(
                rad,
                frame,
                self.proj,
                maxGate=maxGate,
                cbar=self.cbar,
                label=label,
                plot_discreat=discreat,
                p_name=param,
                p_max=p_max,
                p_min=p_min,
            )
        if beams and len(beams) > 0:
            [
                ax.overlay_fov(
                    rad,
                    beamLimits=[b, b + 1],
                    ls="-",
                    lineColor=get_color_by_number(b),
                    lineWidth=0.4,
                )
                for b in beams
            ]
        return

    def add_arc_fov(
        self,
        rad,
        ax=None,
        lineColor="m",
        maxGate=30,
        beamLimits=None,
        lineWidth=0.4,
        ls="-",
    ):
        """
        Generate plot with dataset overlaid
        """
        ax = ax if ax else self.add_axes()
        ax.add_arc_fov(
            rad,
            lineColor=lineColor,
            maxGate=maxGate,
            beamLimits=beamLimits,
            lineWidth=lineWidth,
            ls=ls,
        )
        return ax

    def generate_fovs(
        self,
        fds,
        beams={},
        discreat={},
        laytec=False,
        param="p_l",
        label=r"$P_l$, dB",
        p_max=18,
        p_min=12,
    ):
        """
        Generate plot with dataset overlaid
        """
        ax = self.add_axes()
        for rad in self.rads:
            self.generate_fov(
                rad,
                fds[rad].frame,
                beams[rad],
                discreat[rad],
                ax,
                laytec,
                col=fds[rad].color,
                param=param,
                p_max=p_max,
                p_min=p_min,
                label=label,
            )
        return ax

    def overlay_simulation_fovs(
        self, rad, df, beams=[], gate_filter=[], ax=None, col="b"
    ):
        """
        Generate plot with dataset overlaid
        """
        ax = ax if ax else self.add_axes()
        ax.overlay_radar(rad, font_color=col)
        ax.overlay_fov(rad, lineColor=col)
        if len(df) > 0:
            if len(gate_filter) == 2:
                df = df[(df.slist >= gate_filter[0]) & (df.slist <= gate_filter[1])]
            ax.overlay_data(
                rad,
                df,
                self.proj,
                maxGate=45,
                cbar=self.cbar,
                p_name="p_l",
                scan_time=1,
                p_max=-85,
                p_min=-90,
                label=r"$P_r$, dB",
            )
        if beams and len(beams) > 0:
            [
                ax.overlay_fov(
                    rad,
                    beamLimits=[b, b + 1],
                    ls="-",
                    lineColor=get_color_by_number(b),
                    lineWidth=0.4,
                )
                for b in beams
            ]
        return ax

    def overlay_fovs(self, rad, beams=[], ax=None, col="k"):
        """
        Generate plot with dataset overlaid
        """
        ax = ax if ax else self.add_axes()
        ax.overlay_radar(rad, font_color=col)
        ax.overlay_fov(rad, lineColor=col)
        if beams and len(beams) > 0:
            [
                ax.overlay_fov(
                    rad,
                    beamLimits=[b, b + 1],
                    ls="-",
                    lineColor=get_color_by_number(b),
                    lineWidth=0.3,
                )
                for b in beams
            ]
        return ax

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return


def create_movie(folder, outfile, pat, fps=3):
    """
    Create movies from pngs
    """
    files = glob.glob(f"{folder}/{pat}")
    files.sort()
    fourcc = cv2.VideoWriter_fourcc(*"XVID")
    img = cv2.imread(files[0])
    height, width, layers = img.shape
    size = (width, height)
    out = cv2.VideoWriter(f"{folder}/{outfile}", fourcc, fps, size)
    for idx in range(len(files)):
        img = cv2.imread(files[idx])
        out.write(img)
    out.release()
    return


def create_ovearlay_movies(
    fds, date, rads, ovearlay_tec="data/2022-12-21/WS355.mat", fps=15, clear=False
):
    """
    Create Fov-Fan plots and ovearlay movies
    """

    def plot_fan(d, tec, tec_times, file):
        fov = Fan(rads, d, tec=tec, tec_times=tec_times)
        fov.generate_fovs(fds, laytec=ovearlay_tec is not None)
        fov.save(file)
        fov.close()
        return

    folder = tidUtils.get_folder(date)
    if ovearlay_tec:
        tec, tec_times = tidUtils.read_tec_mat_files(ovearlay_tec)
    time_range = [
        min([fds[rad].scans[0].stime for rad in rads]),
        max([fds[rad].scans[-1].etime for rad in rads]),
    ]
    times = [
        time_range[0] + dt.timedelta(seconds=i * 60)
        for i in range(int((time_range[1] - time_range[0]).total_seconds() / 60))
    ]
    for d in times:
        ipplat, ipplon, dtec = tidUtils.fetch_tec_by_datetime(
            d,
            tec,
            tec_times,
        )
        file = tidUtils.get_folder(d) + f"/Fan,{d.strftime('%H-%M')}.png"
        if clear:
            plot_fan(d, tec, tec_times, file)
        elif not os.path.exists(file):
            plot_fan(d, tec, tec_times, file)
    create_movie(tidUtils.get_folder(date), fps)
    return
