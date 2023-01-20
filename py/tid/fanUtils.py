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

import matplotlib
import matplotlib.pyplot as plt
import numpy as np

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import cartopy
import matplotlib.ticker as mticker
from cartopy.mpl.gridliner import LATITUDE_FORMATTER, LONGITUDE_FORMATTER
from cartoUtils import SDCarto


class Fan(object):
    """
    This class holds plots for all radars FoVs
    """

    def __init__(self, rads, date, fig_title=None, nrows=1, ncols=1, coord="geo"):
        self.rads = rads
        self.date = date
        self.nrows, self.ncols = nrows, ncols
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(3 * ncols, 3 * nrows), dpi=180)
        self.coord = coord
        plt.suptitle(
            f"{self.date_string()} / {fig_title}"
            if fig_title
            else f"{self.date_string()}",
            x=0.1,
            y=0.82,
            ha="left",
            fontweight="bold",
            fontsize=8,
        )
        return

    def add_axes(self):
        """
        Instatitate figure and axes labels
        """
        self._num_subplots_created += 1
        proj = cartopy.crs.Orthographic(central_longitude=-90.0, central_latitude=45.0)
        ax = self.fig.add_subplot(
            100 * self.nrows + 10 * self.ncols + self._num_subplots_created,
            projection="SDCarto",
            map_projection=proj,
            coords=self.coord,
            plot_date=self.date,
        )
        ax.overaly_coast_lakes(lw=0.4, alpha=0.4)
        ax.set_extent([-130, -50, 30, 70], crs=cartopy.crs.PlateCarree())
        plt_lons = np.arange(-180, 181, 15)
        mark_lons = np.arange(-180, 181, 30)
        plt_lats = np.arange(40, 90, 10)
        gl = ax.gridlines(crs=cartopy.crs.PlateCarree(), linewidth=0.5)
        gl.xlocator = mticker.FixedLocator(plt_lons)
        gl.ylocator = mticker.FixedLocator(plt_lats)
        gl.xformatter = LONGITUDE_FORMATTER
        gl.yformatter = LATITUDE_FORMATTER
        gl.n_steps = 90
        ax.mark_latitudes(plt_lats, fontsize="small", color="darkred")
        ax.mark_longitudes(mark_lons, fontsize="small", color="darkblue")
        self.proj = proj
        self.geo = cartopy.crs.PlateCarree()
        ax.text(
            -0.02,
            0.99,
            "Coord: Geo",
            ha="center",
            va="top",
            transform=ax.transAxes,
            fontsize="small",
            rotation=90,
        )
        return ax

    def date_string(self, label_style="web"):
        # Set the date and time formats
        dfmt = "%d/%b/%Y" if label_style == "web" else "%d %b %Y,"
        tfmt = "%H:%M"
        stime = self.date
        date_str = "{:{dd} {tt}} UT".format(stime, dd=dfmt, tt=tfmt)
        return date_str

    def generate_fov(self, fds, beams=[]):
        """
        Generate plot with dataset overlaid
        """
        ax = self.add_axes()
        for rad in self.rads:
            ax.overlay_radar(rad)
            ax.overlay_fov(rad)
            ax.overlay_data(rad, fds[rad].frame, self.proj)
            if beams and len(beams) > 0:
                [ax.overlay_fov(rad, beamLimits=[b, b + 1], ls="--") for b in beams]
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
