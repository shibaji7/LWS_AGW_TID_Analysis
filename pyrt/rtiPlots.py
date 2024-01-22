#!/usr/bin/env python

"""
    rtiPlots.py: module to plot RTI plots with various y-axis transformation
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

import mplstyle
#plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import datetime as dt

import numpy as np
import matplotlib.dates as mdates
import pandas as pd
import utils

class RTIPlots(object):
    """
    Create plots for velocity, width, power, elevation angle, etc.
    """

    def __init__(
        self,
        nGates,
        drange,
        fig_title=None,
        num_subplots=1,
        ylim=[180, 3000],
    ):
        self.nGates = nGates
        self.drange = drange
        self.num_subplots = num_subplots
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(6, 3 * num_subplots), dpi=240)
        if fig_title:
            plt.suptitle(
                fig_title, x=0.075, y=0.99, ha="left", fontweight="bold", fontsize=15
            )
        self.ylim = ylim
        return

    def addParamPlot(
        self,
        beam_soundings_rays,
        xlabel="Time (UT)",
        cmap="jet",
        vlim=[-10, 0],
        alpha=1,
    ):
        o = pd.DataFrame()
        for sound in beam_soundings_rays:
            o = pd.concat([o, sound.ray_power])
        o.lag_power = 10*np.log10(o.lag_power/o.lag_power.max())
        ax = self._add_axis()
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H^{%M}"))
        hours = mdates.HourLocator(byhour=range(0, 24, 1))
        ax.xaxis.set_major_locator(hours)
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim(
                [mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])]
        )
        ax.set_ylim(self.ylim)
        X, Y, Z = utils.get_gridded_parameters(
            o, xparam="date", yparam="srange", zparam="lag_power", rounding=False
        )
        im = ax.pcolormesh(
            X,
            Y,
            #(Y*45)+180,
            Z.T,
            lw=0.01,
            edgecolors="None",
            cmap=cmap,
            snap=True,
            vmax=vlim[1],
            vmin=vlim[0],
            shading="auto",
            alpha=alpha,
        )
        self._add_colorbar(self.fig, ax, im, label=r"$\lambda$ Power (dB)")
        ax.set_ylabel("Slant Range (km)", fontdict={"size": 12, "fontweight": "bold"})
        return ax

    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return ax

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
    
    def _add_colorbar(
        self,
        fig,
        ax,
        im,
        label="",
        xOff=0,
        yOff=0,
    ):
        """
        Add a colorbar to the right of an axis.
        """
        cpos = [1.04 + xOff, 0.1 + yOff, 0.025, 0.8]
        cax = ax.inset_axes(cpos, transform=ax.transAxes)
        cb = fig.colorbar(im, ax=ax, cax=cax)
        cb.set_label(label)
        return
    
def create_RTI(folder, beam_soundings_rays):
    rti = RTIPlots(
        80, 
        [beam_soundings_rays[0].date, beam_soundings_rays[-1].date],
        f"GEMINI2D/{beam_soundings_rays[0].rad}/{beam_soundings_rays[0].beam}"
    )
    rti.addParamPlot(beam_soundings_rays)
    rti.save(folder + "RTI.png")
    rti.close()
    return