#!/usr/bin/env python

"""
    plots.py: module to plot various type of dataset
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

plt.style.use(["science", "ieee"])
import matplotlib.dates as mdates
import numpy as np
import utils
import swifter


class RTI(object):
    """
    Create plots for velocity, width, power, elevation angle, etc.
    """

    def __init__(self, nrang, drange, fig_title="", num_subplots=1):
        self.nrang = nrang
        self.unique_gates = np.linspace(1, nrang, nrang)
        self.num_subplots = num_subplots
        self.drange = drange
        self._num_subplots_created = 0
        self.fig = plt.figure(
            figsize=(6, 3 * num_subplots), dpi=180
        )  # Size for website
        plt.suptitle(
            fig_title, x=0.075, y=0.99, ha="left", fontweight="bold", fontsize=15
        )
        matplotlib.rcParams.update({"font.size": 10})
        return

    def addParamPlot(
        self,
        df,
        beam,
        title,
        p_max=36,
        p_min=0,
        p_step=3,
        xlabel="Time UT",
        zparam="p_l",
        label="Power [dB]",
    ):
        ax = self._add_axis()
        df = df[df.bmnum == beam]
        df["mdates"] = df.time.swifter.apply(lambda x: mdates.date2num(x))
        X, Y, Z = utils.get_gridded_parameters(
            df, xparam="mdates", yparam="slist", zparam=zparam, rounding=False
        )
        bounds = list(range(p_min, p_max + 1, p_step))
        cmap = plt.cm.jet
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        # cmap.set_bad("w", alpha=0.0)
        # Configure axes
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        hours = mdates.HourLocator(byhour=range(0, 24, 4))
        ax.xaxis.set_major_locator(hours)
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])])
        ax.set_ylim([0, self.nrang])
        ax.set_ylabel("Range gate", fontdict={"size": 12, "fontweight": "bold"})
        ax.pcolormesh(
            X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, norm=norm, snap=True
        )
        self._add_colorbar(self.fig, ax, bounds, cmap, label=label)
        ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        return

    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return ax

    def _add_colorbar(self, fig, ax, bounds, colormap, label=""):
        """
        Add a colorbar to the right of an axis.
        """
        import matplotlib as mpl

        pos = ax.get_position()
        cpos = [
            pos.x1 + 0.025,
            pos.y0 + 0.0125,
            0.015,
            pos.height * 0.9,
        ]  # this list defines (left, bottom, width, height
        cax = fig.add_axes(cpos)
        norm = mpl.colors.BoundaryNorm(bounds, colormap.N)
        cb2 = mpl.colorbar.ColorbarBase(
            cax,
            cmap=colormap,
            norm=norm,
            ticks=bounds,
            spacing="uniform",
            orientation="vertical",
        )
        cb2.set_label(label)
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return


class IntervalPlots(object):
    """
    Create Interval Plots
    """

    def __init__(self, o, fq=[], wv=[]):
        """
        Create Series plots
        """
        self.o = o
        self.fq = fq
        self.wv = wv
        self._series_()
        return

    def _series_(self):
        """
        Series plots on frequency and wavenumber
        """
        self.num_subplots = len(self.fq) + len(self.wv)
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(6, 3 * self.num_subplots), dpi=180)
        for f in self.fq:
            ax = self._add_axis()
            x = self.o["freq"][f[0]][f[1]]
            ax.semilogx(x["f"], x["pow"], "k-", lw=0.8)
            ax.set_xlim(1e-6, 1e0)
            ax.set_ylim(0, 1)
            ax.set_xlabel(r"$f_0$, Hz")
            ax.set_ylabel("Power")
        for w in self.wv:
            ax = self._add_axis()
            x = self.o["wv"][f[0]][f[1]]
            ax.semilogx(x["wvn"], x["pow"], "k-", lw=0.8)
            ax.set_xlim(1e-10, 1e-4)
            ax.set_ylim(0, 1)
            ax.set_xlabel(r"$\lambda^{-1}$, $rad/m$")
            ax.set_ylabel("Power")
        self.fig.subplots_adjust(hspace=0.5, wspace=0.5)
        return

    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return ax

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return
