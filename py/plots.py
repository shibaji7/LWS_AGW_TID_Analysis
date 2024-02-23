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

#plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
import matplotlib.dates as mdates
import numpy as np
import utils


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
        xlabel="Time [UT]",
        zparam="p_l",
        label="Power [dB]",
        vlines=[],
        hlines=[],
    ):
        ax = self._add_axis()
        df = df[df.bmnum == beam]
        X, Y, Z = utils.get_gridded_parameters(
            df, xparam="time", yparam="srange", zparam=zparam, rounding=False
        )
        bounds = list(range(p_min, p_max + 1, p_step))
        cmap = plt.cm.jet
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        # cmap.set_bad("w", alpha=0.0)
        # Configure axes
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
        hours = mdates.HourLocator(byhour=range(0, 24, 6))
        ax.xaxis.set_major_locator(hours)
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])])
        ax.set_ylim([0, 4000])
        ax.set_ylabel("Slant Range [km]", fontdict={"size": 12, "fontweight": "bold"})
        ax.pcolormesh(
            X, Y, Z.T, lw=0.01, edgecolors="None", cmap=cmap, norm=norm, snap=True
        )
        X, Y, Z = utils.get_gridded_parameters(
            df, xparam="time", yparam="glat", zparam=zparam, rounding=False
        )
        ax = ax.twinx()
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
        ax.set_ylabel("Latitude [deg]", fontdict={"size": 12, "fontweight": "bold"})
        ax.plot([X[0, 0]] * len(Y[:, 0]), Y[:, 0], "w", ls="None")
        self._add_colorbar(self.fig, ax, bounds, cmap, label=label)
        ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        for vline in vlines:
            ax.axvline(vline, color="green", ls="-", lw=0.7, alpha=0.7)
            ax.axvline(vline, color="k", ls="--", lw=0.6)
        for hline in hlines:
            ax.axhline(hline, color="green", ls="-", lw=0.7, alpha=0.7)
            ax.axhline(hline, color="k", ls="--", lw=0.6)
        return

    def add_series(
        self,
        df,
        beam,
        gate,
        p_lim=[-10, 10],
        xlabel="Time UT",
        ylabel="Power [dB]",
        zparam="p_l",
    ):
        ax = self._add_axis()
        df = df[(df.bmnum == beam) & (df.slist == gate)]
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("%H:%M"))
        hours = mdates.HourLocator(byhour=range(0, 24, 1))
        ax.xaxis.set_major_locator(hours)
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_ylabel(ylabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim([mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])])
        ax.plot(df.time, df[zparam], "r--", lw=0.8)
        ax.set_ylim(p_lim)
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
            pos.x1 + 0.09,
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
            ax.set_xlabel(r"$\lambda^{-1}$, $/m$")
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


def histogram_plots(fqs, filepath):
    """
    Plot histograms
    """
    fig = plt.figure(figsize=(6, 3), dpi=180)
    ax = fig.add_subplot(121)
    ax.hist(np.array(fqs) * 1e3, bins=20, histtype="step", color="r")
    ax.set_yscale("log")
    ax.set_ylabel("Counts")
    ax.set_xlabel(r"$f_0$, mHz")
    #     ax = fig.add_subplot(122)
    #     ax.hist(np.array(wvs) * 1e3, bins=20, histtype="step", color="r")
    #     ax.set_yscale("log")
    #     ax.set_ylabel("Counts")
    #     ax.set_xlabel(r"$\lambda^{-1}$, /km")
    fig.subplots_adjust(hspace=0.5, wspace=0.5)
    fig.savefig(filepath, bbox_inches="tight")
    return


class PlotKarr(object):
    """
    Plot Dlm and Karr
    """

    def __init__(
        self,
        Karr,
        kxVec,
        kyVec,
        maxSignals=None,
        sig_fontsize=24,
        plot_title=True,
        cbar_ticks=None,
        cbar_shrink=1.0,
        cbar_fraction=0.15,
        cbar_gstext_offset=-0.075,
        cbar_gstext_fontsize=None,
        **kwArgs
    ):
        self.kxVec = kxVec
        self.kyVec = kyVec
        self.Karr = Karr
        self.fig = plt.figure(figsize=(5, 5))
        self.axis = self.fig.add_subplot(111, aspect="equal")
        self.draw()
        return

    def draw(self):
        sig_fontsize = 24
        x_labelpad = None
        y_labelpad = None
        cbar_ticks = None
        cbar_shrink = 1.0
        cbar_fraction = 0.15
        cbar_gstext_offset = -0.075
        cbar_gstext_fontsize = None
        cbar_pad = 0.05
        cmap = None
        plot_colorbar = True

        import matplotlib.patheffects as PathEffects
        from matplotlib.collections import PolyCollection
        from scipy import stats

        data = np.abs(self.Karr) - np.min(np.abs(self.Karr))
        sd = np.nanstd(data, axis=None)
        mean = np.nanmean(data, axis=None)
        scMax = mean + 6.5 * sd
        data = data / scMax
        scale = [0.0, 1.0]
        nrL, nrM = np.shape(data)
        verts = []
        scan = []
        # Plot Spectrum
        for ll in range(nrL - 1):
            xx0 = self.kxVec[ll]
            xx1 = self.kxVec[ll + 1]
            for mm in range(nrM - 1):
                scan.append(data[ll, mm])

                yy0 = self.kyVec[mm]
                yy1 = self.kyVec[mm + 1]

                x1, y1 = xx0, yy0
                x2, y2 = xx1, yy0
                x3, y3 = xx1, yy1
                x4, y4 = xx0, yy1
                verts.append(((x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)))
        cmap = matplotlib.cm.jet
        bounds = np.linspace(scale[0], scale[1], 256)
        norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)
        pcoll = PolyCollection(
            np.array(verts),
            edgecolors="face",
            linewidths=0,
            closed=False,
            cmap=cmap,
            norm=norm,
            zorder=99,
        )
        pcoll.set_array(np.array(scan))
        self.axis.add_collection(pcoll, autolim=False)
        self.axis.axvline(color="0.82", lw=2, zorder=150)
        self.axis.axhline(color="0.82", lw=2, zorder=150)
        cbar_label = "Normalized Wavenumber Power"
        if plot_colorbar:
            cbar = self.fig.colorbar(
                pcoll,
                orientation="vertical",
                shrink=cbar_shrink,
                fraction=cbar_fraction,
                pad=cbar_pad,
            )
            cbar.set_label(cbar_label)
            if not cbar_ticks:
                cbar_ticks = np.arange(10) / 10.0
            cbar.set_ticks(cbar_ticks)

        #             if currentData.metadata.has_key('gscat'):
        #                 if currentData.metadata['gscat'] == 1:
        #                     cbar.ax.text(0.5,cbar_gstext_offset,'Ground\nscat\nonly',ha='center',fontsize=cbar_gstext_fontsize)
        self.axis.set_xlim([np.min(self.kxVec), np.max(self.kxVec)])
        self.axis.set_ylim([np.min(self.kyVec), np.max(self.kyVec)])

        ticks = self.axis.get_xticks()
        newLabels = []
        for x in range(len(ticks)):
            tck = ticks[x]
            if tck != 0:
                km = 2 * np.pi / tck
                km_txt = "%i" % km
            else:
                km_txt = ""

            rad_txt = "%.2f" % tck
            txt = "\n".join([rad_txt, km_txt])
            newLabels.append(txt)

        self.axis.set_xticklabels(newLabels)
        self.axis.set_xlabel(
            "kx [rad]\n$\lambda$ [km]", ha="center", labelpad=x_labelpad
        )

        ticks = self.axis.get_yticks()
        newLabels = []
        for y in range(len(ticks)):
            tck = ticks[y]
            if tck != 0:
                km = 2 * np.pi / tck
                km_txt = "%i" % km
            else:
                km_txt = ""

            rad_txt = "%.2f" % tck
            txt = "\n".join([km_txt, rad_txt])
            newLabels.append(txt)
        self.axis.set_yticklabels(newLabels, rotation=90.0)
        self.axis.set_ylabel(
            "ky [rad]\n$\lambda$ [km]", va="center", labelpad=y_labelpad
        )
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight")
        plt.close()
        return


def plot_Dlm(Dlm):
    from matplotlib import pyplot as plt
    from matplotlib.collections import PolyCollection

    fig = plt.figure(figsize=(4, 4), dpi=180)

    import copy

    from scipy import stats

    data = np.abs(Dlm)

    # Determine scale for colorbar.
    sd = np.nanstd(data, axis=None)
    mean = np.nanmean(data, axis=None)
    scMax = mean + 4.0 * sd
    scale = scMax * np.array([0, 1.0])

    # Do plotting here!
    axis = fig.add_subplot(111)

    nrL, nrM = np.shape(data)

    verts = []
    scan = []
    # Plot Spectrum
    for ll in range(nrL):
        xx0 = ll
        xx1 = ll + 1
        for mm in range(nrM):
            scan.append(data[ll, mm])

            yy0 = mm
            yy1 = mm + 1

            x1, y1 = xx0, yy0
            x2, y2 = xx1, yy0
            x3, y3 = xx1, yy1
            x4, y4 = xx0, yy1
            verts.append(((x1, y1), (x2, y2), (x3, y3), (x4, y4), (x1, y1)))

    colors = "lasse"
    if scale is None:
        scale = (np.min(scan), np.max(scan))
    cmap = matplotlib.cm.jet
    bounds = np.linspace(scale[0], scale[1], 256)
    norm = matplotlib.colors.BoundaryNorm(bounds, cmap.N)

    pcoll = PolyCollection(
        np.array(verts),
        edgecolors="face",
        linewidths=0,
        closed=False,
        cmap=cmap,
        norm=norm,
        zorder=99,
    )
    pcoll.set_array(np.array(scan))
    axis.add_collection(pcoll, autolim=False)

    # Colorbar
    cbar = fig.colorbar(pcoll, orientation="vertical")  # ,shrink=.65,fraction=.1)
    cbar.set_label("ABS(Spectral Density)")
    #     if currentData.metadata.has_key('gscat'):
    #         if currentData.metadata['gscat'] == 1:
    #             cbar.ax.text(0.5,-0.075,'Ground\nscat\nonly',ha='center')
    #  labels[-1].set_visible(False)
    axis.set_xlim([0, nrL])
    axis.set_ylim([0, nrM])

    axis.set_xlabel("l")
    axis.set_ylabel("m")

    #     nrTimes, nrBeams, nrGates = np.shape(currentData.data)
    #     ticks   = []
    #     labels  = []
    #     mod = int(np.floor(nrGates / 10))
    #     for x in xrange(nrGates):
    #         if x % mod != 0: continue
    #         ll = nrBeams*x
    #         ticks.append(ll)
    #         txt = '%i\n%i' % (ll, currentData.fov.gates[x])
    #         labels.append(txt)

    #     ticks.append(nrL)
    #     xlabels = copy.copy(labels)
    #     xlabels.append('l\ngate')

    #     axis.set_xticks(ticks)
    #     axis.set_xticklabels(xlabels,ha='left')

    #     ylabels = copy.copy(labels)
    #     ylabels.append('m\ngate')
    #     axis.set_yticks(ticks)
    #     axis.set_yticklabels(ylabels)

    xpos = 0.130
    fig.text(
        xpos, 0.99, "ABS(Cross Spectral Density Matrix Dlm)", fontsize=20, va="top"
    )
    # Get the time limits.
    #     timeLim = (np.min(currentData.time),np.max(currentData.time))
    #     md = currentData.metadata

    #     # Translate parameter information from short to long form.
    #     paramDict = getParamDict(md['param'])
    #     param     = paramDict['param']
    #     cbarLabel = paramDict['label']

    #     text = md['name'] + ' ' + param.capitalize() + timeLim[0].strftime(' (%Y %b %d %H:%M - ') + timeLim[1].strftime('%Y %b %d %H:%M)')

    #     if md.has_key('fir_filter'):
    #         filt = md['fir_filter']
    #         if filt[0] is None:
    #             low = 'None'
    #         else:
    #             low = '%.2f' % (1000. * filt[0])
    #         if filt[1] is None:
    #             high = 'None'
    #         else:
    #             high = '%.2f' % (1000. * filt[1])

    #         text = text + '\n' + 'Digital Filter: [' + low + ', ' + high + '] mHz'

    # fig.text(xpos,0.95,text,fontsize=14,va='top')
    fig.savefig("out.png", bbox_inches="tight")
