#!/usr/bin/env python

"""
    rtiUtils.py: module to plot RTI plots with various y-axis transformation
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
import datetime as dt

import matplotlib.dates as mdates
import model_vheight as mvh
import numpy as np
import utils
from pysolar.solar import get_altitude_fast
from scipy.stats import pearsonr
import tidUtils


class RTI(object):
    """
    Create plots for velocity, width, power, elevation angle, etc.
    """

    def __init__(
        self,
        nGates,
        drange,
        fov,
        xlim,
        fig_title=None,
        num_subplots=1,
        angle_th=100.0,
        vhm=None,
        ylim=[180, 3000],
    ):
        self.nGates = nGates
        self.drange = drange
        self.fov = fov
        self.num_subplots = num_subplots
        self._num_subplots_created = 0
        self.fig = plt.figure(figsize=(6, 3 * num_subplots), dpi=240)
        if fig_title:
            plt.suptitle(
                fig_title, x=0.075, y=0.99, ha="left", fontweight="bold", fontsize=15
            )
        self.angle_th = angle_th
        self.vhm = vhm
        self.ylim = ylim
        self.xlim = xlim
        return

    def addParamPlot(
        self,
        df,
        beam,
        title,
        vlim=[20, 50],
        xlabel="Time (UT)",
        zparam="p_l",
        label="Power (dB)",
        yscale="srange",
        cmap="jet",
        cbar=True,
        plot_fov=False,
        tec_details=None,
        ax=None,
        alpha=1,
    ):
        df = df[df.bmnum == beam]
        if yscale == "srange":
            yrange, ylab, frang = (
                self.nGates * df.rsep.tolist()[0] + df.frang.tolist()[0],
                "Slant Range (km)",
                df.frang.tolist()[0],
            )
        else:
            yrange, ylab, frang = (self.nGates, "Range Gates", 0)
        if self.vhm:
            yscale = "virtual_height"
            df["virtual_height"] = (
                [mvh.standard_vhm(s) for s in df.srange]
                if self.vhm["method"] == "standard"
                else [mvh.chisham_vhm(s) for s in df.srange]
            )
            yrange, ylab = (
                (
                    mvh.standard_vhm(
                        self.nGates * df.rsep.tolist()[0] + df.frang.tolist()[0]
                    )
                    if self.vhm["method"] == "standard"
                    else mvh.chisham_vhm(
                        self.nGates * df.rsep.tolist()[0] + df.frang.tolist()[0]
                    )
                ),
                "Virtual Height (km)",
            )
        if ax is None:
            ax = self._add_axis()
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
            hours = mdates.HourLocator(byhour=range(0, 24, 4))
            ax.xaxis.set_major_locator(hours)
            ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
            ax.set_xlim(
                [mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])]
            )
        ax.set_ylim(0, 3000)
        ax.set_ylabel(ylab, fontdict={"size": 12, "fontweight": "bold"})
        X, Y, Z = tidUtils.get_gridded_parameters(
            df, xparam="time", yparam=yscale, zparam=zparam, rounding=False
        )
        im = ax.pcolormesh(
            X,
            Y,
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
        if cbar:
            self._add_colorbar(self.fig, ax, im, label=label)
        if title:
            ax.set_title(title, loc="left", fontdict={"fontweight": "bold"})
        if plot_fov:
            self.overlay_sza(
                self.fov,
                ax,
                beam,
                [0, self.nGates],
                df.rsep.iloc[0],
                df.frang.iloc[0],
                yscale,
            )
        ax.set_ylim(self.ylim)
        ax.set_xlim(self.xlim)
        return ax, Y[:, 0]

    def overlay_sza(self, fov, ax, beam, gate_range, rsep, frang, yscale):
        """
        Add terminator to the radar
        """
        times = [
            self.drange[0] + dt.timedelta(minutes=i)
            for i in range(int((self.drange[1] - self.drange[0]).total_seconds() / 60))
        ]
        R = 6378.1
        gates = np.arange(gate_range[0], gate_range[1])
        dn_grid = np.zeros((len(times), len(gates)))
        for i, d in enumerate(times):
            d = d.replace(tzinfo=dt.timezone.utc)
            for j, g in enumerate(gates):
                gdlat, glong = fov[0][g, beam], fov[1][g, beam]
                angle = 90.0 - get_altitude_fast(gdlat, glong, d)
                dn_grid[i, j] = angle
        terminator = np.zeros_like(dn_grid)
        terminator[dn_grid > self.angle_th] = 1.0
        terminator[dn_grid <= self.angle_th] = 0.0
        if yscale == "srange":
            gates = frang + (rsep * gates)
        elif yscale == "virtual_height":
            mvh.standard_vhm(self.nGates * df.rsep.tolist()[0] + df.frang.tolist()[0])
        else:
            # TODO
            pass
        times, gates = np.meshgrid(times, gates)
        ax.pcolormesh(
            times.T,
            gates.T,
            terminator,
            lw=0.01,
            edgecolors="None",
            cmap="gray_r",
            vmax=2,
            vmin=0,
            shading="nearest",
            alpha=0.3,
        )
        return

    def ovearlay_TEC(
        self,
        tec,
        tec_times,
        beam,
        cbar=True,
        plot_fov=False,
        gate_range=None,
        rsep=None,
        frang=None,
        yscale="srange",
        label="TEC (TECu)",
        cmap=tidUtils.parula,
        vlim=[-0.2, 0.2],
        ax=None,
        xlabel="Time (UT)",
        cbar_xOff=0.0,
        alpha=1.0,
        tec_param="cdvTECgrid2",
    ):
        """
        Add TEC in the observations panel
        """
        frang = frang if frang else 180
        rsep = rsep if rsep else 45
        gate_range = gate_range if gate_range else [0, self.nGates]
        Yaxis = (
            [frang + i * rsep for i in range(self.nGates + 1)]
            if yscale == "srange"
            else np.arange(self.nGates + 1)
        )
        if ax is None:
            ax = self._add_axis()
            if yscale == "srange":
                yrange, ylab = (
                    self.nGates * rsep + frang,
                    "Slant Range (km)",
                )
            else:
                yrange, ylab = (self.nGates, "Range Gates")
            ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
            hours = mdates.HourLocator(byhour=range(0, 24, 4))
            ax.xaxis.set_major_locator(hours)
            ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
            ax.set_xlim(
                [mdates.date2num(self.drange[0]), mdates.date2num(self.drange[1])]
            )
            ax.set_ylim(frang, yrange)
            ax.set_ylabel(ylab, fontdict={"size": 12, "fontweight": "bold"})

        dset = tidUtils.fetch_tec_by_beam(tec, beam, tec_param)
        im = ax.pcolormesh(
            tec_times[: dset.shape[0]],
            Yaxis,
            dset.T,
            lw=0.01,
            edgecolors="None",
            cmap=cmap,
            snap=True,
            vmax=vlim[1],
            vmin=vlim[0],
            shading="auto",
            alpha=alpha,
        )
        if cbar:
            self._add_colorbar(self.fig, ax, im, label=label, xOff=cbar_xOff)
        if plot_fov:
            self.overlay_sza(
                self.fov,
                ax,
                beam,
                [0, self.nGates],
                rsep,
                frang,
                yscale,
            )
        ax.set_ylim(self.ylim)
        ax.set_xlim(self.xlim)
        return ax

    def _add_axis(self):
        self._num_subplots_created += 1
        ax = self.fig.add_subplot(self.num_subplots, 1, self._num_subplots_created)
        return ax

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

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return


def plot_SDTEC_TS(
    sdx,
    sdy,
    tecx,
    tecy,
    fname,
    xlim,
    sd_ylim=[20, 50],
    tec_ylim=[-0.3, 0.3],
    txt="",
):
    """
    plot SD versus TEC
    """
    fig = plt.figure(figsize=(12, 3), dpi=240)
    spec = fig.add_gridspec(2, 5)
    ax = fig.add_subplot(spec[:, :3])
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
    hours = mdates.HourLocator(byhour=range(0, 24, 4))
    ax.xaxis.set_major_locator(hours)
    ax.set_xlabel("Time (UT)", fontdict={"size": 12, "fontweight": "bold"})
    ax.plot(sdx, sdy, "ko", ls="None", ms=0.8)
    ax.set_ylim(sd_ylim)
    ax.set_xlim(xlim)
    ax.set_ylabel("Power (dB)", fontdict={"size": 12, "fontweight": "bold"})
    ax.text(
        0.05,
        1.05,
        txt,
        ha="left",
        va="center",
        fontdict={"size": 10},
        transform=ax.transAxes,
    )

    ax = ax.twinx()
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter("$%H^{%M}$"))
    hours = mdates.HourLocator(byhour=range(0, 24, 4))
    ax.xaxis.set_major_locator(hours)
    ax.plot(tecx, tecy, "ro", ls="None", ms=0.8)
    ax.set_ylim(tec_ylim)
    ax.set_xlim(xlim)
    ax.set_ylabel(
        "TEC (TECu)", fontdict={"size": 12, "fontweight": "bold", "color": "r"}
    )

    ax = fig.add_subplot(spec[:, 3:])
    from scipy import signal

    N = int(len(tecy) / len(sdy)) + 1
    tec = tidUtils.interpolate_nans(tecy[::N])
    sd = tidUtils.interpolate_nans(sdy[: len(tec)])
    corr, _ = pearsonr(sd, tec)
    ax.set_ylim(tec_ylim)
    ax.set_xlim(sd_ylim)
    ax.plot(sd, tec, "ko", ls="None", ms=0.8)
    ax.set_ylabel("TEC (TECu)", fontdict={"size": 12, "fontweight": "bold"})
    ax.set_xlabel("Power (dB)", fontdict={"size": 12, "fontweight": "bold"})
    ax.plot([0, 1], [0, 1], transform=ax.transAxes, ls="--", color="gray", lw=0.8)
    ax.text(
        0.95,
        1.05,
        r"$\rho=%.2f$" % corr,
        ha="right",
        va="center",
        fontdict={"size": 10},
        transform=ax.transAxes,
    )

    fig.subplots_adjust(wspace=1)
    fig.savefig(fname, bbox_inches="tight", facecolor=(1, 1, 1, 1))
    fig.clf()
    plt.close()
    return


def plot_FFT(o, beam, srange, fname):
    """
    Plot FFT dataset
    """
    fig = plt.figure(figsize=(6, 3), dpi=240)
    ax = fig.add_subplot(111)
    ax.loglog(o.freq, o.power, "b-", lw=0.8)
    ax.set_xlim(1e-6, 1e-2)
    ax.set_ylim(1e-3, 1)
    ax.set_xlabel("Frequency, Hz")
    ax.set_ylabel("Power, Hz")
    ax.text(0.1, 1.05, "Range Cell: [%d,%d]"%(beam, srange), va="center", ha="left", transform=ax.transAxes)
    fig.savefig(fname, bbox_inches="tight", facecolor=(1, 1, 1, 1))
    fig.clf()
    plt.close()
    return