#!/usr/bin/env python

"""rays.py: Calculate all the functions of utility"""

__author__ = "Chakraborty, S."
__copyright__ = "Chakraborty, S."
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "chakras4@erau.edu"
__status__ = "Research"

import datetime as dt
import glob
import os
from types import SimpleNamespace

import matplotlib.colors as colors
import numpy as np
import pandas as pd
import tidUtils
from geopy.distance import great_circle as GC
from loguru import logger
from scipy.io import loadmat


def load_from_file(to_file: str):
    logger.info(f"Load from file {to_file.split('/')[-1]}")
    param = loadmat(to_file)["ne"]
    return param


def load_bearing_mat_file(file_loc: str):
    logger.info(f" Loading bearing file: {file_loc}")
    bearing = SimpleNamespace(**loadmat(file_loc))
    return bearing


def load_rays_mat_file(file_loc: str, limit_elvs=[]):
    logger.info(f" Loading rays file: {file_loc}")
    sim_data = loadmat(file_loc)
    path_data_keys = [
        "ground_range",
        "height",
        "group_range",
        "phase_path",
        "geometric_distance",
        "electron_density",
        "refractive_index",
    ]
    ray_data_keys = [
        "ground_range",
        "group_range",
        "phase_path",
        "geometric_path_length",
        "initial_elev",
        "final_elev",
        "apogee",
        "gnd_rng_to_apogee",
        "plasma_freq_at_apogee",
        "virtual_height",
        "effective_range",
        "deviative_absorption",
        "TEC_path",
        "Doppler_shift",
        "Doppler_spread",
        "frequency",
        "nhops_attempted",
        "ray_label",
    ]
    ray_data, ray_path_data = [], dict()
    for i in range(sim_data["ray_data"].shape[1]):
        r_data, p_data = dict(), dict()
        for key in ray_data_keys:
            r_data[key] = sim_data["ray_data"][0, i][key].ravel()[0]
            if key == "initial_elev":
                e = r_data[key]
            if key == "ray_label":
                rl = r_data[key]
            if key == "plasma_freq_at_apogee":
                pfa = r_data[key]
        for key in path_data_keys:
            p_data[key] = sim_data["ray_path_data"][0, i][key].ravel()
        p_data["elv"] = e
        p_data["plasma_freq_at_apogee"] = pfa
        p_data["ray_label"] = rl
        add = (
            (True if e >= limit_elvs[0] and e <= limit_elvs[1] else False)
            if len(limit_elvs) == 2
            else True
        )
        if add:
            ray_path_data[e] = pd.DataFrame.from_records(p_data)
            ray_data.append(r_data)
    ray_data = pd.DataFrame.from_records(ray_data)
    return ray_data, ray_path_data


def calc_relative_power(ray_data, labels=[1]):
    pwer = pd.DataFrame()
    o = ray_data.copy()
    o = o[o.ray_label.isin(labels)]
    o["weights"] = 1.0 / (o.group_range**3)
    ranges = 180 + 45 * np.arange(76, dtype=int)
    lag_power, bins = np.histogram(
        o.group_range,
        bins=ranges,
        weights=o.weights,
    )
    pwer["p_l"], pwer["srange"], pwer["slist"] = (
        lag_power,
        ranges[:-1],
        range(75),
    )
    pwer.p_l.replace(0, 1e-10, inplace=True)
    pwer["slist"] = (ranges[:-1] - 180) / 45
    px = 10 * np.log10(pwer["p_l"])
    px[px < -95] = np.nan
    pwer["p_l"] = px
    pwer["rsep"], pwer["frang"] = 45, 180
    return pwer


def create_lat_lon_from_routes(
    grange: np.array,
    r_bearing: float,
    olat: float,
    olon: float,
):
    lats, lons = [], []
    p = (olat, olon)
    gc = GC(p, p)
    for d in grange:
        x = gc.destination(p, r_bearing, distance=d)
        lats.append(x[0])
        lons.append(x[1])
    lats, lons = np.array(lats), np.array(lons)
    return lats, lons


def get_datasets_by_beams(
    rad,
    beams=None,
    start_time=None,
    end_time=None,
    limit_elvs=[],
    base_folder="/home/chakras4/OneDrive/trace/outputs/",
    run_name="May2017_gemini_tid_cosmic2",
    model_name="gemini",
):
    """
    Get datasets by beams
    """
    if beams is None:
        beams = [
            int(x.split("/")[-1])
            for x in glob.glob(
                os.path.join(
                    base_folder,
                    run_name,
                    f"{start_time.strftime('%Y-%m-%d')}",
                    f"{rad}",
                    "*",
                )
            )
        ]
    DS = pd.DataFrame()
    for b in beams:
        folder = os.path.join(
            base_folder,
            run_name,
            f"{start_time.strftime('%Y-%m-%d')}",
            f"{rad}",
            "%02d" % b,
            model_name,
        )
        bearing = load_bearing_mat_file(
            os.path.join(
                folder,
                f"bearing.mat",
            )
        )
        for d in range(int((end_time - start_time).total_seconds() / 60)):
            d = start_time + dt.timedelta(minutes=d)
            rays_file_loc = os.path.join(
                folder,
                f"{d.strftime('%H%M')}_rt.mat",
            )
            rays, _ = load_rays_mat_file(rays_file_loc, limit_elvs=limit_elvs)
            powr = calc_relative_power(rays)
            powr["bmnum"] = b
            powr["rad"] = rad
            powr["time"] = d
            DS = pd.concat([DS, powr])
    return DS


class RayTraceObject(object):

    def __init__(
        self,
        event,
        rad,
        beam,
        limit_elvs=[],
        base_folder="/home/chakras4/OneDrive/trace/outputs/",
        run_name="May2017_gemini_tid_cosmic2",
        model_name="gemini",
    ):
        self.rad = rad
        self.event = event
        self.beam = beam

        folder = os.path.join(
            base_folder,
            run_name,
            f"{event.strftime('%Y-%m-%d')}",
            f"{rad}",
            "%02d" % beam,
            model_name,
        )
        logger.debug(f"folder: {folder}")
        self.bearing = load_bearing_mat_file(
            os.path.join(
                folder,
                f"bearing.mat",
            )
        )
        self.edens = load_from_file(
            os.path.join(folder, f"{event.strftime('%H.%M')}.mat")
        )
        self.frequency = self.bearing.freq.ravel()[0] * 1e6
        self.pf = np.sqrt(80.6164e-6 * self.edens)
        self.ref_indx = np.sqrt(1 - (self.pf**2 / self.frequency**2))
        self.rays, self.ray_path = load_rays_mat_file(
            os.path.join(
                folder,
                f"{event.strftime('%H%M')}_rt.mat",
            ),
            limit_elvs=limit_elvs,
        )
        return

    def compile(self, kind="ray_path"):
        df = pd.DataFrame()
        if kind == "ray_path":
            df = pd.concat(
                [getattr(self, kind)[k] for k in list(getattr(self, kind).keys())]
            )
        return df


import matplotlib.pyplot as plt

# import scienceplots

plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = ["Tahoma", "DejaVu Sans", "Lucida Grande", "Verdana"]
plt.rcParams["text.usetex"] = False
import numpy as np


class PlotRays(object):
    def __init__(
        self,
        rto,
        nrows=2,
        ncols=2,
        ylim=[],
        xlim=[],
        xtolim=1700,
        lw=0.2,
        arc=False,
        figsize=(8, 3),
    ):
        self.nrows = nrows
        self.ncols = ncols
        self.rto = rto
        self.set_rto()
        self.xlim = xlim
        self.ylim = ylim
        self.axnum = 0
        self.fig = plt.figure(
            figsize=(figsize[0] * ncols, figsize[1] * nrows), dpi=1000
        )
        self.xtolim = xtolim
        self.lw = lw
        self.arc = arc
        return

    def set_rto(self):
        self.event = self.rto.event
        self.edens = self.rto.edens
        self.pf = self.rto.pf
        self.ref_indx = self.rto.ref_indx
        self.rad = self.rto.rad
        self.beam = self.rto.beam
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return

    def get_parameter(self, kind):
        import matplotlib.colors as colors

        if kind == "pf":
            o, cmap, label, norm = (
                getattr(self, kind),
                "GnBu",
                r"$f_0$ [MHz]",
                colors.Normalize(4.5, 5.5),
            )
        if kind == "edens":
            o, cmap, label, norm = (
                getattr(self, kind),
                "GnBu",
                r"$N_e$ [$/cm^{-3}$]",
                colors.LogNorm(1e5, 1e6),
            )
        if kind == "ref_indx":
            o, cmap, label, norm = (
                getattr(self, kind),
                "GnBu",
                r"$\eta$",
                colors.Normalize(0.8, 1),
            )
        return o, cmap, label, norm

    def get_arc_heights(self, height, dist):
        darc = dist / 6371.0
        true_height = 6371.0 + height
        height = true_height * np.cos(darc) - 6371.0
        return height

    def lay_rays(
        self,
        xlim_max=2000,
        kind="pf",
        zoomed_in=[],
        lcolor="k",
        tag_distance: float = -1,
        ax=None,
        xlabel=r"Ground range, km",
        ylabel=r"Height, km",
        add_time=True,
        add_cbar=True,
        add_tag=True,
        text="(A)",
        rto=None,
        ped_angles=[],
    ):
        self.rto = rto if rto else self.rto
        self.set_rto()
        ax = ax if ax else self.create_figure_pane(xlabel, ylabel)

        o, cmap, label, norm = self.get_parameter(kind)
        ## Check needs interpolation
        if self.xtolim <= xlim_max:
            index_xlim = np.abs(self.xtolim - self.rto.bearing.dist.ravel()).argmin()
            dist, height = (
                self.rto.bearing.dist.ravel(),
                self.rto.bearing.heights.ravel(),
            )
            d = self.rto.bearing.dist.ravel()[1] - self.rto.bearing.dist.ravel()[0]
            dist_ext = np.arange(self.xtolim, xlim_max, d)
            o_ext = np.array(
                [o[:, index_xlim].tolist() for i in range(len(dist_ext))]
            ).T
            dist = np.concatenate((dist[:index_xlim], dist_ext))
            o = np.concatenate((o[:, :index_xlim], o_ext), axis=1)
        else:
            dist, height = (
                self.rto.bearing.dist.ravel(),
                self.rto.bearing.heights.ravel(),
            )
        dist, height = np.meshgrid(dist, height)
        if self.arc:
            height = self.get_arc_heights(height, dist)
        im = ax.pcolormesh(
            dist,
            height,
            o,
            norm=norm,
            cmap=cmap,
            alpha=1.,
            zorder=3,
        )
        ds = dist[0, :]
        hs = []
        for i in range(o.shape[1]):
            hs.append(height[:, i][np.argmax(o[:, i])])
        ax.plot(ds, hs, ls="--", lw=0.3, color="k", zorder=3)
        ax.set_xlim(right=xlim_max)
        if add_cbar:
            pos = ax.get_position()
            cpos = [
                pos.x1 + 0.025,
                pos.y0 + 0.05,
                0.015,
                pos.height * 0.6,
            ]
            cax = self.fig.add_axes(cpos)
            cbax = self.fig.colorbar(
                im, cax, spacing="uniform", orientation="vertical", cmap="plasma"
            )
            _ = cbax.set_label(label, fontsize=11)
            cbax.ax.tick_params(axis="both", labelsize=11)
        rays = self.rto.rays
        self.elvs = rays.initial_elev
        if tag_distance > 100:
            ax.plot(
                [tag_distance, tag_distance],
                (
                    self.get_arc_heights(
                        np.array([0, 100]), np.array([tag_distance, tag_distance])
                    )
                    if self.arc
                    else [0, 100]
                ),
                c="k",
                zorder=4,
                alpha=0.7,
                ls="--",
                lw=1.2,
            )
        for i, elv in enumerate(self.elvs):
            ray_path_data, ray_data = (
                self.rto.ray_path[elv],
                rays[rays.initial_elev == elv],
            )
            th, r = (ray_path_data.ground_range.copy(), ray_path_data.height.copy())
            if self.arc:
                r = self.get_arc_heights(r, th)
                gr, h = self.get_height_range(ray_path_data)
                ax.scatter(gr, h, marker=".", s=0.1, color="k", zorder=6)
            ray_label = ray_data["ray_label"].iloc[0]
            lw = self.lw
            alpha = 0.3
            if ray_path_data.ray_label.iloc[0] == -2:
                lcolor = "r"
            elif ray_path_data.ray_label.iloc[0] == 1:
                lcolor = "k"
            if len(ped_angles) > 0:
                if np.round(ray_path_data.elv.iloc[0], 1) in ped_angles:
                    lcolor, alpha, lw = "darkgreen", 1, 1
            ax.plot(th, r, c=lcolor, zorder=3, alpha=alpha, ls="-", lw=lw)
            col = "k" if ray_label == 1 else "r"
            if ray_label in [-1, 1]:
                ax.scatter([th.iloc[-1]], [r.iloc[-1]], marker="s", s=3, color=col)
        if add_time:
            stitle = "%s UT" % self.event.strftime("%Y-%m-%d %H:%M")
            ax.text(
                0.95,
                1.05,
                stitle,
                ha="right",
                va="center",
                transform=ax.transAxes,
                fontdict={"size": 8, "fontweight": "bold"},
            )
        if add_tag:
            stitle = f"Model: GEMINI / {self.rad}-{'%02d'%self.beam}, $f_0$={self.rto.frequency/1e6} MHz"
            ax.text(
                0.05,
                1.05,
                stitle,
                ha="left",
                va="center",
                transform=ax.transAxes,
                fontdict={"size": 12, "fontweight": "bold"},
            )

        ax.text(
            0.05,
            0.95,
            text,
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict={"size": 12},
        )

        # Create Zoomed in panel
        if len(zoomed_in):
            self.__zoomed_in_panel__(ax, kind, zoomed_in, lcolor, ped_angles)
        ax.tick_params(axis="both", which="both", bottom=True, top=True, left=True, right=True)
        return ax

    def get_height_range(
        self, df, group=np.array([500, 750, 1000, 1250, 1500, 1750, 2000])
    ):
        gr, h = [], []
        for g in group:
            id = (df.group_range - g).abs().idxmin()
            if id < len(df) - 1:
                gr.append(df.ground_range.iloc[id])
                h.append(df.height.iloc[id])
            else:
                gr.append(np.nan)
                h.append(np.nan)
        h, gr = np.array(h), np.array(gr)
        h = self.get_arc_heights(h, gr)
        return gr, h

    def create_figure_pane(self, xlabel=r"Ground range, km", ylabel=r"Height, km"):
        self.axnum += 1
        fignum = 100 * self.nrows + 10 * self.ncols + self.axnum
        ax = self.fig.add_subplot(fignum)
        # Create Arc
        if self.arc:
            R = 6371.0
            theta = np.deg2rad(np.linspace(0, 90, 91))
            x, y = R * np.cos(theta), R * np.sin(theta) - R
            ax.plot(x, y, ls="-", color="k", lw=1)
            ax.text(
                -150,
                200,
                ylabel,
                ha="left",
                va="center",
                fontdict={"size": 12, "fontweight": "bold"},
                rotation=90,
            )
            ax.text(
                1000,
                -150,
                "Ground Range, km",
                ha="center",
                va="top",
                fontdict={"size": 12, "fontweight": "bold"},
            )
            ax.set_facecolor("0.98")
            ax.fill_between(x, -300 * np.ones_like(y), y, color="gray", alpha=0.5)
        else:
            ax.set_ylabel(ylabel, fontdict={"size": 12, "fontweight": "bold"})
            ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_xlim(self.xlim if len(self.xlim) == 2 else [0, 2000])
        ax.set_ylim(self.ylim if len(self.ylim) == 2 else [-300, 400])
        ax.tick_params(axis="both", labelsize=11)
        return ax

    def __zoomed_in_panel__(self, ax, kind, zoomed_in, lcolor="k", ped_angles=[]):
        self.zoom_ax = ax.inset_axes([0.5, 0.7, 0.5, 0.3])
        o, cmap, _, norm = self.get_parameter(kind)
        self.zoom_ax.pcolormesh(
            self.rto.bearing.dist.ravel(),
            self.rto.bearing.heights.ravel(),
            o,
            norm=norm,
            cmap=cmap,
            alpha=0.8,
        )
        rays = self.rto.rays
        self.elvs = rays.initial_elev
        for i, elv in enumerate(self.elvs):
            ray_path_data, ray_data = (
                self.rto.ray_path[elv],
                rays[rays.initial_elev == elv],
            )
            th, r = (ray_path_data.ground_range.copy(), ray_path_data.height.copy())
            lw = self.lw
            alpha = 0.3
            lcolor = "k"
            if len(ped_angles) > 0:
                if np.round(ray_path_data.elv.iloc[0], 1) in ped_angles:
                    lcolor, alpha, lw = "darkgreen", 1, 1
            self.zoom_ax.plot(th, r, c=lcolor, zorder=3, alpha=alpha, ls="-", lw=lw)

        self.zoom_ax.set_xlim(zoomed_in[0])
        self.zoom_ax.set_ylim(zoomed_in[1])
        th_ticklabels, r_ticklabels = (
            self.zoom_ax.get_xticklabels(),
            self.zoom_ax.get_yticklabels(),
        )
        # self.zoom_ax.set_xlabel("Ground Range, km", fontdict={"size": 8})
        # self.zoom_ax.set_ylabel("Height, km", fontdict={"size": 8})
        ax.indicate_inset_zoom(self.zoom_ax)
        return


class PlotChannels(object):
    def __init__(self, rto, nrows=2, ncols=2, ylim=[], xlim=[], xtolim=1700):
        self.nrows = nrows
        self.ncols = ncols
        self.rto = rto
        self.set_rto()
        self.xlim = xlim
        self.ylim = ylim
        self.axnum = 0
        self.fig = plt.figure(figsize=(8 * ncols, 3 * nrows), dpi=1000)
        self.xtolim = xtolim
        return

    def set_rto(self):
        self.event = self.rto.event
        self.edens = self.rto.edens
        self.pf = self.rto.pf
        self.ref_indx = self.rto.ref_indx
        self.rad = self.rto.rad
        self.beam = self.rto.beam
        return

    def save(self, filepath):
        self.fig.savefig(filepath, bbox_inches="tight", facecolor=(1, 1, 1, 1))
        return

    def close(self):
        self.fig.clf()
        plt.close()
        return

    def get_parameter(self, kind):

        if kind == "pf":
            o, cmap, label, norm = (
                getattr(self, kind),
                "Blues",
                r"$f_0$ [MHz]",
                colors.Normalize(4, 6),
            )
        if kind == "edens":
            o, cmap, label, norm = (
                getattr(self, kind),
                "Blues",
                r"$N_e$ [$/cm^{-3}$]",
                colors.LogNorm(1e5, 1e6),
            )
        if kind == "ref_indx":
            o, cmap, label, norm = (
                getattr(self, kind),
                "Blues",
                r"$\eta$",
                colors.Normalize(0.8, 9),
            )
        return o, cmap, label, norm

    def lay_rays(
        self,
        xlim_max=2000,
        kind="pf",
        zoomed_in=[],
        lcolor="k",
        tag_distance: float = -1,
        ax=None,
        xlabel=r"Ground range, km",
        ylabel=r"Height, km",
        add_time=True,
        add_cbar=True,
        add_tag=True,
        text="(A)",
        rto=None,
    ):
        self.rto = rto if rto else self.rto
        self.set_rto()
        df = self.rto.compile()
        ax = ax if ax else self.create_figure_pane(xlabel, ylabel)

        X, Y, Z = tidUtils.get_gridded_parameters(
            df,
            "elv",
            "geometric_distance",
            # "phase_path",
            "refractive_index",
            rounding=False,
        )
        im = ax.scatter(
            X.ravel(),
            Y.ravel(),
            c=Z.T.ravel(),
            s=50,
            marker="s",
            norm=colors.Normalize(0.82, 0.95),
            cmap="plasma",
            alpha=0.7,
        )
        ax.set_xlim(right=xlim_max)
        if add_cbar:
            pos = ax.get_position()
            cpos = [
                pos.x1 + 0.025,
                pos.y0 + 0.05,
                0.015,
                pos.height * 0.6,
            ]
            cax = self.fig.add_axes(cpos)
            cbax = self.fig.colorbar(
                im, cax, spacing="uniform", orientation="vertical", cmap="plasma"
            )
            _ = cbax.set_label(r"$\eta$", fontsize=11)
            cbax.ax.tick_params(axis="both", labelsize=11)
        if add_time:
            stitle = "%s UT" % self.event.strftime("%Y-%m-%d %H:%M")
            ax.text(
                0.95,
                1.05,
                stitle,
                ha="right",
                va="center",
                transform=ax.transAxes,
                fontdict={"size": 8, "fontweight": "bold"},
            )
        if add_tag:
            stitle = f"Model: GEMINI / {self.rad}-{'%02d'%self.beam}, $f_0$={self.rto.frequency/1e6} MHz"
            ax.text(
                0.05,
                1.05,
                stitle,
                ha="left",
                va="center",
                transform=ax.transAxes,
                fontdict={"size": 12, "fontweight": "bold"},
            )

        ax.text(
            0.05,
            0.95,
            text,
            ha="left",
            va="center",
            transform=ax.transAxes,
            fontdict={"size": 12},
        )
        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        ax.set_xlabel(xlabel, fontdict={"size": 12, "fontweight": "bold"})
        ax.set_ylabel(ylabel, fontdict={"size": 12, "fontweight": "bold"})
        # Create Zoomed in panel
        if len(zoomed_in):
            self.__zoomed_in_panel__(ax, kind, zoomed_in, lcolor)
        return ax

    def create_figure_pane(self, xlabel=r"Ground range, km", ylabel=r"Height, km"):
        self.axnum += 1
        fignum = 100 * self.nrows + 10 * self.ncols + self.axnum
        ax = self.fig.add_subplot(fignum)
        ax.tick_params(axis="both", labelsize=11)
        return ax

    def __zoomed_in_panel__(self, ax, kind, zoomed_in, lcolor="k"):
        self.zoom_ax = ax.inset_axes([0.4, 1.3, 0.3, 0.5])
        o, cmap, _, norm = self.get_parameter(kind)
        self.zoom_ax.pcolormesh(
            self.rto.bearing.dist.ravel(),
            self.rto.bearing.heights.ravel(),
            o,
            norm=norm,
            cmap=cmap,
            alpha=0.8,
        )
        rays = self.rto.rays
        self.elvs = rays.initial_elev
        for i, elv in enumerate(self.elvs):
            ray_path_data, ray_data = (
                self.rto.ray_path[elv],
                rays[rays.initial_elev == elv],
            )
            th, r = (ray_path_data.ground_range.copy(), ray_path_data.height.copy())
            self.zoom_ax.plot(th, r, c="k", zorder=3, alpha=0.7, ls="-", lw=0.5)

        self.zoom_ax.set_xlim(zoomed_in[0])
        self.zoom_ax.set_ylim(zoomed_in[1])
        th_ticklabels, r_ticklabels = (
            self.zoom_ax.get_xticklabels(),
            self.zoom_ax.get_yticklabels(),
        )
        self.zoom_ax.set_xlabel("Ground Range, km", fontdict={"size": 8})
        self.zoom_ax.set_ylabel("Height, km", fontdict={"size": 8})
        ax.indicate_inset_zoom(self.zoom_ax)
        return


import matplotlib as mpl

mpl.rc("font", size=15)
from typing import Optional, Sequence

import matplotlib.dates as mdates


class StackPlots:

    def __init__(
        self,
        nrows: int,
        ncols: int,
        dpi: int = 1000,
        datetime: bool = False,
        polar: bool = False,
        figsize: tuple = (8, 3),
    ):
        self.nrows = nrows
        self.ncols = ncols
        self.fig, self.axes = plt.subplots(
            self.nrows,
            self.ncols,
            figsize=(self.ncols * figsize[0], self.nrows * figsize[1]),
            sharex=True,
            dpi=dpi,
            subplot_kw={"projection": "polar"} if polar else {},
        )
        self.axes = self.axes.flatten() if self.nrows > 1 else [self.axes]
        self.fig.subplots_adjust(hspace=0.5, wspace=0.3)
        self.plot_id = 0
        self.datetime = datetime

    def plot_stack_plots(
        self,
        time: Sequence,
        value: Sequence,
        title: Optional[str] = None,
        xlabel: Optional[str] = None,
        ylabel: Optional[str] = None,
        label: Optional[str] = None,
        text: Optional[str] = None,
        ylim: Optional[Sequence] = None,
        xlim: Optional[Sequence] = None,
        color: str = "blue",
        lw: float = 0.8,
        ls: str = "-",
        ax: Optional[plt.Axes] = None,
        ylabel_color: Optional[str] = "k",
    ) -> tuple:
        """
        Plot a stack of plots with the given time and value data.
        :param time: Time data for the x-axis
        :param value: Value data for the y-axis
        :param title: Title for the plot
        :param xlabel: Label for the x-axis
        :param ylabel: Label for the y-axis
        """
        if self.plot_id > len(self.axes):
            raise ValueError("No more axes available for plotting.")
        if ax is None:
            ax = self.axes[self.plot_id]
            self.plot_id += 1
        if title:
            ax.set_title(title)
        if xlabel:
            ax.set_xlabel(xlabel)
        if ylabel:
            ax.set_ylabel(ylabel, color=ylabel_color)
        if ylim:
            ax.set_ylim(ylim)
        else:
            ax.set_ylim([min(value), max(value)])
        if xlim:
            ax.set_xlim(xlim)
        else:
            ax.set_xlim([time[0], time[-1]])
        if self.datetime:
            ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H"))
        ax.plot(time, value, color=color, linewidth=lw, ls=ls, label=label)
        if text:
            ax.text(0.05, 1.05, text, ha="left", va="center", transform=ax.transAxes)

        plt.tight_layout()
        return self.fig, ax

    def save_fig(self, filename: str):
        """
        Save the figure to a file.
        :param filename: Filename to save the figure
        """
        self.fig.savefig(filename, bbox_inches="tight")
        return

    def close(self):
        """
        Close the figure.
        """
        plt.close(self.fig)
        return

    def plot_dirctional_plots(
        self,
        theta: Sequence,
        r: Sequence,
        title: Optional[str] = None,
        text: Optional[str] = None,
        color: str = "black",
        ax: Optional[plt.Axes] = None,
        rlims: Optional[Sequence] = [0, 1],
        rticks: Optional[Sequence] = [0, 0.5, 1.0],
        theta_ticks: Optional[Sequence] = [0, np.pi / 2, np.pi, 3 * np.pi / 2],
        cable_angle: Optional[float] = None,
    ):
        """
        Plot directional plots.
        """
        if self.plot_id > len(self.axes):
            raise ValueError("No more axes available for plotting.")
        if ax is None:
            ax = self.axes[self.plot_id]
            ax.set_theta_zero_location("N")
            ax.set_theta_direction(-1)
            self.plot_id += 1
            if cable_angle is not None:
                ax.plot(
                    np.deg2rad([cable_angle, cable_angle]),
                    [0, 1],
                    lw=1.2,
                    ls="-",
                    color="m",
                )
        if title:
            ax.set_title(title, fontdict=dict(size=12))

        ax.bar(
            np.deg2rad(theta),
            r,
            bottom=0.0,
            color=color,
            width=np.deg2rad(theta[1] - theta[0]),
            alpha=0.5,
        )
        ax.set_rticks(rticks)
        ax.set_xticks(theta_ticks)
        ax.set_rmax(rlims[1])
        ax.set_rmin(rlims[0])
        if text:
            ax.text(
                -0.1,
                1.05,
                text,
                ha="left",
                va="center",
                transform=ax.transAxes,
                color=color,
            )
        return self.fig, ax
