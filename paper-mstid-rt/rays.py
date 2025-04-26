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

from types import SimpleNamespace
import os

import numpy as np
import pandas as pd
from geopy.distance import great_circle as GC
from loguru import logger
from scipy.io import loadmat
import datetime as dt

import glob

def load_bearing_mat_file(file_loc: str):
    logger.info(f" Loading bearing file: {file_loc}")
    bearing = SimpleNamespace(**loadmat(file_loc))
    return bearing


def load_rays_mat_file(file_loc: str):
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
        for key in path_data_keys:
            p_data[key] = sim_data["ray_path_data"][0, i][key].ravel()
        p_data["elv"] = e
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
    pwer["lag_power"], pwer["srange"], pwer["slist"] = (
        lag_power,
        ranges[:-1],
        range(75),
    )
    pwer.replace(0, 1e-10, inplace=True)
    pwer["lag_power"] = 10 * np.log10(pwer["lag_power"])
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
    base_folder="/media/chakras4/Crucial X9/trace/",
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
                    f"{rad}", "*"
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
            "%02d"%b,
            model_name,
        )
        bearing = load_bearing_mat_file(
            os.path.join(
                folder,
                f"bearing.mat",
            )
        )
        for d in range(int((end_time-start_time).total_seconds()/60)):
            d = start_time + dt.timedelta(minutes=d)
            rays_file_loc = os.path.join(
                folder,
                f"{d.strftime('%H%M')}_rt.mat",
            )
            rays, _ = load_rays_mat_file(rays_file_loc)
            powr = calc_relative_power(rays)
            powr["bmnum"] = b
            powr["rad"] = rad
            powr["time"] = d
            DS = pd.concat([DS,powr])
    return DS
