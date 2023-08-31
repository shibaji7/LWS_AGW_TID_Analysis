#!/usr/bin/env python3

"""rt2D.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import argparse
import copy
import datetime as dt
import glob
import os

import numpy as np
import pandas as pd
import pydarn
from dateutil import parser as dparser
from geopy.distance import great_circle as GC
from loguru import logger
from scipy.io import savemat
from density_models_2D import GEMINI2D, IRI2D
from rays import Rays2D
import plots
import rtiPlots


def read_params_2D(fname="cfg/rt2D.json"):
    import json
    from types import SimpleNamespace

    with open(fname, "r") as f:
        param = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    return param

class RayTrace2D(object):
    """
    Ray trace class to trace all the points
    """

    def __init__(
        self,
        event,
        rad,
        beam,
        cfg,
    ):
        self.beam = beam
        self.event = event
        self.rad = rad
        self.folder = "simulation_results/{dn}/{rad}/".format(
            dn=self.event.strftime("%Y.%m.%d"), rad=self.rad
        )
        os.makedirs(self.folder, exist_ok=True)
        self.cfg = cfg
        self.hdw = pydarn.read_hdw_file(self.rad)
        self.fig_name = "{time}_{beam}.png".format(
            time=self.event.strftime("%H%M"),
            beam="%02d" % beam,
        )
        self._estimate_bearing_()
        return

    def _estimate_bearing_(self):
        """Estimate laitude and logitude bearings"""
        fname = self.folder + f"bearing_{'%02d'%self.beam}.mat"
        bearing = self.hdw.boresight.physical + (self.beam - self.hdw.beams/2) * self.hdw.beam_separation
        logger.info(f"Bearing angle of beam {self.beam} is {bearing} deg")
        lat, lon = (self.hdw.geographic.lat, self.hdw.geographic.lon)
        p = (lat, lon)
        gc = GC(p, p)
        dist = np.linspace(
            0, self.cfg.max_ground_range_km, self.cfg.number_of_ground_step_km
        )

        m = {}
        lats, lons = [], []
        for d in dist:
            x = gc.destination(p, bearing, distance=d)
            lats.append(x[0])
            lons.append(x[1])
        m["dist"], m["lat"], m["lon"] = dist, np.array(lats), np.array(lons)
        (
            m["olat"],
            m["olon"],
            m["rb"],
            m["num_range"],
            m["max_range"],
            m["range_inc"],
        ) = (
            lat,
            lon,
            bearing,
            float(len(dist)),
            float(self.cfg.max_ground_range_km),
            float(dist[1] - dist[0]),
        )
        m["ht"] = np.arange(
            self.cfg.start_height_km,
            self.cfg.end_height_km,
            self.cfg.height_incriment_km,
        ).astype(float)
        m["start_height"], m["height_inc"], m["num_heights"], m["heights"] = (
            float(self.cfg.start_height_km),
            float(self.cfg.height_incriment_km),
            float(
                len(
                    np.arange(
                        self.cfg.start_height_km,
                        self.cfg.end_height_km,
                        self.cfg.height_incriment_km,
                    )
                )
            ),
            np.arange(
                self.cfg.start_height_km,
                self.cfg.end_height_km,
                self.cfg.height_incriment_km,
            ),
        )

        m["freq"], m["tol"], m["nhops"] = (
            float(self.cfg.frequency),
            float(1e-7),
            float(self.cfg.nhops),
        )
        m["elev_s"], m["elev_i"], m["elev_e"] = (
            float(self.cfg.start_elevation),
            float(self.cfg.elevation_inctiment),
            float(self.cfg.end_elevation),
        )
        m["radius_earth"] = self.cfg.radius_earth
        savemat(fname, m)
        self.bearing_object = copy.copy(m)
        return

    def compile(self, density):
        """Compute RT using Pharlap"""
        self.density = density
        self.sim_fname = self.folder + "{date}.{bm}_rt.mat".format(
            bm="%02d" % self.beam, date=self.event.strftime("%H%M")
        )
        pwd = os.getcwd() + "/pharlap/pharlap_4.1.3/dat"
        cmd = "export DIR_MODELS_REF_DAT={pwd};\
                cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';fname='{fname}';bm={bm};\
                rt_2D;exit;\"".format(
            pwd=pwd,
            ut=self.event.strftime("%Y %m %d %H %M"),
            rad=self.rad,
            dic=self.folder,
            bm=self.beam,
            fname=self.sim_fname,
        )
        logger.info(f"Running command: {cmd}")
        os.system(cmd)
        logger.info("Data-Model comparison: reading rays....")
        self.rays = Rays2D.read_rays(self.event, self.rad, self.beam, self.cfg, self.folder, self.sim_fname)
        return
    
def execute_gemini2D_simulation(
    args,
):
    """
    Execute GEMINI2D simulation for ray tracing
    """
    cfg = read_params_2D()
    rtobj = RayTrace2D(args.event, args.rad, args.beam, cfg)
    fname = rtobj.folder + "{dn}_{bm}.mat".format(
        dn=args.event.strftime("%H.%M"), bm="%02d" % args.beam
    )
    logger.info(f"Create matlab file: {fname}")
    gem = GEMINI2D.create_object(
        args.event.strftime("%Y%m%d"),
        cfg,
        f"dataset/GEMINI3D/",
        "nsall",
        "grid.mat",
        args.event,
        rtobj.bearing_object["lat"],
        rtobj.bearing_object["lon"],
        rtobj.bearing_object["ht"],
        dlat=0.2, 
        dlon=0.2,
        to_file=fname,
    )
    rtobj.compile(gem.param_val)
    plots.plot_rays(
        rtobj.folder,
        rtobj.fig_name,
        rtobj,
        f"GEMINI2D + {args.rad.upper()}/{str(args.beam)}/{str(cfg.frequency)}",
        maxground = cfg.max_ground_range_km+10,
    )
    return

def execute_gemini2D_simulations(
    args,
):
    """
    Execute GEMINI2D simulation for ray tracing
    """
    cfg = read_params_2D()
    days = GEMINI2D.get_time_keys(args.event.strftime("%Y%m%d"), f"dataset/GEMINI3D/")
    gem = GEMINI2D(
        args.event.strftime("%Y%m%d"),
        cfg,
        f"dataset/GEMINI3D/",
        "nsall",
        "grid.mat",
    )
    folder = "simulation_results/{dn}/{rad}/".format(
            dn=args.event.strftime("%Y.%m.%d"), rad=args.rad
        )
    movie = False
    beam_soundings_rays = []
    for d in days:
        args.event = d
        rtobj = RayTrace2D(args.event, args.rad, args.beam, cfg)
        if not os.path.exists(rtobj.folder + rtobj.fig_name):
            fname = rtobj.folder + "{dn}_{bm}.mat".format(
                dn=args.event.strftime("%H.%M"), bm="%02d" % args.beam
            )
            movie = True
            logger.info(f"Create matlab file: {fname}")
            gem.fetch_dataset_by_locations(
                d, 
                rtobj.bearing_object["lat"],
                rtobj.bearing_object["lon"],
                rtobj.bearing_object["ht"], 
                dlat=0.2, 
                dlon=0.2,
                to_file=fname
            )
            rtobj.compile(gem.param_val)
            plots.plot_rays(
                rtobj.folder,
                rtobj.fig_name,
                rtobj,
                f"GEMINI2D + {args.rad.upper()}/{str(args.beam)}/{str(cfg.frequency)}",
                maxground = cfg.max_ground_range_km+10,
            )
        rt_name = folder + "{date}.{bm}_rt.mat".format(
            bm="%02d" % args.beam, date=d.strftime("%H%M")
        )
        beam_soundings_rays.append(
            Rays2D.read_rays(d, args.rad, args.beam, cfg, folder, rt_name)
        )
        if d>dt.datetime(2016, 7, 8, 3, 13, 12): break
    if movie: plots.create_movie(rtobj.folder, "{rad}_{bm}".format(rad=args.rad, bm="%02d" % args.beam))
    rtiPlots.create_RTI(rtobj.folder, beam_soundings_rays)
    return

def execute_iri2D_simulations(args):
    """
    Execute IRI2D simulation for ray tracing
    """
    cfg = read_params_2D()
    times = [args.event+dt.timedelta(minutes=i) for i in range(args.time_steps_min)]
    folder = "simulation_results/{dn}/{rad}/".format(
            dn=args.event.strftime("%Y.%m.%d"), rad=args.rad
        )
    beam_soundings_rays = []
    tid_prop = dict(
        lamb_x = 1000,
        v_x = 0.1,
        a_lim = [0.1, 1],
        t = 0,
    )
    movie = False
    for ti, d in enumerate(times):
        args.event = d
        tid_prop["t"] = ti*60
        rtobj = RayTrace2D(args.event, args.rad, args.beam, cfg)
        if not os.path.exists(rtobj.folder + rtobj.fig_name):
            fname = rtobj.folder + "{dn}_{bm}.mat".format(
                dn=args.event.strftime("%H.%M"), bm="%02d" % args.beam
            )
            movie = True
            logger.info(f"Create matlab file: {fname}")
            iri = IRI2D(
                d, 
                rtobj.bearing_object["lat"],
                rtobj.bearing_object["lon"],
                rtobj.bearing_object["ht"],
                to_file=fname,
                cfg=cfg,
                tid_prop=tid_prop,
            )
            rtobj.compile(iri.edens)
            rtobj.density = rtobj.density
            plots.plot_rays(
                rtobj.folder,
                rtobj.fig_name,
                rtobj,
                f"IRI2D + {args.rad.upper()}/{str(args.beam)}/{str(cfg.frequency)}",
                maxground = cfg.max_ground_range_km+10,
            )
        rt_name = folder + "{date}.{bm}_rt.mat".format(
            bm="%02d" % args.beam, date=d.strftime("%H%M")
        )
        beam_soundings_rays.append(
            Rays2D.read_rays(d, args.rad, args.beam, cfg, folder, rt_name)
        )
        
    if movie: plots.create_movie(rtobj.folder, "{rad}_{bm}".format(rad=args.rad, bm="%02d" % args.beam))
    rtiPlots.create_RTI(rtobj.folder, beam_soundings_rays)
    return