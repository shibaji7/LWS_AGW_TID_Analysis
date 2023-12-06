#!/usr/bin/env python3

"""rt3D.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import os
from loguru import logger
from scipy.io import savemat
import numpy as np
from density_models_3D import Iono3D
import datetime as dt

import pydarn


def get_point_at_distance(lat, lon, d, bearing, R=6371):
    from math import asin, atan2, cos, degrees, radians, sin
    lat = radians(lat)
    lon = radians(lon)
    a = radians(bearing)
    lat_d = asin(sin(lat) * cos(d/R) + cos(lat) * sin(d/R) * cos(a))
    lon_d = lon + atan2(
        sin(a) * sin(d/R) * cos(lat),
        cos(d/R) - sin(lat) * sin(lat_d)
    )
    return (degrees(lat_d), degrees(lon_d))

def read_params_3D(fname="cfg/rt3D.json"):
    import json
    from types import SimpleNamespace

    with open(fname, "r") as f:
        param = json.load(f, object_hook=lambda x: SimpleNamespace(**x))
    return param

class RayTrace3D(object):
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
        self.folder = "simulation_results_3D/{dn}/{rad}/".format(
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
        self.sim_fname = self.folder + "rt_{date}.{bm}.mat".format(
            bm="%02d" % self.beam, date=self.event.strftime("%H%M")
        )
        fname = self.folder + f"bearing_{'%02d'%self.beam}.mat"
        bearing = self.hdw.boresight.physical + (self.beam - self.hdw.beams/2) * self.hdw.beam_separation
        logger.info(f"Bearing angle of beam {self.beam} is {bearing} deg")
        elevs = np.arange(
            self.cfg.start_elevation, 
            self.cfg.end_elevation, 
            self.cfg.elevation_inctiment
        )
        (lat_start, lon_start, ht_start) = (
            self.hdw.geographic.lat-(self.cfg.lat_inc*self.cfg.num_lat/5),
            self.hdw.geographic.lon-(self.cfg.lon_inc*self.cfg.num_lon/5),
            self.cfg.ht_start,
        )
        lat_d, lon_d = get_point_at_distance(lat_start, lon_start, self.cfg.distance_travel_km, bearing)
        logger.info(f"Difference in degrees: {np.round(lat_start-lat_d,2)}, {np.round(lon_d-lon_start,2)}")
        num_lat, num_lon = (
            int(abs(lat_start-lat_d)/self.cfg.lat_inc) + 1,
            int(abs(lon_start-lon_d)/self.cfg.lon_inc) + 1
        )
        num_lat = num_lat + 1 if np.mod(num_lat, 2) == 0 else num_lat
        num_lon = num_lon + 1 if np.mod(num_lon, 2) == 0 else num_lon
        logger.info(f"Number of latitudes/longitudes: {num_lat}/{num_lon}")
        self.bearing_obj = dict(
            doppler_flag = self.cfg.doppler_flag,
            speed_of_light = self.cfg.speed_of_light,
            R12 = self.cfg.R12,
            radius_earth = self.cfg.radius_earth,
            OX_mode = self.cfg.OX_mode,
            tol = [float(1e-7), 0.01, 25],
            nhops = float(self.cfg.nhops),
            
            origin_lat = self.hdw.geographic.lat,
            origin_long = self.hdw.geographic.lon,
            origin_ht = 0,
            elevs = elevs.astype(float),
            freqs = np.ones(len(elevs))*self.cfg.frequency,
            ray_bears = np.ones(len(elevs))*bearing,
            
            lat_start = lat_start,
            lat_inc = float(self.cfg.lat_inc),
            num_lat = float(num_lat),
            lon_start = lon_start,
            lon_inc = float(self.cfg.lon_inc),
            num_lon = float(num_lon),
            ht_start = float(self.cfg.ht_start),
            ht_inc = float(self.cfg.ht_inc),
            num_ht = float(self.cfg.num_ht),
            
            B_lat_start = lat_start, 
            B_lat_inc = float(self.cfg.B_lat_inc),
            B_lon_start = lon_start,
            B_lon_inc = float(self.cfg.B_lon_inc),
            B_ht_start = float(self.cfg.ht_start),
            B_ht_inc = float(self.cfg.B_ht_inc),
            
            beam = self.beam,
            rad = self.rad,
            fname = self.sim_fname,
        )
        savemat(fname, self.bearing_obj)
        return

    def compile(self, density):
        """Compute RT using Pharlap"""
        self.density = density        
        pwd = os.getcwd() + "/pharlap/pharlap_4.5.3/dat"
        cmd = "cd pharlap;\
                matlab -nodisplay -nodesktop -nosplash -nojvm -r \"UT=[{ut}];rad='{rad}';dic='{dic}';bm={bm};\
                rt_3D;exit;\"".format(
            ut=self.event.strftime("%Y %m %d %H %M"),
            rad=self.rad,
            dic=self.folder,
            bm=self.beam,
        )
        logger.info(f"Running command: {cmd}")
        if not os.path.exists(self.sim_fname):
            os.system(cmd)
        logger.info("Data-Model comparison: reading rays....")
        #self.rays = Rays2D.read_rays(self.event, self.rad, self.beam, self.cfg, self.folder, self.sim_fname)
        return
    
    
def execute_iri3D_simulations(args):
    """
    Run ray trace based on IRI model
    """
    cfg = read_params_3D()
    for i in range(args.time_steps_min):
        tid_prop = dict(
            lamb_x = 1000,
            v_x = 0.1,
            a_lim = [0.1, 1],
            t = i,
            source = "center"
        )
        event = args.event + dt.timedelta(minutes=i)
        rt = RayTrace3D(
            event,
            args.rad,
            args.beam,
            cfg,
        )
        fname = rt.folder + "{dn}_{bm}.mat".format(
            dn=event.strftime("%H.%M"), bm="%02d" % args.beam
        )
        iri = Iono3D(
            event, 
            rt.bearing_obj, 
            to_file=fname, 
            cfg=cfg, 
            tid_prop = tid_prop,
            iri_version=20,
        )
        rt.compile(None)
        del rt
        del iri
        if i==2: break