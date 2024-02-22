#!/usr/bin/env python3

"""simulate.py: simulate python program for RT"""

__author__ = "Chakraborty, S."
__copyright__ = ""
__credits__ = []
__license__ = "MIT"
__version__ = "1.0."
__maintainer__ = "Chakraborty, S."
__email__ = "shibaji7@vt.edu"
__status__ = "Research"

import argparse
import datetime as dt
from loguru import logger
from dateutil import parser as dparser
from pathlib import Path
import glob
import os
from rt2D import (
    execute_gemini2D_simulation, execute_gemini2D_simulations,
    execute_iri2D_simulations
)
from plots import gerenate_fov_plot
from rt3D import execute_iri3D_simulations

def clean():
    files = glob.glob(str(Path.home() / "matlab_crash_dump*"))
    for f in files:
        if os.path.exists(f): os.remove(f)
    return



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--model", default="gemini2D", help="Model name [iri(2D/3D) or gemini(2D/3D)]"
    )
    parser.add_argument(
        "-md", "--method", default="rt", help="Method rt/rti/movie"
    )
    parser.add_argument("-r", "--rad", default="fhe", help="Radar code (default fhe)")
    parser.add_argument(
        "-bm", "--beam", default=11, type=int, help="Radar beam (default 3)"
    )
    parser.add_argument(
        "-tsmin", "--time_steps_min", default=-1, type=int, help="Time steps for the RT simulation"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2016, 7, 8, 1, 30),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    parser.add_argument(
        "-sct", "--scatter_type", default="gs", help="scatter type - gs/is",
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    gerenate_fov_plot(args.rad, args.beam, args.event)
    if args.model == "gemini2D":
        execute_gemini2D_simulations(args)
    if args.model == "iri2D":
        execute_iri2D_simulations(args)
    if args.model == "iri3D":
        execute_iri3D_simulations(args)
