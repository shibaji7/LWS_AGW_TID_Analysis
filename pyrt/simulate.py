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
from dateutil import parser as dparser
from rt2D import *


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--model", default="gemini2D", help="Model name [gitm/waccmx]"
    )
    parser.add_argument("-r", "--rad", default="fhe", help="Radar code (default fhe)")
    parser.add_argument(
        "-bm", "--beam", default=11, type=int, help="Radar beam (default 11)"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2017, 8, 21, 17, 30),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    if args.model == "gemini2D":
        cfg = read_params_2D()
        rtobj = RayTrace2D(args.event, args.rad, args.beam, cfg)
        execute_gemini2D_simulations(rtobj, cfg, args)
