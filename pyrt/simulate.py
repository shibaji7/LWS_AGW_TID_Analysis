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
from rt2D import execute_gemini2D_simulation, execute_gemini2D_simulations


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", "--model", default="gemini2D", help="Model name [gitm/waccmx]"
    )
    parser.add_argument("-r", "--rad", default="fhe", help="Radar code (default fhe)")
    parser.add_argument(
        "-bm", "--beam", default=3, type=int, help="Radar beam (default 7)"
    )
    parser.add_argument(
        "-ev",
        "--event",
        default=dt.datetime(2016, 7, 8, 3, 31, 12),
        help="Event date for simulation [YYYY-mm-ddTHH:MM]",
        type=dparser.isoparse,
    )
    args = parser.parse_args()
    logger.info("\n Parameter list for simulation ")
    for k in vars(args).keys():
        print("     ", k, "->", str(vars(args)[k]))
    if args.model == "gemini2D":
        execute_gemini2D_simulations(args)
