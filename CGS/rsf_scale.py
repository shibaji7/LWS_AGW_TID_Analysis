#!/usr/bin/env python
"""Make radar-style Es panels from a digisonde GRM file."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.dates as mdates
import numpy as np
import pandas as pd
from loguru import logger
import datetime as dt

from pynasonde.digisonde.scaler import normalize_echo_dataframe

# MHJ45 BC840 AL945
DEFAULT_STN = "BC840"
DEFAULT_GRM = Path(f"/media/chakras4/Ionosonde/{DEFAULT_STN}_20170527(147).GRM")
DEFAULT_OUTFIG = Path(__file__).with_name(f"{DEFAULT_STN}_es_digisonde.png")
frame_slice = slice(None, None, 2)


def read_grm_dataframe(grm_file: Path, n_workers: int = 1, frame_slice: slice=None) -> pd.DataFrame:
    """Read a GIRO GRM archive into one concatenated pynasonde dataframe."""
    from pynasonde.digisonde.parsers.grm import GrmSplitter

    splitter = GrmSplitter(grm_file)
    return splitter.load_dataframes(n_workers=n_workers, frame_slice=frame_slice), splitter.fmt



def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate VHF-radar-style Es panels from a digisonde GRM file."
    )
    parser.add_argument("--grm", type=Path, default=DEFAULT_GRM)
    parser.add_argument("--outfig", type=Path, default=DEFAULT_OUTFIG)
    parser.add_argument("--grm-format", default="MMM", choices=["RSF", "SBF", "MMM"])
    parser.add_argument("--all-mode", action="store_true")
    parser.add_argument("--workers", type=int, default=16)
    parser.add_argument("--grm-time-bin", default="5min")
    parser.add_argument("--grm-height-bin-km", type=float, default=5.0)
    parser.add_argument("--log-level", default="WARNING")
    args = parser.parse_args()

    logger.remove()
    # logger.add(sys.stderr, level=args.log_level.upper())

    grm_df, fmt = read_grm_dataframe(args.grm, n_workers=args.workers, frame_slice=frame_slice)
    print(f"GRM dataframe shape: {grm_df.shape}")
    print(f"GRM dataframe columns: {list(grm_df.columns)}")

    


if __name__ == "__main__":
    main()
