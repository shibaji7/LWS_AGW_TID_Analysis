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
DEFAULT_STN = "MHJ45"
DEFAULT_GRM = Path(f"/media/chakras4/Ionosonde/{DEFAULT_STN}_20170527(147).GRM")
DEFAULT_OUTFIG = Path(__file__).with_name(f"{DEFAULT_STN}_es_digisonde.png")
frame_slice = slice(None, None, 5)


def read_grm_dataframe(grm_file: Path, n_workers: int = 1, frame_slice: slice=None) -> pd.DataFrame:
    """Read a GIRO GRM archive into one concatenated pynasonde dataframe."""
    from pynasonde.digisonde.parsers.grm import GrmSplitter

    splitter = GrmSplitter(grm_file)
    return splitter.load_dataframes(n_workers=n_workers, frame_slice=frame_slice), splitter.fmt


def select_o_mode_or_fallback(df: pd.DataFrame, strict_o_mode: bool):
    """Select O-mode rows when present, otherwise use D256-compatible rows."""
    if "polarization" not in df.columns:
        if strict_o_mode:
            raise ValueError("GRM dataframe has no polarization column.")
        return df, "all polarizations"

    polarizations = set(df["polarization"].dropna().astype(str).unique())
    if "O" in polarizations:
        return df[df["polarization"].astype(str) == "O"], "O-mode"

    fallback = {"D256", "D128"}
    if fallback & polarizations:
        if strict_o_mode:
            raise ValueError(
                "No explicit O-mode rows were found. Available polarizations: "
                f"{sorted(polarizations)}"
            )
        return df[df["polarization"].astype(str).isin(fallback)], "D256-compatible"

    if strict_o_mode:
        raise ValueError(
            "No explicit O-mode rows were found. Available polarizations: "
            f"{sorted(polarizations)}"
        )
    return df, "all polarizations"

def _grid_edges(centers):
    """Return pcolormesh cell edges from sorted cell centers."""
    import numpy as np

    centers = np.asarray(centers, dtype=float)
    if len(centers) == 1:
        return np.array([centers[0] - 0.5, centers[0] + 0.5])

    midpoints = (centers[1:] + centers[:-1]) / 2
    first = centers[0] - (midpoints[0] - centers[0])
    last = centers[-1] + (centers[-1] - midpoints[-1])
    return np.concatenate([[first], midpoints, [last]])


def parse_day_time(day: pd.Timestamp, hhmm: str) -> pd.Timestamp:
    """Build a timestamp from the data day and an HH:MM string."""
    hour, minute = [int(part) for part in hhmm.split(":")]
    return day.replace(hour=hour, minute=minute, second=0, microsecond=0)

def interpolate_grid_for_display(grid: pd.DataFrame, height_step_km: float) -> pd.DataFrame:
    """Interpolate height axis for display only; time cadence is unchanged."""
    if grid.empty:
        return grid
    new_heights = np.arange(grid.index.min(), grid.index.max() + height_step_km, height_step_km)
    return grid.reindex(grid.index.union(new_heights)).sort_index().interpolate(
        method="index", limit_direction="both"
    ).reindex(new_heights)


def _select_o_mode_or_fallback(df: pd.DataFrame, strict_o_mode: bool):
    """Select O-mode rows when present, otherwise fall back to D256-compatible rows."""
    if "polarization" not in df.columns:
        if strict_o_mode:
            raise ValueError("GRM dataframe has no polarization column.")
        return df, "all polarizations"

    polarizations = set(df["polarization"].dropna().astype(str).unique())
    if "O" in polarizations:
        return df[df["polarization"].astype(str) == "O"], "O-mode"

    d256_labels = {"D256", "D128"}
    if d256_labels & polarizations:
        if strict_o_mode:
            raise ValueError(
                "No explicit O-mode rows were found. Available polarizations: "
                f"{sorted(polarizations)}"
            )
        return df[df["polarization"].astype(str).isin(d256_labels)], "D256-compatible"

    if strict_o_mode:
        raise ValueError(
            "No explicit O-mode rows were found. Available polarizations: "
            f"{sorted(polarizations)}"
        )
    return df, "all polarizations"

def plot_grm_o_mode_power_bins(
    grm_df: pd.DataFrame,
    outfile: Path,
    time_bin: str = "15min",
    height_bin_km: float = 10.0,
    strict_o_mode: bool = False,
    frequency_bins: tuple[tuple[float, float], ...] = (
        (1.5, 3.5),
        (3.0, 5.0),
    ),
    switch_param:bool = False,
) -> Path:
    """Plot GRM amplitude in frequency bands on a 2x2 pcolormesh."""
    import matplotlib.pyplot as plt
    import numpy as np
    from matplotlib import colormaps
    from pynasonde.digisonde.digi_utils import setsize

    setsize(14)

    grm_df = normalize_echo_dataframe(grm_df)
    print(grm_df.columns)
    if switch_param:
        grm_df["range_km"] = grm_df["height_km"]
        grm_df["amplitude_dB"] = grm_df["amplitude_db"]
        grm_df["polarization"] = grm_df["mode"]
        grm_df["datetime"] = grm_df["date"]

    required = {"datetime", "frequency_mhz", "range_km", "amplitude_dB"}
    missing = required.difference(grm_df.columns)
    if missing:
        raise KeyError(f"GRM dataframe is missing required columns: {sorted(missing)}")

    df = grm_df.copy()[list(required)]
    df["datetime"] = pd.to_datetime(df["datetime"])

    df, pol_label = _select_o_mode_or_fallback(df, strict_o_mode)

    df["time_bin"] = df["datetime"].dt.floor(time_bin)
    df["range_bin_km"] = (
        (df["range_km"] / height_bin_km).round() * height_bin_km
    ).astype(float)
    day_start = df["datetime"].min().normalize()
    day_end = day_start + pd.Timedelta(days=1)
    day_start += dt.timedelta(hours=12)

    grids = []
    values_for_scale = []
    for fmin, fmax in frequency_bins:
        band = df[(df["frequency_mhz"] >= fmin) & (df["frequency_mhz"] < fmax)]
        if band.empty:
            grids.append(None)
            continue

        grid = (
            band.groupby(["time_bin", "range_bin_km"], observed=True)["amplitude_dB"]
            .median()
            .unstack("range_bin_km")
            .sort_index()
            .sort_index(axis=1)
        )
        grids.append(grid)
        values_for_scale.append(grid.to_numpy().ravel())

    if not values_for_scale:
        raise ValueError("No GRM rows found inside the requested frequency bins.")

    all_values = np.concatenate(values_for_scale)
    all_values = all_values[np.isfinite(all_values)]
    vmin, vmax = np.nanpercentile(all_values, [2, 98])
    # vmin, vmax = 0, 10

    cmap = colormaps["jet"].copy()
    cmap.set_bad("0.88")

    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharex=True, sharey=True)
    fig.subplots_adjust(top=0.72, wspace=0.18)
    axes = np.atleast_1d(axes)
    mesh = None
    for ax, (fmin, fmax), grid in zip(axes.flat, frequency_bins, grids):
        ax.set_title(f"{fmin:g}-{fmax:g} MHz")
        if grid is None or grid.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            continue

        time_edges = _grid_edges(mdates.date2num(grid.index.to_pydatetime()))
        range_edges = _grid_edges(grid.columns.to_numpy(dtype=float))
        mesh = ax.pcolormesh(
            time_edges,
            range_edges,
            grid.to_numpy().T,
            shading="auto",
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            edgecolors="none",
            linewidth=0,
            antialiased=False,
        )
        ax.xaxis_date()
        ax.set_xlim(day_start, day_end)
        ax.xaxis.set_major_locator(mdates.HourLocator(interval=6))
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    for ax in axes.flat:
        ax.set_ylim(90, 300)
        ax.set_xlabel("Time, UT")
    axes.flat[0].set_ylabel("Virtual height, km")

    start_time = df["datetime"].min().strftime("%Y-%m-%d %H:%M:%S")
    end_time = df["datetime"].max().strftime("%Y-%m-%d %H:%M:%S")
    fig.suptitle(
        f"{DEFAULT_STN} GRM {pol_label} power by frequency band\n"
        f"{df['datetime'].min().strftime('%Y-%m-%d')}"
    )
    if mesh is not None:
        cbar = fig.colorbar(mesh, ax=axes.ravel().tolist(), pad=0.02)
        cbar.set_label("Amplitude, dB")

    fig.savefig(outfile, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return outfile


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

    grm_power_path = plot_grm_o_mode_power_bins(
        grm_df,
        args.outfig,
        time_bin=args.grm_time_bin,
        height_bin_km=args.grm_height_bin_km,
        strict_o_mode=not args.all_mode,
        switch_param=not (fmt.upper()=="MMM")
    )
    print(f"Saved GRM power plot: {grm_power_path}")


if __name__ == "__main__":
    main()
