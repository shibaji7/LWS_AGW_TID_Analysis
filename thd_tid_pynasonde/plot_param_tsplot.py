"""Time-series plot of one scaled parameter across a full date range.

Plots the complete time series (native sounding cadence, no resampling) for
one CONUS station over an arbitrary run of days, reading the daily DIDBASE
CSVs downloaded by ``didbase_pipeline.py``
(``Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.csv``). The plotted
parameter defaults to foF2 but any column in those CSVs can be selected.

Run from this directory, for example:

    RSF_TSPLOT_STATION=AU930 \
    RSF_TSPLOT_PARAM=foF2 \
    RSF_TSPLOT_START_DATE=2024-04-01 \
    RSF_TSPLOT_N_DAYS=183 \
    /home/chakras4/anaconda3/envs/apep/bin/python -u plot_param_tsplot.py

Daily batch mode -- one figure per station per day per parameter (from
``config.toml``'s ``[tsplot_all]``), over every downloaded day -- also
``make plot-ts``:

    RSF_TSPLOT_DAILY=1 /home/chakras4/anaconda3/envs/apep/bin/python -u plot_param_tsplot.py

Full-season batch mode -- one figure per station per parameter, spanning the
whole downloaded date range instead of one-figure-per-day:

    RSF_TSPLOT_ALL=1 /home/chakras4/anaconda3/envs/apep/bin/python -u plot_param_tsplot.py
"""

from __future__ import annotations

import os
from datetime import date, timedelta
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

import utils

_CFG = utils.load_config()["tsplot"]
_ALL_CFG = utils.load_config()["tsplot_all"]
_SMOOTH_CFG = utils.load_config()["smoothing"]

STATION = os.environ.get("RSF_TSPLOT_STATION", "AU930")
DATA_DIR = Path(os.environ.get("RSF_TSPLOT_DATA_DIR", str(Path(__file__).resolve().parent / "Data" / "didbase_giro_scaled")))
OUTPUT_DIR = Path(os.environ.get("RSF_TSPLOT_OUTPUT_DIR", str(Path(__file__).resolve().parent / _CFG["output_dir"])))
START_DATE = pd.Timestamp(os.environ.get("RSF_TSPLOT_START_DATE", "2024-04-01"))
N_DAYS = int(os.environ.get("RSF_TSPLOT_N_DAYS", "30"))
PARAM = os.environ.get("RSF_TSPLOT_PARAM", _CFG["parameter"])

RUN_ALL = os.environ.get("RSF_TSPLOT_ALL", "0").strip().lower() in {"1", "true", "yes", "y", "on"}
RUN_DAILY = os.environ.get("RSF_TSPLOT_DAILY", "0").strip().lower() in {"1", "true", "yes", "y", "on"}
ALL_PARAMS = [p.strip() for p in os.environ.get(
    "RSF_TSPLOT_ALL_PARAMS", ",".join(_ALL_CFG["parameters"])).split(",") if p.strip()]
ALL_OUTPUT_DIR = Path(os.environ.get(
    "RSF_TSPLOT_ALL_OUTPUT_DIR", str(Path(__file__).resolve().parent / _ALL_CFG["output_dir"])))

# Pandas time-offset alias, e.g. "30min", "15min", "1h", "2h30min". Any
# arbitrary window works for any parameter -- this is not foF2/hmF2-specific.
# Empty string disables the smoothed overlay.
SMOOTH_WINDOW = os.environ.get("RSF_TSPLOT_SMOOTH_WINDOW", _SMOOTH_CFG["window"]).strip()


def _units_for(param: str) -> str:
    """DIDBASE naming convention: f*/foX = frequency [MHz], h*/hmX = height [km]."""
    p = param.lower()
    if p.startswith("h"):
        return "km"
    if p.startswith("f"):
        return "MHz"
    return ""


def smooth_series(df: pd.DataFrame, param: str, window: str = SMOOTH_WINDOW) -> pd.Series:
    """Time-based rolling mean of ``param`` over ``window``.

    ``window`` is any pandas time-offset alias (``"30min"``, ``"15min"``,
    ``"1h"``, ``"2h30min"``, ...) -- a *duration*, not a sample count, so
    gaps in the native sounding cadence (missing soundings, station
    outages) don't distort the average the way a fixed-N rolling window
    would. Works for any numeric parameter, not just foF2/hmF2.

    Returns an empty/NaN-dropped series indexed by ``datetime`` when
    ``window`` is falsy (smoothing disabled) or there is no valid data.
    """
    sub = df.dropna(subset=[param]).set_index("datetime").sort_index()
    if not window or sub.empty:
        return sub[param].iloc[0:0]
    return sub[param].rolling(window, min_periods=1).mean()


def _date_range(start: pd.Timestamp, n_days: int) -> list[pd.Timestamp]:
    start = pd.Timestamp(start).normalize()
    return [start + timedelta(days=i) for i in range(n_days)]


def _csv_path(day: pd.Timestamp, station: str) -> Path:
    return DATA_DIR / str(day.year) / station / f"{day.strftime('%Y-%m-%d')}.csv"


def load_range(start: pd.Timestamp, n_days: int, station: str = STATION, param: str = PARAM) -> pd.DataFrame:
    """Load and concatenate every available daily CSV over the date range."""
    frames = []
    for day in _date_range(start, n_days):
        path = _csv_path(day, station)
        if not path.exists():
            continue
        df = pd.read_csv(path, parse_dates=["datetime"])
        if param not in df.columns:
            raise ValueError(f"{path} has no {param!r} column; available: {sorted(df.columns)}")
        if df.empty:
            continue  # header-only "empty" day; an empty frame here would still
            # trip pandas' all-NA-column FutureWarning on concat below.
        day_frame = df[["datetime", param]].copy()
        day_frame[param] = pd.to_numeric(day_frame[param], errors="coerce")
        frames.append(day_frame)
    if not frames:
        logger.warning(f"No CSVs found for {station} in [{start.date()}, +{n_days}d)")
        return pd.DataFrame(columns=["datetime", param])
    combined = pd.concat(frames, ignore_index=True)
    combined = combined.sort_values("datetime").reset_index(drop=True)
    logger.info(f"Loaded {len(combined)} rows for {station}, {combined[param].notna().sum()} valid {param}")
    return combined


def plot_ts(df: pd.DataFrame, out_path: Path | None = None, station: str = STATION,
            param: str = PARAM, font_size: int = 9, figsize: tuple = (10, 3.5),
            smooth_window: str = SMOOTH_WINDOW) -> plt.Figure:
    """Plot the full native-cadence time series for one station/parameter,
    with an optional smoothed overlay (see ``smooth_series``)."""
    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

    sub = df.dropna(subset=[param])
    ax.plot(sub["datetime"], sub[param], color="tab:blue", lw=0.5, alpha=0.35,
             marker=".", ms=1.5, label="raw", zorder=2)

    smoothed = smooth_series(df, param, smooth_window)
    if not smoothed.empty:
        ax.plot(smoothed.index, smoothed.values, color="tab:red", lw=1.2, alpha=0.9,
                 label=f"{smooth_window} mean", zorder=3)

    ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m-%d"))
    ax.set_xlabel("UT", fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size - 1)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    ax.grid(axis="x", which="major", ls=":", lw=0.4, alpha=0.5)
    ax.grid(axis="y", which="major", ls=":", lw=0.4, alpha=0.5)
    units = _units_for(param)
    ax.set_ylabel(f"{param}" + (f" [{units}]" if units else ""), fontsize=font_size)
    ax.set_title(f"{station} - {param}", fontsize=font_size, loc="left")
    if not smoothed.empty:
        ax.legend(fontsize=font_size - 2, loc="best", framealpha=0.75, edgecolor="none")

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved: {out_path}")

    return fig


def discover_stations(data_dir: Path = DATA_DIR) -> list[str]:
    """Every station with at least one downloaded year under ``data_dir``."""
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return []
    stations = set()
    for year_dir in data_dir.iterdir():
        if year_dir.is_dir() and year_dir.name.isdigit():
            stations.update(p.name for p in year_dir.iterdir() if p.is_dir())
    return sorted(stations)


def _all_requested_days() -> list[date]:
    """Every day in config.toml's [download] years/month range, capped at today."""
    cfg = utils.load_config()["download"]
    today = date.today()
    days: list[date] = []
    for year in sorted(cfg["years"]):
        start = date(year, cfg["month_start"], 1)
        month_end = cfg["month_end"]
        last = (date(year, month_end + 1, 1) - timedelta(days=1)) if month_end < 12 else date(year, 12, 31)
        day = start
        while day <= min(last, today):
            days.append(day)
            day += timedelta(days=1)
    return days


def _downloaded_day_csvs(data_dir: Path = DATA_DIR):
    """Yield (station, day, csv_path) for every downloaded per-day CSV.

    Discovers directly from what's on disk (not config.toml's season window),
    so it also picks up ``--raw``-only or partially-downloaded stations/years.
    """
    data_dir = Path(data_dir)
    if not data_dir.is_dir():
        return
    for year_dir in sorted(data_dir.iterdir()):
        if not (year_dir.is_dir() and year_dir.name.isdigit()):
            continue
        for station_dir in sorted(year_dir.iterdir()):
            if not station_dir.is_dir():
                continue
            for csv_path in sorted(station_dir.glob("????-??-??.csv")):
                try:
                    day = date.fromisoformat(csv_path.stem)
                except ValueError:
                    continue
                yield station_dir.name, day, csv_path


def plot_day_ts(df: pd.DataFrame, day: date, out_path: Path | None = None,
                 station: str = STATION, param: str = PARAM,
                 font_size: int = 9, figsize: tuple = (6, 3.5),
                 smooth_window: str = SMOOTH_WINDOW) -> plt.Figure:
    """Plot one station's one-day time series, UT-of-day on the x-axis,
    with an optional smoothed overlay (see ``smooth_series``)."""
    fig, ax = plt.subplots(1, 1, figsize=figsize, constrained_layout=True)

    sub = df.dropna(subset=[param])
    ax.plot(sub["datetime"], sub[param], color="tab:blue", lw=0.6, alpha=0.4,
             marker=".", ms=3, label="raw", zorder=2)

    smoothed = smooth_series(df, param, smooth_window)
    if not smoothed.empty:
        ax.plot(smoothed.index, smoothed.values, color="tab:red", lw=1.6, alpha=0.9,
                 label=f"{smooth_window} mean", zorder=3)

    day_start = pd.Timestamp(day)
    ax.set_xlim(day_start, day_start + pd.Timedelta(days=1))
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax.set_xlabel("UT", fontsize=font_size)
    ax.tick_params(axis="both", labelsize=font_size - 1)
    ax.grid(axis="x", which="major", ls=":", lw=0.4, alpha=0.5)
    ax.grid(axis="y", which="major", ls=":", lw=0.4, alpha=0.5)
    units = _units_for(param)
    ax.set_ylabel(f"{param}" + (f" [{units}]" if units else ""), fontsize=font_size)
    ax.set_title(f"{station} - {param} - {day.isoformat()}", fontsize=font_size, loc="left")
    if not smoothed.empty:
        ax.legend(fontsize=font_size - 2, loc="best", framealpha=0.75, edgecolor="none")

    if out_path is not None:
        out_path = Path(out_path)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved: {out_path}")

    return fig


def plot_all_daily(params: list[str] = ALL_PARAMS, output_dir: Path = ALL_OUTPUT_DIR,
                    data_dir: Path = DATA_DIR) -> None:
    """One figure per (station) x (day) x (param), for every downloaded day."""
    saved, skipped = 0, 0
    for station, day, csv_path in _downloaded_day_csvs(data_dir):
        raw = pd.read_csv(csv_path, parse_dates=["datetime"])
        for param in params:
            if param not in raw.columns:
                logger.warning(f"{csv_path} has no {param!r} column; skipping")
                skipped += 1
                continue
            sub = raw[["datetime", param]].copy()
            sub[param] = pd.to_numeric(sub[param], errors="coerce")
            if sub[param].notna().sum() == 0:
                skipped += 1
                continue
            out = output_dir / station / f"{param}_ts_{station}_{day.strftime('%Y%m%d')}.png"
            fig = plot_day_ts(sub, day, out_path=out, station=station, param=param)
            plt.close(fig)
            saved += 1
    logger.info(f"plot_all_daily: saved {saved} figures, skipped {skipped} (no data) -> {output_dir}")


def plot_all(params: list[str] = ALL_PARAMS, output_dir: Path = ALL_OUTPUT_DIR,
             data_dir: Path = DATA_DIR) -> None:
    """One figure per (discovered station) x (param), full downloaded range."""
    days = _all_requested_days()
    if not days:
        logger.warning("No requested dates in config.toml [download]; nothing to plot.")
        return
    start, end = pd.Timestamp(days[0]), pd.Timestamp(days[-1])
    n_days = (end - start).days + 1

    stations = discover_stations(data_dir)
    if not stations:
        logger.warning(f"No stations found under {data_dir}; run didbase_pipeline.py first.")
        return

    date0, date1 = start.strftime("%Y%m%d"), end.strftime("%Y%m%d")
    saved, skipped = 0, 0
    for station in stations:
        for param in params:
            df = load_range(start, n_days, station, param)
            if df.empty or df[param].notna().sum() == 0:
                logger.warning(f"{station} {param}: no data in range, skipping plot")
                skipped += 1
                continue
            out = output_dir / f"{param}_ts_{station}_{date0}_{date1}.png"
            fig = plot_ts(df, out_path=out, station=station, param=param)
            plt.close(fig)
            saved += 1
    logger.info(f"plot_all: saved {saved} figures, skipped {skipped} (no data) -> {output_dir}")


if __name__ == "__main__":
    if RUN_DAILY:
        plot_all_daily()
    elif RUN_ALL:
        plot_all()
    else:
        df = load_range(START_DATE, N_DAYS, STATION, PARAM)
        if df.empty or df[PARAM].notna().sum() == 0:
            raise SystemExit("No data loaded; nothing to plot.")

        date0 = START_DATE.strftime("%Y%m%d")
        date1 = (START_DATE + timedelta(days=N_DAYS - 1)).strftime("%Y%m%d")
        fig = plot_ts(df, out_path=OUTPUT_DIR / f"{PARAM}_ts_{STATION}_{date0}_{date1}.png")
        plt.close(fig)
