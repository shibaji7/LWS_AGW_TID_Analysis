"""Diagnostic scaled-parameter versus solar zenith angle plots.

Reads daily DIDBASE CSVs downloaded by ``didbase_pipeline.py``
(``Data/didbase_giro_scaled/<year>/<station>/<YYYY-MM-DD>.csv``) and plots
any one scaled parameter (foE, foF2, hF, ...) against SZA. Any parameter
present in the CSV can be selected via ``RSF_PARAM_SZA_PARAM`` -- not just foE
-- to diagnose whether questionable picks cluster near the day/twilight/night
transition.

Run from this directory, for example:

    RSF_PARAM_SZA_STATION=AU930 \
    RSF_PARAM_SZA_PARAM=foF2 \
    RSF_PARAM_SZA_START_DATE=2024-04-01 \
    RSF_PARAM_SZA_N_DAYS=3 \
    conda run -n apep python plot_param_vs_sza.py
"""

from __future__ import annotations

import os
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

import utils

_CFG = utils.load_config()["foe_sza_plot"]

STATION = os.environ.get("RSF_PARAM_SZA_STATION", "AU930")
DATA_DIR = Path(os.environ.get("RSF_PARAM_SZA_DATA_DIR", str(Path(__file__).resolve().parent / "Data" / "didbase_giro_scaled")))
OUTPUT_DIR = Path(os.environ.get("RSF_PARAM_SZA_OUTPUT_DIR", str(Path(__file__).resolve().parent / _CFG["output_dir"])))
START_DATE = pd.Timestamp(os.environ.get("RSF_PARAM_SZA_START_DATE", "2024-04-01"))
N_DAYS = int(os.environ.get("RSF_PARAM_SZA_N_DAYS", "3"))
# Any scaled parameter column present in the DIDBASE CSVs, e.g. foE, foF2, hF.
PARAM = os.environ.get("RSF_PARAM_SZA_PARAM", _CFG["parameter"]).removesuffix("_mhz")

# Station coordinates for the SZA calculation; falls back to config.toml's
# extra-station metadata, else raises so a station isn't silently misplaced.
_STATION_LOCATIONS = {
    "KAFB": (35.04, -106.61), "KR835": (35.04, -106.61),
    "WSMR": (32.40, -106.40), "WS833": (32.40, -106.40),
    "AU930": (30.40, -97.70),
}


def _station_location(station: str) -> tuple[float, float]:
    if station in _STATION_LOCATIONS:
        return _STATION_LOCATIONS[station]
    for row in utils.load_config()["stations"].get("extra", []):
        if row["ursi"] == station:
            return row["lat"], row["long"]
    raise KeyError(f"No lat/lon known for station {station!r}; add it to _STATION_LOCATIONS "
                    "or config.toml's [[stations.extra]]")


def _date_range(start: pd.Timestamp, n_days: int) -> list[pd.Timestamp]:
    start = pd.Timestamp(start).normalize()
    return [start + timedelta(days=i) for i in range(n_days)]


def _csv_path(day: pd.Timestamp, station: str) -> Path:
    return DATA_DIR / str(day.year) / station / f"{day.strftime('%Y-%m-%d')}.csv"


def load_day(day: pd.Timestamp, station: str, param: str = PARAM) -> pd.DataFrame:
    path = _csv_path(day, station)
    if not path.exists():
        logger.warning(f"Missing scaled CSV: {path}")
        return pd.DataFrame()

    df = pd.read_csv(path, parse_dates=["datetime"])
    if param not in df.columns:
        raise ValueError(f"{path} has no {param!r} column; available: {sorted(df.columns)}")

    lat, lon = _station_location(station)
    df["sza_deg"] = df["datetime"].apply(lambda t: utils._solar_zenith_angle(t.to_pydatetime(), lat, lon))
    df[param] = pd.to_numeric(df[param], errors="coerce")
    df = df.sort_values("datetime").reset_index(drop=True)
    logger.info(f"{day.strftime('%Y-%m-%d')} {station}: rows={len(df)}, {param}={df[param].notna().sum()}")
    return df


def plot_param_vs_sza(daily: dict[str, pd.DataFrame], out_path: Path, param: str = PARAM, font_size: int = 10) -> plt.Figure:
    if not daily:
        raise ValueError("No daily data loaded")

    fig, axes = plt.subplots(1, len(daily), figsize=(4.0 * len(daily), 3.6), sharey=True, constrained_layout=True)
    if len(daily) == 1:
        axes = [axes]

    last_scatter = None
    for ax, (label, df) in zip(axes, daily.items()):
        sub = df.dropna(subset=["sza_deg", param]).copy()
        if not sub.empty:
            last_scatter = ax.scatter(
                sub["sza_deg"], sub[param],
                c=sub["datetime"].dt.hour + sub["datetime"].dt.minute / 60.0,
                cmap="viridis", s=14, alpha=0.85, edgecolors="none", zorder=2,
            )
        ax.set_title(label, fontsize=font_size, loc="left")
        ax.set_xlabel("Solar zenith angle [deg]", fontsize=font_size)
        ax.set_xlim(45, 135)
        ax.grid(True, ls=":", lw=0.5, alpha=0.55)
        ax.tick_params(labelsize=font_size - 1)

    axes[0].set_ylabel(f"{param} [MHz]", fontsize=font_size)
    if last_scatter is not None:
        cbar = fig.colorbar(last_scatter, ax=axes, shrink=0.82, pad=0.02)
        cbar.set_label("UT [hours]", fontsize=font_size)
        cbar.ax.tick_params(labelsize=font_size - 1)

    fig.suptitle(f"{STATION}: raw {param} versus SZA", fontsize=font_size + 1)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    logger.info(f"Saved: {out_path}")
    return fig


if __name__ == "__main__":
    daily: dict[str, pd.DataFrame] = {}
    for day in _date_range(START_DATE, N_DAYS):
        df = load_day(day, STATION)
        if not df.empty:
            daily[day.strftime("%Y-%m-%d")] = df

    if not daily:
        raise SystemExit("No data loaded; nothing to plot.")

    date0 = START_DATE.strftime("%Y%m%d")
    date1 = (START_DATE + timedelta(days=N_DAYS - 1)).strftime("%Y%m%d")
    out = OUTPUT_DIR / f"{PARAM}_vs_sza_{date0}_{date1}_{STATION}.png"
    fig = plot_param_vs_sza(daily, out)
    plt.close(fig)
