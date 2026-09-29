"""Viterbi ionogram diagnostics for an arbitrary UT window.

Reruns the same scaling logic as ``rsf_batch_scaler.py`` (SZA regime + any
station/year override from ``config.toml``) for selected raw RSF files and
saves raw/cleaned/Viterbi overlay figures. Meant for diagnosing suspicious
foE/foF2 picks without changing the main scaled CSV outputs.

Raw files come from ``didbase_pipeline.py --raw``:
    Data/didbase_giro_scaled/<year>/<station>/raw/<YYYY-MM-DD>/*.RSF

Run from this directory, for example:

    RSF_VITERBI_STATION=AU930 \
    RSF_VITERBI_START_DATE=2024-04-08 \
    RSF_VITERBI_N_DAYS=1 \
    RSF_VITERBI_UT_START=15 RSF_VITERBI_UT_END=20 \
    conda run -n apep python plot_viterbi.py
"""

from __future__ import annotations

import os
from datetime import timedelta
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from loguru import logger

import rsf_batch_scaler as scaler
import utils

_VITERBI_CFG = utils.load_config()["viterbi_plot"]

STATION = os.environ.get("RSF_VITERBI_STATION", scaler.STATION)
DATA_ROOT = Path(os.environ.get("RSF_VITERBI_DATA_ROOT", str(scaler.DATA_ROOT)))
START_DATE = pd.Timestamp(os.environ.get("RSF_VITERBI_START_DATE", "2024-04-08"))
N_DAYS = int(os.environ.get("RSF_VITERBI_N_DAYS", "1"))
UT_START = float(os.environ.get("RSF_VITERBI_UT_START", _VITERBI_CFG["ut_start"]))
UT_END = float(os.environ.get("RSF_VITERBI_UT_END", _VITERBI_CFG["ut_end"]))
KEEP_EVERY = int(os.environ.get("RSF_VITERBI_KEEP_EVERY", _VITERBI_CFG["keep_every"]))
OUTPUT_DIR = Path(os.environ.get("RSF_VITERBI_OUTPUT_DIR", str(Path(__file__).resolve().parent / _VITERBI_CFG["output_dir"])))


def _date_range(start: pd.Timestamp, n_days: int) -> list[pd.Timestamp]:
    start = pd.Timestamp(start).normalize()
    return [start + timedelta(days=i) for i in range(n_days)]


def _files_for_day(day: pd.Timestamp, station: str = STATION) -> list[Path]:
    date_dir = DATA_ROOT / str(day.year) / station / "raw" / day.strftime("%Y-%m-%d")
    return sorted(date_dir.glob("*.RSF"))


def _decimal_ut(t: pd.Timestamp) -> float:
    return t.hour + t.minute / 60.0 + t.second / 3600.0


def _selected_files_for_day(day: pd.Timestamp, station: str = STATION) -> list[Path]:
    """Files within [UT_START, UT_END) UT, subsampled every KEEP_EVERY-th."""
    from pynasonde.digisonde.parsers.rsf import RsfExtractor

    selected: list[tuple[pd.Timestamp, Path]] = []
    for path in _files_for_day(day, station):
        try:
            ext = RsfExtractor(path, extract_time_from_name=True, extract_stn_from_name=True)
        except Exception as exc:
            logger.warning(f"Could not parse timestamp from {path.name}: {exc}")
            continue
        t = pd.Timestamp(ext.date)
        if UT_START <= _decimal_ut(t) < UT_END:
            selected.append((t, path))

    selected = sorted(selected, key=lambda item: item[0])
    if KEEP_EVERY > 1:
        selected = selected[::KEEP_EVERY]
    return [path for _, path in selected]


def _plot_one(result: dict, title: str, out_path: Path) -> None:
    raw_o = result["normalized"]
    raw_o = raw_o[raw_o["mode"] == "O"].copy() if "mode" in raw_o.columns else raw_o

    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)
    if not raw_o.empty:
        plot_raw = raw_o[raw_o["height_km"].between(70, 600) & raw_o["frequency_mhz"].between(1, 14)]
        ax.scatter(
            plot_raw["frequency_mhz"], plot_raw["height_km"], c=plot_raw["amplitude_db"],
            cmap="Greys", s=1.0, alpha=0.20, linewidths=0, rasterized=True, label="raw O-mode",
        )

    cleaned = result["cleaned"]
    if not cleaned.empty and {"frequency_mhz", "height_km"}.issubset(cleaned.columns):
        band = cleaned[cleaned["height_km"].between(70, 600) & cleaned["frequency_mhz"].between(1, 14)]
        ax.scatter(band["frequency_mhz"], band["height_km"], s=4.0, color="tab:orange",
                   alpha=0.22, linewidths=0, rasterized=True, label="cleaned echoes")

    for layer, color in [("E", "blue"), ("F2", "red")]:
        trace = result["traces"].get(layer, pd.DataFrame())
        if not trace.empty:
            trace = trace.sort_values("frequency_mhz")
            ax.plot(trace["frequency_mhz"], trace["height_km"], color=color, lw=2.0,
                    label=f"{layer} Viterbi ({len(trace)} pts)", zorder=5)

    row = result["row"]
    text = [
        f"SZA={row.get('sza_deg', float('nan')):.1f} deg, {row.get('sza_regime', '')}",
        f"E scaling={'on' if row.get('e_layer_enabled', False) else 'off'}",
        f"foE={row['foE_mhz']:.2f} MHz" if pd.notna(row.get("foE_mhz")) else "foE=NaN",
        f"foF2={row['foF2_mhz']:.2f} MHz" if pd.notna(row.get("foF2_mhz")) else "foF2=NaN",
    ]
    ax.text(0.02, 0.98, "\n".join(text), transform=ax.transAxes, va="top", ha="left",
            fontsize=9, bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.85))
    ax.set_title(title, loc="left", fontsize=11)
    ax.set_xlabel("Frequency [MHz]")
    ax.set_ylabel("Virtual height [km]")
    ax.set_xlim(1, 14)
    ax.set_ylim(70, 600)
    ax.grid(True, ls=":", lw=0.45, alpha=0.55)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.75, edgecolor="none")

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


def process_file(path: Path, station: str) -> None:
    result = scaler._process_one_file(str(path))
    if result is None:
        logger.warning(f"Skipping unscalable RSF: {path}")
        return
    t = pd.Timestamp(result["row"]["time"])
    date_dir = OUTPUT_DIR / station / t.strftime("%Y%m%d")
    title = f"{station} {t:%Y-%m-%d %H:%M:%S} UT"
    out_path = date_dir / f"viterbi_{station}_{t:%Y%m%d_%H%M%S}.png"
    _plot_one(result, title, out_path)


if __name__ == "__main__":
    total = 0
    for day in _date_range(START_DATE, N_DAYS):
        files = _selected_files_for_day(day, STATION)
        logger.info(
            f"{STATION} {day.strftime('%Y-%m-%d')}: selected {len(files)} RSF files "
            f"from {UT_START:g}-{UT_END:g} UT with KEEP_EVERY={KEEP_EVERY}"
        )
        for path in files:
            process_file(path, STATION)
            total += 1
    logger.info(f"Completed {total} Viterbi diagnostic plots")
