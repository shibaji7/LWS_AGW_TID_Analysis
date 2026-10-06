"""Quality control for downloaded DIDBASE scaled parameters, plus the TID band-pass.

Per (station, day, parameter), from the daily CSVs written by
``didbase_pipeline.py``:

1. Sample level: drop NaN, values outside a physical range, and single-point
   spikes (Hampel filter: deviation from the rolling median > k scaled MADs).
2. Segment level: split valid samples wherever the time gap exceeds
   ``max(max_gap_min, 2.5 * cadence)`` and keep the *longest* segment.
3. Day level: the day passes for that parameter only if the kept segment is
   long enough, densely sampled, has few spikes, and its sample-to-sample noise
   is not much above that station's own median noise.

``python quality.py`` scans ``Data/didbase_giro_scaled`` and writes
``Data/quality/day_quality.csv`` (one row per station/day/parameter, with the
reason for every failure). Raw downloads are never modified.

``bandpass`` is a zero-phase Butterworth filter over the kept segment, used by
the plots. Thresholds all live in config.toml ``[quality]`` / ``[bandpass]``.
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger
from scipy.signal import butter, sosfiltfilt

import utils

ROOT = Path(__file__).resolve().parent
DATA_DIR = Path(os.environ.get("QC_DATA_DIR", str(ROOT / "Data" / "didbase_giro_scaled")))
QC_PATH = Path(os.environ.get("QC_OUTPUT", str(ROOT / "Data" / "quality" / "day_quality.csv")))


def cfg() -> dict:
    return utils.load_config()["quality"]


def load_day(csv_path: Path, params) -> pd.DataFrame:
    """Read one daily CSV: datetime + params, numeric, duplicate timestamps averaged."""
    df = pd.read_csv(csv_path, usecols=lambda c: c == "datetime" or c in params, parse_dates=["datetime"])
    for p in params:
        df[p] = pd.to_numeric(df[p], errors="coerce") if p in df else np.nan
    return df.groupby("datetime", as_index=False).mean().sort_values("datetime").reset_index(drop=True)


def cadence_min(times: pd.Series) -> float:
    """Median sounding interval [min], rounded to 0.5 min; NaN if undeterminable."""
    d = times.sort_values().diff().dt.total_seconds().dropna() / 60.0
    d = d[d > 0]
    return round(float(d.median()) * 2) / 2 if len(d) else float("nan")


def hampel_mask(values: np.ndarray, window: int, k: float, min_scale: float) -> np.ndarray:
    """True where a sample is a spike relative to its rolling median.

    ``min_scale`` floors the MAD so quantised/flat runs don't flag every wiggle.
    """
    s = pd.Series(values)
    med = s.rolling(window, center=True, min_periods=1).median()
    mad = (s - med).abs().rolling(window, center=True, min_periods=1).median() * 1.4826
    return ((s - med).abs() > k * np.maximum(mad, min_scale)).to_numpy()


def split_segments(times: pd.Series, max_gap_min: float) -> list[np.ndarray]:
    """Index arrays (into ``times``) of runs with no gap above ``max_gap_min``."""
    if len(times) == 0:
        return []
    gaps = times.diff().dt.total_seconds().to_numpy() / 60.0
    cuts = np.flatnonzero(gaps > max_gap_min)
    return np.split(np.arange(len(times)), cuts)


def noise_level(values: np.ndarray) -> float:
    """Robust sample-to-sample noise: MAD of the second difference / sqrt(6).

    A second difference removes slow trends *and* smooth TID-scale waves, so
    real wave activity is not counted as noise.
    """
    if len(values) < 5:
        return float("nan")
    d2 = np.diff(values, 2)
    return float(1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6.0))


def _units(param: str) -> str:
    return "km" if param.lower().startswith("h") else "MHz"


def clean_day(df: pd.DataFrame, param: str, c: dict | None = None):
    """Clean one day of one parameter.

    Returns ``(segment, metrics)``: ``segment`` is the longest clean segment
    (columns ``datetime``, ``param``; empty if none), ``metrics`` a dict of QC
    numbers. ``metrics['reason']`` is '' until ``evaluate`` fills it in.
    """
    c = c or cfg()
    lo, hi = c["range"][param]
    cad = cadence_min(df["datetime"])
    m = dict(cadence_min=cad, n_rows=len(df), n_valid=0, n_spikes=0, spike_frac=np.nan,
             seg_start=pd.NaT, seg_hours=0.0, seg_n=0, seg_coverage=0.0, noise=np.nan)
    empty = pd.DataFrame(columns=["datetime", param])
    if not np.isfinite(cad):
        return empty, m
    v = df[["datetime", param]].dropna()
    v = v[v[param].between(lo, hi)].reset_index(drop=True)
    m["n_valid"] = len(v)
    if len(v) < 3:
        return empty, m
    window = max(3, int(round(c["hampel_window_min"] / cad)) | 1)
    spike = hampel_mask(v[param].to_numpy(), window, c["hampel_k"], c["min_scale"][_units(param)])
    m["n_spikes"], m["spike_frac"] = int(spike.sum()), float(spike.mean())
    v = v[~spike].reset_index(drop=True)
    segs = split_segments(v["datetime"], max(c["max_gap_min"], 2.5 * cad))
    if not segs:
        return empty, m
    best = max(segs, key=len)
    seg = v.iloc[best].reset_index(drop=True)
    span_h = (seg["datetime"].iloc[-1] - seg["datetime"].iloc[0]).total_seconds() / 3600.0
    m.update(seg_start=seg["datetime"].iloc[0], seg_hours=span_h, seg_n=len(seg),
             seg_coverage=min(1.0, len(seg) / (span_h * 60.0 / cad + 1)),
             noise=noise_level(seg[param].to_numpy()))
    return seg, m


def evaluate(metrics: dict, station_noise: float | None, c: dict | None = None) -> str:
    """Return '' if the day passes, else a ';'-joined list of failure reasons."""
    c = c or cfg()
    why = []
    if metrics["seg_hours"] < c["min_segment_h"]:
        why.append(f"segment<{c['min_segment_h']}h")
    if metrics["seg_coverage"] < c["min_coverage"]:
        why.append(f"coverage<{c['min_coverage']}")
    if np.isfinite(metrics["spike_frac"]) and metrics["spike_frac"] > c["max_spike_frac"]:
        why.append(f"spikes>{c['max_spike_frac']}")
    noise = metrics["noise"]
    if station_noise and np.isfinite(noise) and noise > c["max_noise_rel"] * station_noise:
        why.append(f"noise>{c['max_noise_rel']}x station median")
    return ";".join(why)


def day_csvs(data_dir: Path = DATA_DIR):
    """Yield (station, date, path) for every non-empty daily CSV."""
    for path in sorted(Path(data_dir).glob("[0-9][0-9][0-9][0-9]/*/????-??-??.csv")):
        if path.stat().st_size > 0:
            yield path.parent.name, pd.Timestamp(path.stem), path


def build_quality_table(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """QC every station/day/parameter and return the table (see module docstring)."""
    c = cfg()
    params = list(c["parameters"])
    rows = []
    for station, day, path in day_csvs(data_dir):
        df = load_day(path, params)
        if df.empty:
            continue
        for p in params:
            _, m = clean_day(df, p, c)
            rows.append(dict(station=station, date=day.date().isoformat(), param=p, **m))
    table = pd.DataFrame(rows)
    if table.empty:
        return table
    # Noise reference per station/parameter: median over days that otherwise have a usable segment.
    usable = table[(table.seg_hours >= c["min_segment_h"]) & table.noise.notna()]
    ref = usable.groupby(["station", "param"]).noise.median()
    table["station_noise"] = [ref.get((s, p), np.nan) for s, p in zip(table.station, table.param)]
    table["reason"] = [evaluate(r, r["station_noise"]) for r in table.to_dict("records")]
    table["passed"] = table.reason.eq("")
    return table


def load_quality(path: Path = QC_PATH) -> pd.DataFrame:
    return pd.read_csv(path, parse_dates=["seg_start"])


# ---------------------------------------------------------------------------
# Band-pass for TID-scale variability
# ---------------------------------------------------------------------------

def effective_band(cadence: float, c: dict | None = None) -> tuple[float, float]:
    """(short, long) period limits [min]. The short limit is raised to
    ``nyquist_margin * cadence`` so the band stays below the Nyquist period
    (e.g. 15-min cadence cannot resolve 15-30 min periods)."""
    b = (c or utils.load_config())["bandpass"]
    return max(b["period_min_short"], b["nyquist_margin"] * cadence), b["period_min_long"]


def bandpass(seg: pd.DataFrame, param: str, cadence: float, c: dict | None = None):
    """Zero-phase band-pass of one clean segment.

    Linearly resamples onto a uniform ``cadence`` grid (gaps in a clean segment
    are at most 2.5 samples), then applies a Butterworth band-pass over
    ``effective_band``. Returns ``(Series indexed by time, (short, long))`` or
    ``(None, band)`` when the segment is too short for the filter to be valid.
    Output is NaN inside gaps wider than ``bandpass.mask_gap_min`` (default 30 min).
    """
    c = c or utils.load_config()
    b = c["bandpass"]
    band = effective_band(cadence, c)
    if len(seg) < 4 or band[0] >= band[1]:
        return None, band
    s = seg.set_index("datetime")[param]
    grid = pd.date_range(s.index[0].ceil(f"{cadence}min"), s.index[-1], freq=f"{cadence}min")
    # sosfiltfilt needs len > padlen (3*(2*sections+1)); also want >= ~3 long periods.
    if len(grid) <= max(3 * (2 * b["order"] + 1) + 1, 3 * band[1] / cadence):
        return None, band
    t = s.index.astype("int64").to_numpy()
    y = np.interp(grid.astype("int64"), t, s.to_numpy())
    fs = 1.0 / cadence  # cycles per minute
    sos = butter(b["order"], [1.0 / band[1], 1.0 / band[0]], btype="band", fs=fs, output="sos")
    out = pd.Series(sosfiltfilt(sos, y - y.mean()), index=grid)
    # Mask output where it rests on interpolation: grid points lying inside a gap
    # between real samples wider than mask_gap_min carry no real TID-band information.
    mask_ns = b.get("mask_gap_min", 30.0) * 60e9
    g = grid.astype("int64").to_numpy()
    nxt = np.searchsorted(t, g, side="right")  # first sample after each grid point
    nxt_c, prv_c = np.minimum(nxt, len(t) - 1), np.maximum(nxt - 1, 0)
    out[((t[nxt_c] - t[prv_c]) > mask_ns) & (g > t[prv_c])] = np.nan  # strictly inside the gap
    return out, band


if __name__ == "__main__":
    table = build_quality_table()
    QC_PATH.parent.mkdir(parents=True, exist_ok=True)
    table.to_csv(QC_PATH, index=False)
    logger.info(f"Wrote {len(table)} rows -> {QC_PATH}")
    if not table.empty:
        print(table.groupby(["param"]).passed.agg(["sum", "count", "mean"]).round(3))
        print(table.groupby(["station", "param"]).passed.mean().unstack().round(2))
        fails = table[~table.passed].reason.str.split(";").explode().str.replace(r"[<>][\d.]+.*", "", regex=True)
        print("\nFailure reasons:\n", fails.value_counts())
