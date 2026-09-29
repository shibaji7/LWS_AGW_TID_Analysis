"""Shared SZA regimes and station/year RSF scaling overrides.

All numeric settings are read from ``config.toml`` (same directory) so
tuning changes do not require editing code. This module is imported by
``rsf_batch_scaler.py`` and ``plot_viterbi.py``.
"""
from __future__ import annotations

import datetime as dt
import math
import tomllib
from functools import lru_cache
from pathlib import Path

import pandas as pd

CONFIG_PATH = Path(__file__).resolve().parent / "config.toml"


@lru_cache(maxsize=1)
def load_config() -> dict:
    with open(CONFIG_PATH, "rb") as handle:
        return tomllib.load(handle)


# ---------------------------------------------------------------------------
# Solar zenith angle + SZA-based parameter sets
# ---------------------------------------------------------------------------

def _sza_regimes_config() -> dict:
    return load_config()["sza_regimes"]


def _solar_zenith_angle(datetime_utc: dt.datetime, lat_deg: float, lon_deg: float) -> float:
    """Return approximate solar zenith angle [degrees].

    Uses the standard declination/hour-angle formula (Spencer 1971), accurate
    to +/-1 deg -- sufficient for day/twilight/night classification.

    Parameters
    ----------
    datetime_utc : UTC datetime of the sounding.
    lat_deg      : Station geodetic latitude  [deg N].
    lon_deg      : Station longitude [deg E]; negative for West.

    Returns
    -------
    float : SZA in degrees (0 deg = overhead, 90 deg = horizon, >90 deg = below horizon).
    """
    doy = datetime_utc.timetuple().tm_yday
    ut_h = datetime_utc.hour + datetime_utc.minute / 60.0 + datetime_utc.second / 3600.0
    decl = math.radians(23.45 * math.sin(math.radians(360.0 / 365.0 * (doy - 81))))
    ha = math.radians(15.0 * (ut_h - 12.0) + lon_deg)
    lat_r = math.radians(lat_deg)
    cos_sza = (math.sin(lat_r) * math.sin(decl) + math.cos(lat_r) * math.cos(decl) * math.cos(ha))
    return math.degrees(math.acos(max(-1.0, min(1.0, cos_sza))))


def _regime_for_sza(sza: float) -> str:
    """Map SZA [degrees] -> 'day' | 'twilight' | 'night'."""
    cfg = _sza_regimes_config()
    if sza < cfg["day_max_deg"]:
        return "day"
    elif sza < cfg["twilight_max_deg"]:
        return "twilight"
    return "night"


class _SzaParamSets:
    """Lazily built, dict-like view of config.toml's [sza_regimes.<regime>] tables."""

    def __getitem__(self, regime: str) -> dict:
        cfg = _sza_regimes_config()
        return {k: v for k, v in cfg[regime].items()}


_SZA_PARAM_SETS = _SzaParamSets()


# ---------------------------------------------------------------------------
# Per-station, per-year RSF scaler overrides (config.toml [rsf_scaler.*])
# ---------------------------------------------------------------------------

def rsf_scaler_config() -> dict:
    return load_config()["rsf_scaler"]


def station_year_overrides(station: str, year: int) -> dict:
    """Return the override dict for one station/year, or {} if unlisted.

    Unlike the retired eclipse workflow (which required every station/year
    combination to be listed explicitly), this now covers all CONUS stations
    across three seasons, so an unlisted combination silently falls back to
    the plain SZA-regime defaults rather than raising.
    """
    overrides = rsf_scaler_config().get("station_overrides", {})
    return overrides.get(station, {}).get(str(year), {})


def adaptive_min_amp_db(normalized: pd.DataFrame, margin_db: float | None = None) -> float | None:
    """Data-driven amplitude floor: this file's own noise-floor median + margin."""
    if margin_db is None:
        margin_db = rsf_scaler_config()["adaptive_amp_margin_db"]
    amp = normalized.get("amplitude_db")
    if amp is None or amp.empty:
        return None
    return float(amp.median()) + margin_db


def station_amplitude_threshold(normalized: pd.DataFrame, station: str, year: int, layer: str) -> float | None:
    threshold = station_year_overrides(station, year).get("prefilter_min_amp_dB", {}).get(layer)
    if threshold == "adaptive":
        return adaptive_min_amp_db(normalized)
    return threshold


def station_amplitude_prefilter(
    raw: pd.DataFrame,
    normalized: pd.DataFrame,
    threshold: float | None,
    masked_bands: list[tuple[float, float]] | None = None,
) -> pd.DataFrame:
    """Drop raw echoes below a resolved amplitude floor and inside masked bands."""
    if threshold is None and not masked_bands:
        return raw
    keep = pd.Series(True, index=normalized.index)
    if threshold is not None:
        keep &= normalized["amplitude_db"] > threshold
    for lo, hi in masked_bands or []:
        keep &= ~normalized["frequency_mhz"].between(lo, hi)
    return raw[keep.to_numpy()]


def trim_trace_by_gap(trace: pd.DataFrame, max_gap_mhz: float, bridge: bool = False) -> pd.DataFrame:
    """Keep the largest frequency-contiguous segment of a Viterbi trace.

    A fixed-height multipath/artifact echo (or a stray point) can sit apart
    from the real trace with a frequency gap wider than ``max_gap_mhz``, so
    split on every gap and keep the largest segment. With ``bridge=True``,
    neighbouring segments that are both at least
    ``rsf_scaler.bridge_min_segment_points`` long and whose heights differ by
    at most ``rsf_scaler.bridge_max_height_step_km`` across the gap are
    merged first, so a real trace with an internal dropout isn't reduced to
    its longer half.
    """
    if trace is None or trace.empty or len(trace) < 2:
        return trace
    cfg = rsf_scaler_config()
    trace = trace.sort_values("frequency_mhz").reset_index(drop=True)
    gaps = trace["frequency_mhz"].diff()
    cut_positions = gaps[gaps > max_gap_mhz].index.tolist()
    if not cut_positions:
        return trace
    bounds = [0] + cut_positions + [len(trace)]
    segments = [trace.iloc[bounds[i]:bounds[i + 1]] for i in range(len(bounds) - 1)]
    if bridge:
        min_points = cfg["bridge_min_segment_points"]
        max_step = cfg["bridge_max_height_step_km"]
        merged = [segments[0]]
        for seg in segments[1:]:
            prev = merged[-1]
            step = abs(float(seg["height_km"].iloc[0]) - float(prev["height_km"].iloc[-1]))
            if len(prev) >= min_points and len(seg) >= min_points and step <= max_step:
                merged[-1] = pd.concat([prev, seg])
            else:
                merged.append(seg)
        segments = merged
    return max(segments, key=len)


def trim_trailing_gap_contamination(
    row: dict, details: dict, layer: str, fo_key: str, hprime_key: str,
    max_gap_mhz: float, bridge: bool = False,
) -> None:
    trace = details.get("traces", {}).get(layer)
    if trace is None or trace.empty:
        return
    trimmed = trim_trace_by_gap(trace, max_gap_mhz, bridge)
    if len(trimmed) == len(trace):
        return
    if trimmed.empty:
        row[fo_key] = float("nan")
        row[hprime_key] = float("nan")
        row[f"{layer}_trace_n_points"] = 0
    else:
        last = trimmed.iloc[-1]
        row[fo_key] = float(last["frequency_mhz"])
        row[hprime_key] = float(last["height_km"])
        row[f"{layer}_trace_n_points"] = len(trimmed)
        row[f"{layer}_trace_freq_max_mhz"] = float(last["frequency_mhz"])


def apply_reliability_gate(row: dict, layer: str, fo_key: str, hprime_key: str, flag_key: str) -> None:
    """Null out a layer's foX/h'X when its trace reliability is too low.

    The Viterbi tracer can lock onto a handful of noise points and still
    report a foX value; ``{layer}_trace_reliability_score`` is a much better
    quality signal than "a trace was found at all", so gate on it instead of
    trusting every detection.
    """
    cfg = rsf_scaler_config()
    score = row.get(f"{layer}_trace_reliability_score")
    npts = row.get(f"{layer}_trace_n_points")
    insufficient_points = npts is None or pd.isna(npts) or npts < cfg["min_trace_points"]
    if score is None or pd.isna(score) or score < cfg["min_trace_reliability"] or insufficient_points:
        row[fo_key] = float("nan")
        row[hprime_key] = float("nan")
        row[flag_key] = False
