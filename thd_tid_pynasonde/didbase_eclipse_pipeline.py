"""Build DIDBASE reference medians and eclipse-day locally scaled products.

For each configured station and eclipse event this script:

1. Queries DIDBASE scaled characteristics for the full calendar month.
2. Excludes the eclipse date and writes a time-of-day median reference.
3. Downloads the eclipse-day raw ionograms, trying RSF and then MMM when
   ``--raw-format auto`` is used.
4. Parses and locally scales the eclipse-day raw files with pynasonde.

Outputs are written under ``Data/`` relative to this file:

    Data/<station_code>_<eclipse_year>_median.csv
    Data/<station_code>_<eclipse_year>_scaled.csv

The raw DIDBASE downloads are cached below ``Data/didbase_raw/``. DIDBASE
credentials are read from ``DIDBASE_USER`` and ``DIDBASE_PASSWORD``; the URL
and role follow the pynasonde defaults unless overridden by environment
variables.
"""

from __future__ import annotations

import argparse
import datetime as dt
import os
import sys
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import matplotlib.patheffects as pe
import matplotlib.pyplot as plt
import pandas as pd
from joblib import Parallel, delayed
from loguru import logger

sys.path.insert(0, str(Path(__file__).resolve().parent))

from pynasonde.digisonde.didbase import DidBaseConnector
from pynasonde.digisonde.parsers.mmm import ModMaxExtractor
from pynasonde.digisonde.parsers.rsf import RsfExtractor
from pynasonde.digisonde.scaler import normalize_echo_dataframe, scale_dataframe
from pynasonde.digisonde.scaler.rfi_filter import VerticalRfiFilter

import utils
from rsf_batch_scaler import (
    E_SCALER_KWARGS,
    F2_SCALER_KWARGS,
    _add_trace_metrics,
    _layer_from_sza_params,
    _normal_e_layer_enabled,
    _scaler_kwargs_for_layer,
)


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "Data"
RAW_DIR = DATA_DIR / "didbase_raw"

DIDBASE_URL = os.environ.get(
    "DIDBASE_URL", "jdbc:firebirdsql://didbase.giro.uml.edu/didb"
)
DIDBASE_ROLE = os.environ.get("DIDBASE_ROLE", "COMMON")
DIDBASE_USER = os.environ.get("DIDBASE_USER", "ERAU")
DIDBASE_PASSWORD = os.environ.get("DIDBASE_PASSWORD", "eraucar")

DIDBASE_WORKERS = int(os.environ.get("DIDBASE_WORKERS", "1"))
RSF_SCALE_JOBS = int(os.environ.get("RSF_SCALE_JOBS", "4"))
RSF_RFI_ENABLED = os.environ.get("RSF_SCALER_RFI_ENABLED", "1") == "1"
RSF_RFI_MAX_ECHOES_PER_BIN = int(
    os.environ.get("RSF_SCALER_RFI_MAX_ECHOES_PER_BIN", "100")
)
RSF_RFI_FREQ_BIN_MHZ = float(
    os.environ.get("RSF_SCALER_RFI_FREQ_BIN_MHZ", "0.05")
)
RSF_MAX_TRACE_FREQ_GAP_MHZ = float(
    os.environ.get("RSF_SCALER_MAX_TRACE_FREQ_GAP_MHZ", "1.0")
)
RSF_TUNING_PRESET = os.environ.get("RSF_SCALER_TUNING_PRESET", "balanced")
E_LAYER_MAX_SZA = float(os.environ.get("RSF_SCALER_E_LAYER_MAX_SZA", "100"))
RSF_MIN_TRACE_RELIABILITY = float(
    os.environ.get("RSF_SCALER_MIN_TRACE_RELIABILITY", "0.4")
)
# A trace with very few points can be trivially "self-consistent" and score
# a misleadingly high reliability (e.g. BC840 2023-10-14 17:35 UT: a 3-point
# fragment scored 0.906 while the real ~100-point trace was clearly visible
# in the raw data), so gate on point count too, not reliability alone.
RSF_MIN_TRACE_POINTS = int(os.environ.get("RSF_SCALER_MIN_TRACE_POINTS", "10"))
# Gap bridging (opt-in per station/year): two gap-separated trace segments are
# re-joined only when both are substantial and heights are continuous across
# the gap, so a stray fragment or a fixed-height artifact is never bridged.
# 8, not 10: AU930 2024-04-08 12:55/13:00 UT have real-trace halves of 9 and 8 points.
RSF_BRIDGE_MIN_SEGMENT_POINTS = int(os.environ.get("RSF_SCALER_BRIDGE_MIN_POINTS", "8"))
RSF_BRIDGE_MAX_HEIGHT_STEP_KM = float(os.environ.get("RSF_SCALER_BRIDGE_MAX_STEP_KM", "50"))
VITERBI_UT_START = float(os.environ.get("RSF_VITERBI_UT_START", "15"))
VITERBI_UT_END = float(os.environ.get("RSF_VITERBI_UT_END", "20"))
VITERBI_KEEP_EVERY = int(os.environ.get("RSF_VITERBI_KEEP_EVERY", "1"))
VITERBI_DIR = DATA_DIR / "viterbi_diagnostics"

EVENTS = {
    2023: dt.date(2023, 10, 14),
    2024: dt.date(2024, 4, 8),
}

STATIONS = {
    # "AL945": {"name": "Alpena", "lat": 45.07, "lon": -83.56},
    # "BC840": {"name": "Boulder", "lat": 40.00, "lon": -105.30},
    "AU930": {"name": "Austin", "lat": 30.40, "lon": -97.70},
}

# Per-station, per-year overrides for scaler behavior that doesn't generalize
# across stations or even across sessions for the same station. AL945's raw
# echo density is high enough that the default RFI_MAX_ECHOES_PER_BIN=100
# threshold hard-deletes the real dense ionospheric trace (hundreds of
# echoes/bin) while a sparse fixed-height (~557 km) multipath artifact
# survives untouched and gets traced instead (confirmed by manual inspection
# of AL945 raw echoes on 2023-10-14). A flat "no_default" would silently
# apply one year's tuning to another; every station/year combo must be
# listed explicitly, or a per-file KeyError is raised (see
# _station_year_overrides) rather than guessing.
RSF_ADAPTIVE_AMP_MARGIN_DB = float(
    os.environ.get("RSF_SCALER_ADAPTIVE_AMP_MARGIN_DB", "10.0")
)

_BC840_COMMON_OVERRIDES = {
    "rfi_enabled": True,
    "e_skip_sza_f_hi_clamp": True,
    "max_trace_freq_gap_mhz": 0.4,
}

_AU930_COMMON_OVERRIDES = {
    "rfi_enabled": True,
    "e_skip_sza_f_hi_clamp": True,
    # 2023-10-14 16:00 UT: real trace has a 0.4 MHz internal dropout (6.0->6.4)
    # that must NOT split it, while the RFI-stripe tail starts after a >=1.0 MHz
    # jump (9.6->10.6) that must. 0.4 drops the 4-6 MHz part; 1.0 keeps the tail.
    "max_trace_freq_gap_mhz": 0.8,
}

_AL945_COMMON_OVERRIDES = {
    "rfi_enabled": False,
    # rsf_batch_scaler's SZA-based E f_hi clamp (f_hi = min(f_hi, 6*cos(SZA)),
    # floored at f_lo+2.0) is non-monotonic: it allows up to 4.0 MHz around
    # SZA~74-87 deg but drops to ~2.0-3.5 MHz for SZA~48-70 deg -- narrower
    # right when the sun is higher and foE should be rising, not falling.
    # AL945's real E edge sits at ~3.5-3.9 MHz in that range (visually
    # confirmed), so skip the clamp and use the regime's base f_hi instead.
    "e_skip_sza_f_hi_clamp": True,
}

STATION_OVERRIDES: dict[str, dict[int, dict]] = {
    "AL945": {
        # Visually validated against Viterbi diagnostics for 2023-10-14: a
        # flat 50 dB floor cleanly isolates the real trace on this specific
        # day without truncating it early.
        2023: {
            **_AL945_COMMON_OVERRIDES,
            "prefilter_min_amp_dB": {"E": 50.0, "F2": 50.0},
        },
        # A flat 50 dB floor deleted most of the real (weaker) trace on
        # 2024-04-08 (peaks only 45-66 dB vs 2023's 78-81 dB), causing
        # isolated single-sample Viterbi dropouts through the afternoon.
        # "adaptive" -> this file's own noise-floor median + margin instead
        # of a number tuned on a different day; self-calibrates per session.
        2024: {
            **_AL945_COMMON_OVERRIDES,
            "prefilter_min_amp_dB": {"E": 40, "F2": 40},
        },
    },
    "BC840": {
        2023: {
            **_BC840_COMMON_OVERRIDES,
            "prefilter_min_amp_dB": {"E": 55.0, "F2": 55.0},
        },
        2024: {
            **_BC840_COMMON_OVERRIDES,
            "prefilter_min_amp_dB": {"E": 55.0, "F2": 55.0},
        },
    },
    "AU930": {
        2023: {
            **_AU930_COMMON_OVERRIDES,
            # AU930's echoes are much weaker than AL945/BC840 (median O-mode
            # amplitude ~18-19 dB vs AL945's ~30 dB); a flat 45 dB floor left
            # only a handful of candidate points -- too sparse for DBSCAN to
            # find a cluster (F2 kept only 18/287 rows). Adaptive (median+10)
            # recovers a clean ~50-point trace on 2023-10-14 15:00 UT.
            "prefilter_min_amp_dB": {"E": 20, "F2": 20},
            # Below ~4 MHz the ionogram is dominated by vertical RFI stripes
            # spanning the full height window (visually confirmed on
            # 2023-10-14 15:00 UT); the real F2 trace only starts at f=4 MHz.
            "f2_f_lo_mhz": 4.0,
            # Persistent full-height RFI stripes inside the F2 band (>=18 of 30
            # sampled 2023-10-14 files); their stray points tail the real trace
            # and inflate foF2 (e.g. 14:50 UT read 11.3 MHz vs neighbors 7.2).
            "masked_freq_bands_mhz": [(5.15, 5.25), (10.55, 10.75), (11.15, 11.35)],
            # Real trace can fade out for ~1 MHz mid-curve (14:45/14:55/15:15 UT,
            # 7.2 -> 8.2+ MHz) and then resume; rejoin the two halves.
            "bridge_trace_gaps": True,
        },
        # 2024-04-08: same MMM format and 18 dB median noise floor as 2023;
        # stripes in the F2 band measured from 30 sampled files (present in
        # >=50%): 5.2, 10.6 and 11.0-11.2 MHz (upper stripe sits ~0.2 MHz lower
        # than in 2023). Bridging left off until a run shows it is needed.
        2024: {
            **_AU930_COMMON_OVERRIDES,
            "prefilter_min_amp_dB": {"E": 20, "F2": 20},
            "f2_f_lo_mhz": 4.0,
            "masked_freq_bands_mhz": [(5.15, 5.25), (10.55, 10.65), (10.95, 11.25)],
            # 18:00 UT: a 19-pt noise fragment (4.1-6.7 MHz) outnumbered the real
            # 17-pt trace (8.3-9.9 MHz); 12:55/13:00 UT: real trace with a >1 MHz
            # dropout was cut to its shorter half. Rejoin height-continuous halves.
            "bridge_trace_gaps": True,
        },
    }
}


def _station_year_overrides(station: str, year: int) -> dict:
    by_year = STATION_OVERRIDES.get(station)
    if by_year is None:
        # Station has no recorded quirks at all yet -- fall back to defaults.
        return {}
    if year not in by_year:
        # Station is known to need per-year tuning; a missing year here means
        # nobody has checked this session's Viterbi diagnostics yet, so don't
        # silently inherit another year's tuning.
        raise KeyError(
            f"No STATION_OVERRIDES tuning recorded for {station} {year}; "
            "add one (even an empty {}) after checking Viterbi diagnostics "
            "for this session rather than inheriting another year's tuning."
        )
    return by_year[year]


def _station_rfi_enabled(station: str, year: int) -> bool:
    return _station_year_overrides(station, year).get("rfi_enabled", RSF_RFI_ENABLED)


def _station_bridge_gaps(station: str, year: int) -> bool:
    return bool(_station_year_overrides(station, year).get("bridge_trace_gaps", False))


def _station_max_trace_gap(station: str, year: int) -> float:
    return _station_year_overrides(station, year).get(
        "max_trace_freq_gap_mhz", RSF_MAX_TRACE_FREQ_GAP_MHZ
    )


def _adaptive_min_amp_db(
    normalized: pd.DataFrame, margin_db: float = RSF_ADAPTIVE_AMP_MARGIN_DB
) -> float | None:
    """Data-driven amplitude floor: this file's own noise-floor median + margin.

    Self-calibrates per session instead of assuming one day's SNR generalizes
    to another (confirmed: AL945 2023-10-14 and 2024-04-08 share the same
    30.0 dB median noise floor despite very different real-echo strength).
    """
    amp = normalized.get("amplitude_db")
    if amp is None or amp.empty:
        return None
    return float(amp.median()) + margin_db


def _station_amplitude_threshold(
    normalized: pd.DataFrame, station: str, year: int, layer: str
) -> float | None:
    threshold = (
        _station_year_overrides(station, year)
        .get("prefilter_min_amp_dB", {})
        .get(layer)
    )
    if threshold == "adaptive":
        return _adaptive_min_amp_db(normalized)
    return threshold


def _station_amplitude_prefilter(
    raw: pd.DataFrame,
    normalized: pd.DataFrame,
    threshold: float | None,
    masked_bands: list[tuple[float, float]] | None = None,
) -> pd.DataFrame:
    """Drop raw echoes below a resolved amplitude floor and inside masked bands.

    Applied to the raw echoes handed to the scaler for a given layer, before
    any clustering/RFI/Viterbi logic runs. A higher floor near the true
    critical frequency will truncate the trace early (weaker there by
    physics), so tune per station only after checking against a known-good
    trace, not blindly. ``masked_bands`` are persistent fixed-frequency RFI
    stripes (MHz, inclusive) that would otherwise get traced as part of the
    real echo.
    """
    if threshold is None and not masked_bands:
        return raw
    keep = pd.Series(True, index=normalized.index)
    if threshold is not None:
        keep &= normalized["amplitude_db"] > threshold
    for lo, hi in masked_bands or []:
        keep &= ~normalized["frequency_mhz"].between(lo, hi)
    return raw[keep.to_numpy()]


def _station_masked_bands(station: str, year: int) -> list[tuple[float, float]]:
    return _station_year_overrides(station, year).get("masked_freq_bands_mhz", [])


def _month_bounds(event_date: dt.date) -> tuple[dt.datetime, dt.datetime]:
    start = dt.datetime(event_date.year, event_date.month, 1)
    if event_date.month == 12:
        end = dt.datetime(event_date.year + 1, 1, 1)
    else:
        end = dt.datetime(event_date.year, event_date.month + 1, 1)
    return start, end - dt.timedelta(seconds=1)


def _parameter_names(db: DidBaseConnector) -> list[str]:
    """Return all usable PARAMCAT names, or an explicit env override."""
    override = os.environ.get("DIDBASE_PARAMS", "").strip()
    if override:
        return [name.strip() for name in override.split(",") if name.strip()]

    catalog = db.scaled_params_catalog()
    names = catalog["NAME"].dropna().astype(str).str.strip()
    return sorted({name for name in names if 0 < len(name) <= 15})


def _write_monthly_median(
    db: DidBaseConnector,
    station: str,
    event_date: dt.date,
    params: list[str],
    output_path: Path,
    bin_minutes: int = 15,
) -> pd.DataFrame:
    """Fetch a month of DIDBASE scaled values and exclude the event date."""
    start, end = _month_bounds(event_date)
    monthly = db.fetch_scaled_params(
        station=station,
        start=start,
        end=end,
        params=params,
        include_quality=True,
        include_provenance=True,
    )
    if monthly.empty:
        raise RuntimeError(f"No DIDBASE scaled data for {station} in {start:%Y-%m}")

    monthly["datetime"] = pd.to_datetime(monthly["datetime"], errors="coerce")
    monthly = monthly.dropna(subset=["datetime"])
    monthly = monthly[monthly["datetime"].dt.date != event_date].copy()
    if monthly.empty:
        raise RuntimeError(f"No non-eclipse baseline data for {station} in {start:%Y-%m}")

    time_of_day = monthly["datetime"] - monthly["datetime"].dt.normalize()
    bin_size = pd.Timedelta(minutes=bin_minutes)
    monthly["time_of_day"] = (time_of_day / bin_size).round().astype(int) * bin_size

    numeric_params = []
    for param in params:
        if param in monthly.columns and pd.api.types.is_numeric_dtype(monthly[param]):
            numeric_params.append(param)

    if not numeric_params:
        raise RuntimeError(f"No numeric DIDBASE parameters returned for {station}")

    grouped = monthly.groupby("time_of_day", sort=True)
    median = grouped[numeric_params].median().add_suffix("_median")
    count = grouped[numeric_params].count().add_suffix("_n")
    result = pd.concat([median, count], axis=1).reset_index()
    result.insert(0, "station", station)
    result.insert(1, "event_year", event_date.year)
    result["time_of_day"] = result["time_of_day"].dt.total_seconds() / 3600.0
    result.to_csv(output_path, index=False, float_format="%g")
    logger.info(f"Saved baseline median: {output_path} ({len(result)} bins)")
    return result


def _median_path(station: str, year: int) -> Path:
    return DATA_DIR / f"{station}_{year}_median.csv"


def _scaled_path(station: str, year: int) -> Path:
    return DATA_DIR / f"{station}_{year}_scaled.csv"


def _is_nonempty_file(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0


def _cached_raw_files(station: str, event_date: dt.date) -> list[Path]:
    """Return cached DIDBASE raw files, ignoring unrelated cache artifacts."""
    event_dir = RAW_DIR / station / event_date.isoformat()
    if not event_dir.is_dir():
        return []
    return sorted(
        path for path in event_dir.iterdir()
        if path.is_file() and path.suffix.upper() in {".RSF", ".MMM"}
        and path.stat().st_size > 0
    )


def _download_event_raw(
    db: DidBaseConnector,
    station: str,
    event_date: dt.date,
    raw_format: str,
    overwrite: bool,
) -> list[Path]:
    event_dir = RAW_DIR / station / event_date.isoformat()
    event_dir.mkdir(parents=True, exist_ok=True)
    existing = _cached_raw_files(station, event_date)
    if existing and not overwrite:
        logger.info(f"Using cached raw files: {station} {event_date} ({len(existing)})")
        return existing

    start = dt.datetime.combine(event_date, dt.time.min)
    end = start + dt.timedelta(days=1) - dt.timedelta(seconds=1)
    formats = [raw_format.upper()] if raw_format.lower() != "auto" else ["RSF", "MMM"]
    for fmt in formats:
        paths = db.fetch_files(
            station=station,
            start=start,
            end=end,
            fmt=fmt,
            out_dir=event_dir,
        )
        if paths:
            logger.info(f"Downloaded {len(paths)} {fmt} files for {station} {event_date}")
            return paths
        logger.warning(f"No {fmt} files available for {station} {event_date}")
    return []


def _extract_scaled_file(path: Path, station: str) -> dict | None:
    suffix = path.suffix.upper()
    extractor_cls = {".RSF": RsfExtractor, ".MMM": ModMaxExtractor}.get(suffix)
    if extractor_cls is None:
        return None

    extractor = extractor_cls(
        str(path),
        extract_time_from_name=True,
        extract_stn_from_name=True,
    )
    extractor.extract()
    raw = extractor.to_pandas()
    if raw.empty:
        return None
    normalized = normalize_echo_dataframe(raw)

    location = STATIONS[station]
    timestamp = getattr(extractor, "date", None)
    if timestamp is None:
        return None

    year = timestamp.year
    sza = utils._solar_zenith_angle(
        timestamp, location["lat"], location["lon"]
    )
    regime = utils._regime_for_sza(sza)
    params = utils._SZA_PARAM_SETS[regime]
    if _station_year_overrides(station, year).get("e_skip_sza_f_hi_clamp"):
        e_layers = _layer_from_sza_params("E", params, sza=None)
    else:
        e_layers = _layer_from_sza_params("E", params, sza=sza)
    f2_layers = _layer_from_sza_params("F2", params, sza=sza)
    f2_f_lo_override = _station_year_overrides(station, year).get("f2_f_lo_mhz")
    if f2_f_lo_override is not None:
        f2_layers[0]["f_lo"] = f2_f_lo_override
    rfi_filter = (
        VerticalRfiFilter(
            max_echoes_per_bin=RSF_RFI_MAX_ECHOES_PER_BIN,
            freq_bin=RSF_RFI_FREQ_BIN_MHZ,
        )
        if _station_rfi_enabled(station, year)
        else None
    )

    masked_bands = _station_masked_bands(station, year)
    e_amp_threshold = _station_amplitude_threshold(normalized, station, year, "E")
    f2_amp_threshold = _station_amplitude_threshold(normalized, station, year, "F2")

    if _normal_e_layer_enabled(sza):
        e_result, e_details = scale_dataframe(
            _station_amplitude_prefilter(raw, normalized, e_amp_threshold, masked_bands),
            mode="O",
            method="viterbi",
            layers=e_layers,
            scaler_kwargs=_scaler_kwargs_for_layer(
                "E", params, E_SCALER_KWARGS, rfi_filter=rfi_filter
            ),
            return_details=True,
        )
    else:
        e_result = None
        e_details = {"trace_metrics": {}}

    f2_result, f2_details = scale_dataframe(
        _station_amplitude_prefilter(raw, normalized, f2_amp_threshold, masked_bands),
        mode="O",
        method="viterbi",
        layers=f2_layers,
        scaler_kwargs=_scaler_kwargs_for_layer(
            "F2", params, F2_SCALER_KWARGS, rfi_filter=rfi_filter
        ),
        return_details=True,
    )

    row = f2_result.to_series().to_dict()
    if e_result is not None:
        row["foE_mhz"] = e_result.foE_mhz
        row["h_prime_E_km"] = e_result.h_prime_E_km
        # f2_result (an F2-only Viterbi result) always reports flag_E_detected=False
        # since it never scaled E; reflect the actual E scaling outcome here.
        row["flag_E_detected"] = pd.notna(e_result.foE_mhz)
    else:
        row["foE_mhz"] = float("nan")
        row["h_prime_E_km"] = float("nan")
        row["flag_E_detected"] = False
    row.update(
        {
            "time": timestamp,
            "station": station,
            "source_file": path.name,
            "raw_format": suffix.lstrip("."),
            "n_raw_echoes": len(raw),
            "sza_deg": sza,
            "sza_regime": regime,
            "tuning_preset": RSF_TUNING_PRESET,
            "e_layer_enabled": _normal_e_layer_enabled(sza),
            "E_f_hi_effective": e_layers[0]["f_hi"],
            "F2_f_hi_effective": f2_layers[0]["f_hi"],
            "E_prefilter_min_amp_dB": e_amp_threshold,
            "F2_prefilter_min_amp_dB": f2_amp_threshold,
        }
    )
    for details in (e_details, f2_details):
        for key, value in details.get("trace_metrics", {}).items():
            row[key] = value

    max_trace_gap = _station_max_trace_gap(station, year)
    bridge_gaps = _station_bridge_gaps(station, year)
    _trim_trailing_gap_contamination(row, e_details, "E", "foE_mhz", "h_prime_E_km", max_trace_gap, bridge_gaps)
    _trim_trailing_gap_contamination(row, f2_details, "F2", "foF2_mhz", "h_prime_F2_km", max_trace_gap, bridge_gaps)

    _apply_reliability_gate(row, "E", "foE_mhz", "h_prime_E_km", "flag_E_detected")
    _apply_reliability_gate(row, "F2", "foF2_mhz", "h_prime_F2_km", "flag_F2_detected")

    keep_diagnostic = VITERBI_UT_START <= _decimal_ut(pd.Timestamp(timestamp)) < VITERBI_UT_END
    return {
        "row": row,
        "normalized": normalized if keep_diagnostic else None,
        "e_result": e_result,
        "f2_result": f2_result,
        "e_details": e_details if keep_diagnostic else {"trace_metrics": {}},
        "f2_details": f2_details if keep_diagnostic else {"trace_metrics": {}},
    }


def _decimal_ut(timestamp: pd.Timestamp) -> float:
    return timestamp.hour + timestamp.minute / 60.0 + timestamp.second / 3600.0


def _trim_trace_by_gap(
    trace: pd.DataFrame,
    max_gap_mhz: float = RSF_MAX_TRACE_FREQ_GAP_MHZ,
    bridge: bool = False,
) -> pd.DataFrame:
    """Keep the largest frequency-contiguous segment of a trace.

    A fixed-height multipath/artifact echo (or a single isolated stray
    point) can sit apart from the real trace with a frequency gap wider
    than ``max_gap_mhz`` between them, on either side -- e.g. AL945's
    artifact trails the real trace, while BC840 2023-10-14 17:35/17:40 UT
    had a single isolated point *before* a 1.35 MHz gap followed by the
    real 113-point trace. Always keeping "everything before the first gap"
    silently discarded the real trace in that second case, so split on
    every gap and keep the largest segment instead.

    With ``bridge=True``, neighbouring segments that are both at least
    ``RSF_BRIDGE_MIN_SEGMENT_POINTS`` long and whose heights differ by at
    most ``RSF_BRIDGE_MAX_HEIGHT_STEP_KM`` across the gap are merged first,
    so a real trace with an internal dropout (AU930 2023-10-14 15:15 UT)
    isn't reduced to its longer half.
    """
    if trace is None or trace.empty or len(trace) < 2:
        return trace
    trace = trace.sort_values("frequency_mhz").reset_index(drop=True)
    gaps = trace["frequency_mhz"].diff()
    cut_positions = gaps[gaps > max_gap_mhz].index.tolist()
    if not cut_positions:
        return trace
    bounds = [0] + cut_positions + [len(trace)]
    segments = [trace.iloc[bounds[i]:bounds[i + 1]] for i in range(len(bounds) - 1)]
    if bridge:
        merged = [segments[0]]
        for seg in segments[1:]:
            prev = merged[-1]
            step = abs(float(seg["height_km"].iloc[0]) - float(prev["height_km"].iloc[-1]))
            if (
                len(prev) >= RSF_BRIDGE_MIN_SEGMENT_POINTS
                and len(seg) >= RSF_BRIDGE_MIN_SEGMENT_POINTS
                and step <= RSF_BRIDGE_MAX_HEIGHT_STEP_KM
            ):
                merged[-1] = pd.concat([prev, seg])
            else:
                merged.append(seg)
        segments = merged
    return max(segments, key=len)


def _trim_trailing_gap_contamination(
    row: dict, details: dict, layer: str, fo_key: str, hprime_key: str, max_gap_mhz: float,
    bridge: bool = False,
) -> None:
    trace = details.get("traces", {}).get(layer)
    if trace is None or trace.empty:
        return
    trimmed = _trim_trace_by_gap(trace, max_gap_mhz, bridge)
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


def _apply_reliability_gate(
    row: dict, layer: str, fo_key: str, hprime_key: str, flag_key: str
) -> None:
    """Null out a layer's foX/h'X when its trace reliability is too low.

    The Viterbi tracer can lock onto a handful of noise points and still
    report a foX value; ``{layer}_trace_reliability_score`` is a much better
    quality signal than "a trace was found at all", so gate on it instead of
    trusting every detection.
    """
    score = row.get(f"{layer}_trace_reliability_score")
    npts = row.get(f"{layer}_trace_n_points")
    insufficient_points = npts is None or pd.isna(npts) or npts < RSF_MIN_TRACE_POINTS
    if score is None or pd.isna(score) or score < RSF_MIN_TRACE_RELIABILITY or insufficient_points:
        row[fo_key] = float("nan")
        row[hprime_key] = float("nan")
        row[flag_key] = False


def _plot_viterbi_diagnostic(item: dict, station: str, output_path: Path) -> None:
    """Plot raw O-mode, cleaned E/F2 echoes, and separate Viterbi traces."""
    raw = item["normalized"]
    raw_o = raw[raw["mode"] == "O"].copy() if "mode" in raw.columns else raw
    fig, ax = plt.subplots(figsize=(7.2, 5.2), constrained_layout=True)

    if not raw_o.empty:
        plot_raw = raw_o[
            raw_o["height_km"].between(70, 600)
            & raw_o["frequency_mhz"].between(1, 14)
        ]
        pivot = plot_raw.pivot_table(
            index="height_km", columns="frequency_mhz",
            values="amplitude_db", aggfunc="max",
        ).sort_index().sort_index(axis=1)
        if not pivot.empty:
            mesh = ax.pcolormesh(
                pivot.columns.values, pivot.index.values, pivot.values,
                cmap="viridis", shading="auto", rasterized=True,
            )
            fig.colorbar(mesh, ax=ax, pad=0.01, fraction=0.045, label="Amplitude [dB]")

    plot_year = pd.Timestamp(item["row"]["time"]).year
    max_trace_gap = _station_max_trace_gap(station, plot_year)
    bridge_gaps = _station_bridge_gaps(station, plot_year)
    for details, color, label in [
        (item["e_details"], "white", "E Viterbi"),
        (item["f2_details"], "black", "F2 Viterbi"),
    ]:
        trace = details.get("traces", {}).get("E" if label.startswith("E") else "F2", pd.DataFrame())
        trace = _trim_trace_by_gap(trace, max_trace_gap, bridge_gaps)
        if not trace.empty:
            trace = trace.sort_values("frequency_mhz")
            ax.plot(
                trace["frequency_mhz"], trace["height_km"],
                color=color, lw=2.2, zorder=5,
                path_effects=[pe.withStroke(linewidth=3.6, foreground="black" if color == "white" else "white")],
                label=f"{label} ({len(trace)} pts)",
            )

    timestamp = pd.Timestamp(item["row"]["time"])
    foE = item["row"].get("foE_mhz", float("nan"))
    foF2 = item["row"].get("foF2_mhz", float("nan"))
    annotation = [
        f"SZA={item['row'].get('sza_deg', float('nan')):.1f} deg, "
        f"{item['row'].get('sza_regime', '')}",
        f"E scaling={'on' if item['row'].get('e_layer_enabled', False) else 'off'}",
        f"foE={foE:.2f} MHz" if pd.notna(foE) else "foE=NaN",
        f"foF2={foF2:.2f} MHz" if pd.notna(foF2) else "foF2=NaN",
    ]
    ax.text(
        0.02, 0.98, "\n".join(annotation), transform=ax.transAxes,
        va="top", ha="left", fontsize=9,
        bbox=dict(facecolor="white", edgecolor="0.8", alpha=0.85),
    )
    ax.set_title(f"{station} {timestamp:%Y-%m-%d %H:%M:%S} UT", loc="left", fontsize=11)
    ax.set_xlabel("Frequency [MHz]")
    ax.set_ylabel("Virtual height [km]")
    ax.set_xlim(1, 14)
    ax.set_ylim(70, 600)
    ax.grid(True, ls=":", lw=0.45, alpha=0.55)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.75, edgecolor="none")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=170, bbox_inches="tight")
    plt.close(fig)
    logger.info(f"Saved Viterbi diagnostic: {output_path}")


def _write_event_scaled(paths: list[Path], station: str, event_date: dt.date, output_path: Path) -> pd.DataFrame:
    results = Parallel(n_jobs=RSF_SCALE_JOBS, backend="loky")(
        delayed(_scale_one_file_safe)(path, station) for path in paths
    )
    rows = [item["row"] for item in results if item is not None]

    scaled = pd.DataFrame(rows)
    if not scaled.empty and "time" in scaled.columns:
        scaled["time"] = pd.to_datetime(scaled["time"], errors="coerce")
        scaled = scaled.sort_values("time").reset_index(drop=True)
    scaled.to_csv(output_path, index=False, float_format="%g")
    logger.info(f"Saved eclipse-day scaled data: {output_path} ({len(scaled)} rows)")

    diagnostic_items = [
        item for item in results
        if item is not None
        and VITERBI_UT_START <= _decimal_ut(pd.Timestamp(item["row"]["time"])) < VITERBI_UT_END
    ]
    diagnostic_items = sorted(diagnostic_items, key=lambda item: item["row"]["time"])
    if VITERBI_KEEP_EVERY > 1:
        diagnostic_items = diagnostic_items[::VITERBI_KEEP_EVERY]
    for item in diagnostic_items:
        timestamp = pd.Timestamp(item["row"]["time"])
        diagnostic_path = (
            VITERBI_DIR / f"{station}_{timestamp.year}"
            / f"viterbi_{station}_{timestamp:%Y%m%d_%H%M%S}.png"
        )
        _plot_viterbi_diagnostic(item, station, diagnostic_path)
    return scaled


def _scale_one_file_safe(path: Path, station: str) -> dict | None:
    try:
        return _extract_scaled_file(path, station)
    except Exception as exc:
        logger.warning(f"Failed to scale {path.name}: {exc}")
        return None


def _run_station_event(
    station: str,
    year: int,
    raw_format: str,
    overwrite: bool,
) -> None:
    event_date = EVENTS[year]
    scaled_path = _scaled_path(station, year)

    if not overwrite and _is_nonempty_file(scaled_path):
        logger.info(f"Using cached scaled output: {scaled_path}")
        return

    cached_raw = _cached_raw_files(station, event_date)
    if not cached_raw:
        raise RuntimeError(
            f"No cached raw files available for {station} {event_date}; "
            "the DIDBASE preflight did not produce them."
        )
    raw_paths = cached_raw
    logger.info(f"Using cached raw files: {station} {event_date} ({len(raw_paths)})")

    _write_event_scaled(raw_paths, station, event_date, scaled_path)


def run(stations: list[str], years: list[int], raw_format: str, overwrite: bool) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    for station in stations:
        if station not in STATIONS:
            raise ValueError(f"Unsupported station code: {station}")

    jobs = [(station, year) for year in years for station in stations]
    missing_medians = [
        (station, year) for station, year in jobs
        if overwrite or not _is_nonempty_file(_median_path(station, year))
    ]
    missing_raw = [
        (station, year) for station, year in jobs
        if overwrite or not _cached_raw_files(station, EVENTS[year])
    ]
    needs_db = bool(missing_medians or missing_raw)

    params: list[str] | None = None
    if needs_db:
        if not DIDBASE_USER or not DIDBASE_PASSWORD:
            raise RuntimeError(
                "DIDBASE credentials are required because cached medians/raw files "
                "are missing. Set DIDBASE_USER and DIDBASE_PASSWORD."
            )
        db_context = DidBaseConnector(
            DIDBASE_URL, DIDBASE_USER, DIDBASE_PASSWORD, DIDBASE_ROLE
        )
        with db_context as db:
            if missing_medians:
                params = _parameter_names(db)
                logger.info(
                    f"Using {len(params)} DIDBASE scaled parameters for "
                    f"{len(missing_medians)} missing median file(s)"
                )
            for station, year in jobs:
                event_date = EVENTS[year]
                median_path = _median_path(station, year)
                if (station, year) in missing_medians:
                    _write_monthly_median(db, station, event_date, params or [], median_path)
                if (station, year) in missing_raw:
                    _download_event_raw(db, station, event_date, raw_format, overwrite)
    else:
        logger.info("All median and event-day raw files are cached; skipping DIDBASE connection.")

    logger.info(
        f"DIDBASE_WORKERS={DIDBASE_WORKERS}, RSF_SCALE_JOBS={RSF_SCALE_JOBS}, "
        f"overwrite={overwrite}"
    )

    with ThreadPoolExecutor(max_workers=max(1, DIDBASE_WORKERS)) as executor:
        futures = [
            executor.submit(
                _run_station_event,
                station,
                year,
                raw_format,
                overwrite,
            )
            for station, year in jobs
        ]
        for (station, year), future in zip(jobs, futures):
            try:
                future.result()
            except Exception:
                logger.exception(f"Failed pipeline for {station} {year}")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stations", nargs="+", default=list(STATIONS))
    parser.add_argument("--years", nargs="+", type=int, default=sorted(EVENTS))
    parser.add_argument("--raw-format", choices=["auto", "RSF", "MMM"], default="auto")
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args.stations, args.years, args.raw_format, args.overwrite)
