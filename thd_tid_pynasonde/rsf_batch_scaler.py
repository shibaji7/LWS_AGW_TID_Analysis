"""Batch-scale RSF ionograms with the generic Digisonde scaler.

DataFrame-native scaler workflow, driven by SZA regime + per-station/year
overrides (``config.toml``):

1. Read one or more RSF files with ``RsfExtractor``.
2. Convert each parsed file to a pandas DataFrame.
3. Run the generic Viterbi scaler on O-mode echoes for the E and F2 layers,
   with amplitude/band prefiltering, trailing-gap trimming, and a trace
   reliability gate.
4. Save a scaled-parameter CSV.

By default this script loops over daily raw RSF folders downloaded by
``didbase_pipeline.py --raw`` under
``Data/didbase_giro_scaled/<year>/<station>/raw/<date>/``. To run on a custom
directory, set ``RSF_SCALER_GLOB``, for example:

    RSF_SCALER_GLOB="/path/to/day/*.RSF" python rsf_batch_scaler.py

Outputs:
    Data/processing/scaled/rsf_scaled_YYYYMMDD_STATION.csv
"""

from __future__ import annotations

import glob
import math
import os
from datetime import timedelta
from pathlib import Path

import pandas as pd
from joblib import Parallel, delayed
from loguru import logger

from pynasonde.digisonde.parsers.rsf import RsfExtractor
from pynasonde.digisonde.scaler import normalize_echo_dataframe, scale_dataframe
from pynasonde.digisonde.scaler.rfi_filter import VerticalRfiFilter

import utils

ROOT = Path(__file__).resolve().parent
STATION = os.environ.get("RSF_SCALER_STATION", "KAFB")
DATA_ROOT = Path(os.environ.get("RSF_SCALER_DATA_ROOT", str(ROOT / "Data" / "didbase_giro_scaled")))
DATA_GLOB = os.environ.get("RSF_SCALER_GLOB")
START_DATE = pd.Timestamp(os.environ.get("RSF_SCALER_START_DATE", "2024-04-01"))
N_DAYS = int(os.environ.get("RSF_SCALER_N_DAYS", "60"))
N_JOBS = int(os.environ.get("RSF_SCALER_N_JOBS", "20"))
BATCH_SIZE = int(os.environ.get("RSF_SCALER_BATCH_SIZE", "32"))
MAX_FILES = int(os.environ.get("RSF_SCALER_MAX_FILES", "0"))
_RSF_CFG = utils.rsf_scaler_config()
TUNING_PRESET = os.environ.get("RSF_SCALER_TUNING_PRESET", _RSF_CFG["tuning_preset"]).lower()
E_LAYER_MAX_SZA = float(os.environ.get("RSF_SCALER_E_LAYER_MAX_SZA", _RSF_CFG["e_layer_max_sza"]))
RSF_SCALER_OVERWRITE = os.environ.get("RSF_SCALER_OVERWRITE", "0").strip().lower() in {
    "1", "true", "yes", "y", "on"
}
REQUIRED_OUTPUT_COLUMNS = (
    "time",
    "foE_mhz",
    "foF2_mhz",
    "source_file",
    "F2_trace_reliability_score",
)
OUT_DIR_ENV = os.environ.get("RSF_SCALER_OUT_DIR")
if OUT_DIR_ENV:
    OUT_DIR = Path(OUT_DIR_ENV)
elif TUNING_PRESET == "balanced":
    OUT_DIR = ROOT / "Data" / "processing" / "scaled"
else:
    OUT_DIR = ROOT / "Data" / "processing" / "scaled_tuning" / TUNING_PRESET


STATION_LOCATIONS = {
    # Kirtland AFB / KR835
    "KAFB": (35.04, -106.61),
    "KR835": (35.04, -106.61),
    # White Sands Missile Range / WS833
    "WSMR": (32.40, -106.40),
    "WS833": (32.40, -106.40),
}
_LAT, _LON = STATION_LOCATIONS.get(STATION, STATION_LOCATIONS["KAFB"])

E_SCALER_KWARGS = {
    "dbscan_eps": 0.035,
    "dbscan_min_samples": 10,
    "min_echoes": 3,
    "min_amp_dB": 12.0,
    "sigma_trans_km": 10.0,
    "consecutive_gap": 8,
    "max_jump_km": 25.0,
}

F2_SCALER_KWARGS = {
    "dbscan_eps": 0.02,
    "dbscan_min_samples": 5,
    "min_echoes": 5,
    "min_amp_dB": 24.0,
    "sigma_trans_km": 18.0,
    "max_jump_km": 22.0,
}

TUNING_PRESETS = {
    "balanced": {
        "description": "Current production behavior with SZA-driven parameters.",
        "E": {},
        "F2": {},
    },
    "e_sensitive": {
        "description": "More permissive E-region detection for weak/lower echoes.",
        "E": {
            "scaler_kwargs": {
                "dbscan_eps": 0.020,
                "dbscan_min_samples": 3,
                "min_amp_dB": 18.0,
                "min_echoes": 2,
                "consecutive_gap": 10,
                "max_jump_km": 28.0,
                "sigma_trans_km": 12.0,
            },
            "layer": {
                "h_lo": 90.0,
                "h_hi": 170.0,
                "f_lo": 1.5,
                "f_hi": 5.5,
                "sigma_trans_km": 12.0,
                "max_jump_km": 28.0,
            },
        },
        "F2": {},
    },
    "f2_smooth": {
        "description": "More continuous F2 tracing with slightly stronger smoothing.",
        "E": {},
        "F2": {
            "scaler_kwargs": {
                "dbscan_eps": 0.015,
                "dbscan_min_samples": 5,
                "min_amp_dB": 25.0,
                "min_echoes": 5,
                "sigma_trans_km": 28.0,
                "max_jump_km": 45.0,
            },
            "layer": {
                "sigma_trans_km": 32.0,
                "max_jump_km": 50.0,
            },
        },
    },
    "strict": {
        "description": "Higher-confidence traces; useful for rejecting noisy false locks.",
        "E": {
            "scaler_kwargs": {
                "dbscan_eps": 0.010,
                "dbscan_min_samples": 5,
                "min_amp_dB": 27.0,
                "min_echoes": 4,
                "max_jump_km": 18.0,
            },
            "layer": {
                "h_lo": 100.0,
                "h_hi": 155.0,
                "f_lo": 2.0,
                "f_hi": 5.0,
                "max_jump_km": 18.0,
            },
        },
        "F2": {
            "scaler_kwargs": {
                "dbscan_eps": 0.010,
                "dbscan_min_samples": 6,
                "min_amp_dB": 28.0,
                "min_echoes": 6,
                "sigma_trans_km": 20.0,
                "max_jump_km": 30.0,
            },
            "layer": {
                "sigma_trans_km": 22.0,
                "max_jump_km": 35.0,
            },
        },
    },
    "permissive": {
        "description": "Broad capture pass for diagnosing missed E/F2 traces.",
        "E": {
            "scaler_kwargs": {
                "dbscan_eps": 0.030,
                "dbscan_min_samples": 3,
                "min_amp_dB": 16.0,
                "min_echoes": 2,
                "consecutive_gap": 12,
                "sigma_trans_km": 14.0,
                "max_jump_km": 32.0,
            },
            "layer": {
                "h_lo": 90.0,
                "h_hi": 180.0,
                "f_lo": 1.5,
                "f_hi": 6.0,
                "sigma_trans_km": 14.0,
                "max_jump_km": 32.0,
            },
        },
        "F2": {
            "scaler_kwargs": {
                "dbscan_eps": 0.025,
                "dbscan_min_samples": 4,
                "min_amp_dB": 22.0,
                "min_echoes": 4,
                "sigma_trans_km": 35.0,
                "max_jump_km": 60.0,
            },
            "layer": {
                "h_lo": 170.0,
                "h_hi": 560.0,
                "f_lo": 1.8,
                "f_hi": 13.0,
                "sigma_trans_km": 40.0,
                "max_jump_km": 65.0,
            },
        },
    },
}


def _preset_for_layer(name: str) -> dict:
    if TUNING_PRESET not in TUNING_PRESETS:
        raise ValueError(
            f"Unknown RSF_SCALER_TUNING_PRESET={TUNING_PRESET!r}. "
            f"Available presets: {', '.join(sorted(TUNING_PRESETS))}"
        )
    return TUNING_PRESETS[TUNING_PRESET].get(name, {})


def _effective_layer_f_hi(layer: dict, sza: float | None, skip_sza_clamp: bool = False) -> float:
    """Return layer f_hi after optional daylight SZA frequency limiting."""
    f_lo = float(layer["f_lo"])
    f_hi = float(layer["f_hi"])

    name = layer["name"]
    if (
        name == "E"
        and not skip_sza_clamp
        and sza is not None
        and pd.notna(sza)
        and 0.0 < float(sza) < 90.0
    ):
        f_hi = min(f_hi, 6.0 * math.cos(math.radians(float(sza))))

    if f_hi < f_lo:
        f_hi = f_lo + 2.0

    return f_hi


def _layer_from_sza_params(name: str, params: dict, sza: float | None = None, skip_sza_clamp: bool = False) -> list[dict]:
    layer = dict(params[name])
    layer.update(_preset_for_layer(name).get("layer", {}))
    layer["name"] = name
    layer["f_hi"] = _effective_layer_f_hi(layer, sza, skip_sza_clamp=skip_sza_clamp)
    return [layer]


def _scaler_kwargs_from_sza_params(params: dict) -> dict:
    return {
        "dbscan_eps": params["dbscan_eps"],
        "dbscan_min_samples": params["dbscan_min_samples"],
        "min_amp_dB": params["min_amp_dB"],
    }


def _scaler_kwargs_for_layer(
    name: str,
    params: dict,
    base_kwargs: dict,
    rfi_filter: VerticalRfiFilter | None = None,
) -> dict:
    kwargs = {
        **base_kwargs,
        **_scaler_kwargs_from_sza_params(params),
    }
    kwargs.update(_preset_for_layer(name).get("scaler_kwargs", {}))
    if rfi_filter is not None:
        kwargs["rfi_filter"] = rfi_filter
    return kwargs


def _normal_e_layer_enabled(sza: float) -> bool:
    """Return True when normal E-layer scaling is allowed by solar illumination."""
    return pd.notna(sza) and sza < E_LAYER_MAX_SZA


def _add_trace_metrics(row: pd.Series, details: dict) -> None:
    """Copy scaler trace reliability metrics into one output row."""
    for key, value in details.get("trace_metrics", {}).items():
        row[key] = value


def _station_from_path(fpath: str) -> str:
    return Path(fpath).name.split("_")[0]


def _process_one_file(fpath: str) -> dict | None:
    """Parse and scale one RSF file, applying station/year overrides."""
    try:
        ext = RsfExtractor(fpath, extract_time_from_name=True, extract_stn_from_name=True)
        ext.extract()
        raw_df = ext.to_pandas()
        if raw_df.empty:
            logger.warning(f"Skipping empty RSF file: {fpath}")
            return None

        normalized_all = normalize_echo_dataframe(raw_df)
        station = getattr(ext, "stn_code", None) or _station_from_path(fpath)
        year = ext.date.year
        overrides = utils.station_year_overrides(station, year)

        sza = utils._solar_zenith_angle(ext.date, _LAT, _LON)
        regime = utils._regime_for_sza(sza)
        p = utils._SZA_PARAM_SETS[regime]
        logger.info(
            f"{station} {ext.date}: SZA = {sza:.1f} deg -> regime = '{regime}'  "
            f"(eps={p['dbscan_eps']}, min_samples={p['dbscan_min_samples']}, "
            f"min_amp={p['min_amp_dB']} dB, tuning={TUNING_PRESET})"
        )

        skip_e_clamp = bool(overrides.get("e_skip_sza_f_hi_clamp"))
        e_layer_enabled = _normal_e_layer_enabled(sza)
        e_layers = _layer_from_sza_params("E", p, sza=sza, skip_sza_clamp=skip_e_clamp)
        f2_layers = _layer_from_sza_params("F2", p, sza=sza)
        f2_f_lo_override = overrides.get("f2_f_lo_mhz")
        if f2_f_lo_override is not None:
            f2_layers[0]["f_lo"] = f2_f_lo_override

        rfi_enabled = overrides.get("rfi_enabled", _RSF_CFG["rfi_enabled"])
        rfi_filter = (
            VerticalRfiFilter(
                max_echoes_per_bin=_RSF_CFG["rfi_max_echoes_per_bin"],
                freq_bin=_RSF_CFG["rfi_freq_bin_mhz"],
            )
            if rfi_enabled
            else None
        )

        masked_bands = [tuple(band) for band in overrides.get("masked_freq_bands_mhz", [])]
        e_amp_threshold = utils.station_amplitude_threshold(normalized_all, station, year, "E")
        f2_amp_threshold = utils.station_amplitude_threshold(normalized_all, station, year, "F2")
        e_raw = utils.station_amplitude_prefilter(raw_df, normalized_all, e_amp_threshold, masked_bands)
        f2_raw = utils.station_amplitude_prefilter(raw_df, normalized_all, f2_amp_threshold, masked_bands)

        if e_layer_enabled:
            e_result, e_details = scale_dataframe(
                e_raw,
                mode="O",
                method="viterbi",
                layers=e_layers,
                scaler_kwargs=_scaler_kwargs_for_layer(
                    "E", p, E_SCALER_KWARGS, rfi_filter=rfi_filter
                ),
                return_details=True,
            )
        else:
            logger.info(
                f"Skipping normal E-layer scaling because SZA={sza:.1f} deg "
                f">= RSF_SCALER_E_LAYER_MAX_SZA={E_LAYER_MAX_SZA:.1f} deg"
            )
            e_result = None
            e_details = {"cleaned": pd.DataFrame(), "traces": {}}

        f2_result, f2_details = scale_dataframe(
            f2_raw,
            mode="O",
            method="viterbi",
            layers=f2_layers,
            scaler_kwargs=_scaler_kwargs_for_layer(
                "F2", p, F2_SCALER_KWARGS, rfi_filter=rfi_filter
            ),
            return_details=True,
        )

        row = f2_result.to_series()
        _add_trace_metrics(row, e_details)
        _add_trace_metrics(row, f2_details)
        row["foE_mhz"] = e_result.foE_mhz if e_result is not None else float("nan")
        row["h_prime_E_km"] = (
            e_result.h_prime_E_km if e_result is not None else float("nan")
        )
        row["time"] = getattr(ext, "date", pd.NaT)
        row["source_file"] = Path(fpath).name
        row["n_raw_echoes"] = len(raw_df)
        row["n_cleaned_echoes_E"] = len(e_details["cleaned"])
        row["n_cleaned_echoes_F2"] = len(f2_details["cleaned"])
        row["station"] = station
        row["sza_deg"] = sza
        row["sza_regime"] = regime
        row["tuning_preset"] = TUNING_PRESET
        row["e_layer_enabled"] = e_layer_enabled
        row["e_layer_max_sza"] = E_LAYER_MAX_SZA
        row["E_f_hi_effective"] = e_layers[0]["f_hi"]
        row["F2_f_hi_effective"] = f2_layers[0]["f_hi"]
        row["E_prefilter_min_amp_dB"] = e_amp_threshold
        row["F2_prefilter_min_amp_dB"] = f2_amp_threshold

        max_trace_gap = overrides.get("max_trace_freq_gap_mhz", _RSF_CFG["max_trace_freq_gap_mhz"])
        bridge_gaps = bool(overrides.get("bridge_trace_gaps", False))
        row = row.to_dict()
        utils.trim_trailing_gap_contamination(row, e_details, "E", "foE_mhz", "h_prime_E_km", max_trace_gap, bridge_gaps)
        utils.trim_trailing_gap_contamination(row, f2_details, "F2", "foF2_mhz", "h_prime_F2_km", max_trace_gap, bridge_gaps)
        utils.apply_reliability_gate(row, "E", "foE_mhz", "h_prime_E_km", "flag_E_detected")
        utils.apply_reliability_gate(row, "F2", "foF2_mhz", "h_prime_F2_km", "flag_F2_detected")

        if "height_km" in e_details["cleaned"].columns:
            e_plot = e_details["cleaned"][
                e_details["cleaned"]["height_km"].between(90.0, 190.0)
            ]
        else:
            e_plot = pd.DataFrame()
        f2_plot = f2_details["cleaned"][
            f2_details["cleaned"]["height_km"].between(170.0, 560.0)
        ]

        return {
            "row": row,
            "normalized": normalized_all,
            "cleaned": pd.concat([e_plot, f2_plot], ignore_index=True),
            "traces": {
                "E": e_details["traces"].get("E", pd.DataFrame()),
                "F2": f2_details["traces"].get("F2", pd.DataFrame()),
            },
        }
    except Exception as exc:
        logger.warning(f"Failed to scale {fpath}: {exc}")
        return None


def _date_range(start: pd.Timestamp, n_days: int) -> list[pd.Timestamp]:
    """Return midnight-normalized daily timestamps."""
    start = pd.Timestamp(start).normalize()
    return [start + timedelta(days=i) for i in range(n_days)]


def _files_for_date(day: pd.Timestamp, station: str = STATION) -> list[str]:
    """Select RSF files from one raw-download date directory for one station."""
    date_dir = DATA_ROOT / str(day.year) / station / "raw" / day.strftime("%Y-%m-%d")
    files = sorted(str(path) for path in date_dir.glob("*.RSF"))
    if MAX_FILES > 0:
        files = files[:MAX_FILES]
    return files


def _station_from_scaled(scaled: pd.DataFrame) -> str:
    """Return parsed URSI station code when available."""
    station = STATION
    if not scaled.empty and "station" in scaled.columns:
        parsed = scaled["station"].dropna()
        if not parsed.empty:
            station = str(parsed.mode().iloc[0])
    return station


def _station_from_files(files: list[str]) -> str:
    """Return parsed URSI station code from an RSF file name when available."""
    if not files:
        return STATION
    return _station_from_path(files[0])


def _existing_output_is_valid(out_csv: Path, expected_rows: int) -> bool:
    """Return whether an existing daily output is complete enough to reuse."""
    if not out_csv.is_file() or out_csv.stat().st_size == 0:
        return False

    try:
        header = pd.read_csv(out_csv, nrows=0)
        if not set(REQUIRED_OUTPUT_COLUMNS).issubset(header.columns):
            return False
        n_rows = sum(
            len(chunk)
            for chunk in pd.read_csv(out_csv, usecols=["time"], chunksize=4096)
        )
    except (OSError, ValueError, pd.errors.ParserError):
        return False

    return n_rows == expected_rows


def _output_path_for_day(day: pd.Timestamp, scaled: pd.DataFrame | None = None, station: str | None = None) -> Path:
    """Build daily output path using the parsed URSI station code when present."""
    station = station or (_station_from_scaled(scaled) if scaled is not None else STATION)
    if MAX_FILES > 0:
        return (
            OUT_DIR / "test" / f"rsf_scaled_{day.strftime('%Y%m%d')}_{station}_max{MAX_FILES}.csv"
        )
    return OUT_DIR / f"rsf_scaled_{day.strftime('%Y%m%d')}_{station}.csv"


def _scale_file_batch(files: list[str]) -> pd.DataFrame:
    """Scale one file batch and return sorted rows."""
    results = Parallel(n_jobs=N_JOBS, backend="loky", verbose=5)(
        delayed(_process_one_file)(fpath) for fpath in files
    )
    scaled = pd.DataFrame(
        [result["row"] for result in results if result is not None]
    )
    if not scaled.empty and "time" in scaled.columns:
        scaled = scaled.sort_values("time").reset_index(drop=True)
    return scaled


def _scale_files(files: list[str], out_csv: Path | None = None) -> pd.DataFrame:
    """Scale a file list, optionally writing incremental batch results."""
    logger.info(
        f"Processing {len(files)} RSF files with N_JOBS={N_JOBS}, "
        f"BATCH_SIZE={BATCH_SIZE}"
    )
    if BATCH_SIZE <= 0:
        batches = [files]
    else:
        batches = [
            files[start : start + BATCH_SIZE]
            for start in range(0, len(files), BATCH_SIZE)
        ]

    all_scaled: list[pd.DataFrame] = []
    wrote_header = False
    for batch_no, batch in enumerate(batches, start=1):
        logger.info(f"Scaling batch {batch_no}/{len(batches)} ({len(batch)} RSF files)")
        batch_scaled = _scale_file_batch(batch)
        if batch_scaled.empty:
            logger.warning(f"Batch {batch_no}/{len(batches)} produced no scaled rows")
            continue

        all_scaled.append(batch_scaled)
        if out_csv is not None:
            out_csv.parent.mkdir(parents=True, exist_ok=True)
            batch_scaled.to_csv(
                out_csv,
                mode="w" if not wrote_header else "a",
                header=not wrote_header,
                index=False,
                float_format="%g",
            )
            wrote_header = True
            logger.info(
                f"Saved incremental rows through batch {batch_no}/{len(batches)}: "
                f"{out_csv}"
            )

    if not all_scaled:
        return pd.DataFrame()

    scaled = pd.concat(all_scaled, ignore_index=True)
    if "time" in scaled.columns:
        scaled = scaled.sort_values("time").reset_index(drop=True)
    if out_csv is not None:
        scaled.to_csv(out_csv, index=False, float_format="%g")
        logger.info(f"Saved final sorted {len(scaled)} scaled rows: {out_csv}")
    return scaled


def _log_scaled_summary(label: str, scaled: pd.DataFrame) -> None:
    """Log compact E/F2 yield and range diagnostics for tuning comparisons."""
    if scaled.empty:
        logger.warning(f"{label}: no scaled rows")
        return

    pieces = [f"{label}: rows={len(scaled)}"]
    for col, name in [("foE_mhz", "foE"), ("foF2_mhz", "foF2")]:
        if col not in scaled.columns:
            continue
        values = pd.to_numeric(scaled[col], errors="coerce").dropna()
        if values.empty:
            pieces.append(f"{name}=0 valid")
            continue
        pieces.append(
            f"{name}={len(values)} valid "
            f"({values.quantile(0.05):.2f}-{values.quantile(0.95):.2f} MHz p05-p95)"
        )
    for col, name in [("h_prime_E_km", "h'E"), ("h_prime_F2_km", "h'F2")]:
        if col not in scaled.columns:
            continue
        values = pd.to_numeric(scaled[col], errors="coerce").dropna()
        if values.empty:
            continue
        pieces.append(
            f"{name}={values.quantile(0.05):.0f}-{values.quantile(0.95):.0f} km p05-p95"
        )
    logger.info(" | ".join(pieces))


if __name__ == "__main__":
    logger.info(
        f"RSF scaler station={STATION} lat={_LAT:.2f}, lon={_LON:.2f}, "
        f"N_JOBS={N_JOBS}, data_root={DATA_ROOT}, "
        f"overwrite={int(RSF_SCALER_OVERWRITE)}"
    )

    OUT_DIR.mkdir(parents=True, exist_ok=True)

    if DATA_GLOB:
        files = sorted(glob.glob(DATA_GLOB, recursive=True))
        if MAX_FILES > 0:
            files = files[:MAX_FILES]
        if not files:
            raise FileNotFoundError(
                f"No RSF files matched {DATA_GLOB!r}. Set RSF_SCALER_GLOB to a "
                "directory containing Digisonde .RSF files."
            )
        if MAX_FILES > 0:
            out_csv = OUT_DIR / "test" / f"rsf_scaled_{STATION}_max{MAX_FILES}.csv"
        else:
            out_csv = OUT_DIR / f"rsf_scaled_{STATION}.csv"
        out_csv.parent.mkdir(parents=True, exist_ok=True)
        if not RSF_SCALER_OVERWRITE and _existing_output_is_valid(out_csv, len(files)):
            logger.info(
                f"Existing valid output found ({len(files)} rows); "
                f"skipping scaling: {out_csv}"
            )
            raise SystemExit(0)
        scaled = _scale_files(files, out_csv=out_csv)
        _log_scaled_summary(f"{STATION} {TUNING_PRESET}", scaled)
        logger.info(f"Completed {len(scaled)} scaled rows: {out_csv}")
    else:
        saved = 0
        skipped = 0
        reused = 0
        for day in _date_range(START_DATE, N_DAYS):
            files = _files_for_date(day)
            if not files:
                skipped += 1
                logger.warning(
                    f"No RSF files for {day.strftime('%Y-%m-%d')} "
                    f"under {DATA_ROOT / str(day.year) / STATION / 'raw'}; skipping."
                )
                continue

            logger.info(
                f"Scaling {day.strftime('%Y-%m-%d')} with {len(files)} RSF files"
            )
            out_csv = _output_path_for_day(day, station=_station_from_files(files))
            out_csv.parent.mkdir(parents=True, exist_ok=True)
            if not RSF_SCALER_OVERWRITE and _existing_output_is_valid(out_csv, len(files)):
                reused += 1
                logger.info(
                    f"Existing valid output found ({len(files)} rows); "
                    f"skipping {day.strftime('%Y-%m-%d')}: {out_csv}"
                )
                continue
            scaled = _scale_files(files, out_csv=out_csv)
            _log_scaled_summary(
                f"{day.strftime('%Y-%m-%d')} {_station_from_files(files)} {TUNING_PRESET}",
                scaled,
            )
            logger.info(f"Completed {len(scaled)} scaled rows: {out_csv}")
            saved += 1

        logger.info(
            f"Daily scaling complete: saved {saved} day files, "
            f"reused {reused} valid files, skipped {skipped} days"
        )
