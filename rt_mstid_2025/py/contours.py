from pathlib import Path
import sys
sys.path.extend([
    str(Path(__file__).resolve().parents[1]),
    str(Path(__file__).resolve().parents[2]),
])

import datetime as dt
import math

from typing import Dict, Iterable, List, Tuple

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.dates import DateFormatter
from matplotlib.lines import Line2D
from matplotlib import cm

import iricore
import utils
import pandas as pd
from pynasonde.digisonde.digi_utils import get_digisonde_info
from pynasonde.digisonde.parsers.sao import SaoExtractor

size = 15
import matplotlib as mpl
import scienceplots
plt.style.use(["science", "ieee"])
plt.rcParams["font.family"] = "sans-serif"
plt.rcParams["font.sans-serif"] = [
    "Tahoma",
    "DejaVu Sans",
    "Lucida Grande",
    "Verdana",
]
mpl.rcParams.update(
    {"xtick.labelsize": size, "ytick.labelsize": size, "font.size": size}
)

# Colours chosen to keep frequencies and heights visually distinct.
COLORS = {
    "foF2": "#1f77b4",  # blue
    "foE": "#d62728",  # red
    "hmF2": "#2ca02c",  # green
    "hEs": "#ff7f0e",  # orange
}

FREQ_BIN_EDGES = np.arange(0.1, 10.1, 1.0)  # 0-1, 1-2, ... 7-8 MHz
FREQ_BIN_CENTERS = FREQ_BIN_EDGES[:-1] + 0.5


def _to_utc_datetime(value) -> dt.datetime:
    """Convert StartTimeUTC strings of the form 'YYYY-MM-DD offset HH:MM:SS.sss' to naive UTC datetimes."""
    if isinstance(value, dt.datetime):
        return value
    try:
        date_part, offset_minutes, time_part = value.split()
    except ValueError:
        raise ValueError(f"Unexpected StartTimeUTC format: {value!r}")

    base_time = dt.datetime.strptime(f"{date_part} {time_part}", "%Y-%m-%d %H:%M:%S.%f")
    tzinfo = dt.timezone(dt.timedelta(minutes=int(offset_minutes)))
    local = base_time.replace(tzinfo=tzinfo).astimezone(dt.timezone.utc).replace(tzinfo=None)
    return base_time


def _extract_peak(freqs: np.ndarray, heights: np.ndarray, mask: np.ndarray) -> Tuple[float, float]:
    """Return peak frequency and associated height within the masked region."""
    if freqs.size == 0 or not np.any(mask):
        return np.nan, np.nan

    idx = np.where(mask)[0]
    sub_freqs = freqs[idx]
    sub_heights = heights[idx]
    valid = ~np.isnan(sub_freqs)
    if not np.any(valid):
        return np.nan, np.nan
    valid_idx = idx[valid]
    best_local_idx = np.argmax(sub_freqs[valid])
    peak_idx = valid_idx[best_local_idx]
    return freqs[peak_idx], heights[peak_idx]


def compute_iri_series(
    times: Iterable[dt.datetime],
    lat: float,
    lon: float,
    alt_profile: Tuple[int, int, int] = (90, 450, 5),
    version: int = 20,
) -> Dict[str, np.ndarray]:
    """Sample IRI and derive foF2/foE (MHz) and hmF2/hEs (km) for each timestamp."""
    heights = np.arange(alt_profile[0], alt_profile[1] + alt_profile[2], alt_profile[2], dtype=float)
    times = list(times)

    foF2, foE, hmF2, hEs = [], [], [], []
    for event in times:
        iri_output = iricore.iri(event, list(alt_profile), lat, lon, version)
        edens = np.asarray(iri_output.edens, dtype=float)
        freqs = utils.plasma_freq_mhz(edens)

        f2_freq, f2_height = _extract_peak(freqs, heights, heights >= 160)
        es_freq, es_height = _extract_peak(freqs, heights, (heights >= 90) & (heights <= 140))

        foF2.append(f2_freq if f2_freq > 0 else np.nan)
        hmF2.append(f2_height if f2_height > 0 else np.nan)
        foE.append(es_freq if es_freq > 0 else np.nan)
        hEs.append(es_height if es_height > 0 else np.nan)

    return {
        "foF2": np.array(foF2, dtype=float),
        "foE": np.array(foE, dtype=float),
        "hmF2": np.array(hmF2, dtype=float),
        "hEs": np.array(hEs, dtype=float),
    }


def prepare_profile_df(extractor: SaoExtractor)->pd.DataFrame:
    O = pd.DataFrame()
    for rec in extractor.sao.SAORecord:
        o = pd.DataFrame()
        if hasattr(rec.ProfileList, "Profile"):
            h, f = (
                np.array(rec.ProfileList.Profile[0].Tabulated.AltitudeList),
                np.array(rec.ProfileList.Profile[0].Tabulated.ProfileValueList[0].values)
            )
            arg_max = np.argmax(f)
            o["h"], o["f"] = h[:arg_max+1], f[:arg_max+1]
            o["datetime"] =  _to_utc_datetime(rec.StartTimeUTC)
            o.dropna(inplace=True)
            O = pd.concat([O, o])
    return O


def compute_binned_median_heights(profile: pd.DataFrame) -> pd.DataFrame:
    if profile.empty:
        return pd.DataFrame(columns=FREQ_BIN_CENTERS)

    def _median_per_bin(group: pd.DataFrame) -> np.ndarray:
        medians = []
        for low, high in zip(FREQ_BIN_EDGES[:-1], FREQ_BIN_EDGES[1:]):
            mask = (group["f"] >= low) & (group["f"] < high)
            if mask.any():
                medians.append(group.loc[mask, "h"].median())
            else:
                medians.append(np.nan)
        return np.array(medians)

    grouped = (
        profile.groupby("datetime", sort=True)
        .apply(_median_per_bin)
        .apply(pd.Series)
    )
    print(grouped)
    grouped.columns = FREQ_BIN_CENTERS
    grouped.sort_index(inplace=True)
    return grouped

def prepare_scaled_dataframe(extractor: SaoExtractor) -> np.ndarray:
    scaled = extractor.get_scaled_datasets_xml()
    datetimes = [_to_utc_datetime(rec.StartTimeUTC) for rec in extractor.sao.SAORecord]
    scaled["datetime"] = datetimes
    scaled.sort_values("datetime", inplace=True)
    return scaled

def discover_station_files(year: int, limit: int = 6) -> List[Path]:
    base_dir = Path(f"obs/")
    if not base_dir.exists():
        return []
    station_files = sorted(base_dir.glob("*_SAO.XML"))
    return station_files[:limit]


def plot_station_panel(ax, scaled, profile, iri_series, stn_info, time_limits):
    times = scaled["datetime"].to_list()

    binned_heights = compute_binned_median_heights(profile)
    cmap = cm.get_cmap("plasma", len(FREQ_BIN_CENTERS))
    if not binned_heights.empty:
        base = np.zeros(len(binned_heights.index))
        for idx, center in enumerate(FREQ_BIN_CENTERS):
            if center not in binned_heights.columns:
                continue
            heights = binned_heights[center].to_numpy(dtype=float)
            layer_times = binned_heights.index.to_list()
            mask = ~np.isnan(heights)
            if not np.any(mask):
                continue
            lower = base.copy()
            upper = base.copy()
            upper[mask] = heights[mask]
            ax.fill_between(
                layer_times,
                lower,
                upper,
                where=mask,
                color=cmap(idx),
                alpha=0.35,
                linewidth=0,
            )
            ax.plot(layer_times, heights, color=cmap(idx), lw=1.0)
            base[mask] = heights[mask]

    ax.scatter(times, scaled["hmF2"], s=14, marker="s", color=COLORS["hmF2"], ls="None", alpha=0.4, label="hmF2 obs")
    ax.plot(times, iri_series["hmF2"], color=COLORS["hmF2"], ls="-", lw=1.2, label="hmF2 IRI")
    
    ax.set_yscale("log")
    # ax.set_ylim(100, 500)
    ax.set_xlim(time_limits)
    ax.xaxis.set_major_locator(mdates.AutoDateLocator())
    ax.xaxis.set_major_formatter(DateFormatter("%H"))

    ax.set_ylabel("Height (km)", color=COLORS["hmF2"])

    title = f"{stn_info['URSI']}/{stn_info['STATIONNAME']}"
    ax.text(0.02, 1.02, title, transform=ax.transAxes, ha="left", va="bottom", fontsize=15, fontweight="bold")
    ax.axvline(dt.datetime(2025, 9, 4, 18), ls="--", lw=0.8, color="k")
    ax.axvline(dt.datetime(2025, 9, 5), ls="--", lw=0.8, color="k")
    ax.axvspan(dt.datetime(2025, 9, 4, 18), dt.datetime(2025, 9, 4, 20), alpha=0.2, color="k")

def main():
    target_date = dt.date(2025, 9, 2)
    time_limits = (
        dt.datetime.combine(target_date, dt.time.min),
        dt.datetime.combine(target_date + dt.timedelta(days=4), dt.time.min),
    )

    station_files = discover_station_files(target_date.year)
    if not station_files:
        raise FileNotFoundError(f"No SAO XML files found under data/{target_date.year}")

    fig, axes = plt.subplots(
        figsize=(8, 4 * 2),
        nrows=2, ncols=1,
        sharex=True
    )


    bin_handles = []

    for ax, sao_path in zip(axes, station_files):
        code = sao_path.name.split("_")[0]
        extractor = SaoExtractor(str(sao_path), extract_time_from_name=False, extract_stn_from_name=False)
        extractor.extract_xml()
        extractor.date = target_date
        extractor.local_time = target_date
        extractor.stn_info = get_digisonde_info(code.upper())

        for rec in extractor.sao.SAORecord:
            rec.StartTimeUTC = _to_utc_datetime(rec.StartTimeUTC)

        scaled = prepare_scaled_dataframe(extractor)
        profile = prepare_profile_df(extractor)
        iri_series = compute_iri_series(
            scaled["datetime"].tolist(),
            extractor.stn_info["LAT"],
            extractor.stn_info["LONG"],
        )

        plot_station_panel(ax, scaled, profile, iri_series, extractor.stn_info, time_limits)
        if ax.get_subplotspec().is_last_row():
            ax.set_xlabel("Time (UTC)")
        else:
            ax.set_xlabel("")
            ax.tick_params(labelbottom=False)

    cmap = cm.get_cmap("plasma", len(FREQ_BIN_CENTERS))
    for idx, center in enumerate(FREQ_BIN_CENTERS):
        bin_handles.append(
            Line2D(
                [0],
                [0],
                color=cmap(idx),
                lw=4,
                label=f"{center - 0.5:.0f}-{center + 0.5:.0f} MHz",
            )
        )

    leg = fig.legend(
        handles=bin_handles,
        loc="center left",
        ncol=1,
        frameon=False,
        fontsize=14,
        bbox_to_anchor=(0.9, 0.3),
        borderaxespad=0.0,
        handlelength=2.5,
    )

    fig.suptitle("Contour plots during 2-5 Spet 2025 SW/TID", fontsize=18, y=0.86, fontweight="bold")
    fig.tight_layout(rect=[0.03, 0.05, 0.92, 0.9])

    output_dir = Path("figures/2025")
    output_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_dir / "ionosonde_contours.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
