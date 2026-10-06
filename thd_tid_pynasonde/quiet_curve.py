"""Statistical quiet-day curves from QC-passed DIDBASE parameters.

Pipeline (all thresholds in ``DEFAULTS`` below; override any of them with an
``[quiet_curve]`` table in config.toml):

1. Good data: station/day/parameter rows with ``passed == True`` in
   ``Data/quality/day_quality.csv`` (run ``quality.py`` first); samples come
   from the cleaned longest segment (``quality.clean_day``).
2. Quiet days: UT days where, on the day *and* the previous day, Kp, SYM-H and
   AE stay below storm levels. Kp comes from GFZ Potsdam (cached in
   ``Data/analysis/``); SYM-H/AE from OMNI via pyomnidata (local ``~/omnidata``).
   Days with too little OMNI coverage or missing Kp are *not* called quiet.
3. Curve: for each station x parameter x calendar month, every quiet QC-passed
   day is averaged into UT time-of-day bins (one value per day per bin, so each
   day has equal weight); the per-bin distribution across days gives the
   median (the quiet-day curve), quartiles and the number of days. Bins with
   fewer than ``min_days`` days are NaN.

Outputs under ``Data/analysis/``: ``quiet_days.csv``,
``quiet_curve_<station>_<param>.csv`` and figures in
``Data/figures/quiet_curve/``. Run from this directory:

    python quiet_curve.py                      # all stations
    ANALYSIS_STATIONS=AL945 python quiet_curve.py
    ANALYSIS_STATIONS=AL945 ANALYSIS_BY_YEAR=1 python quiet_curve.py   # one curve per month per year
"""
from __future__ import annotations

import os
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import quality
import utils

ROOT = Path(__file__).resolve().parent
OUT_DIR = Path(os.environ.get("ANALYSIS_OUT_DIR", str(ROOT / "Data" / "analysis")))
FIG_DIR = Path(os.environ.get("ANALYSIS_FIG_DIR", str(ROOT / "Data" / "figures" / "quiet_curve")))
KP_URL = "https://kp.gfz.de/app/files/Kp_ap_Ap_SN_F107_since_1932.txt"

DEFAULTS = dict(
    kp_max=3.0,            # max 3-hourly Kp on the day and the previous day (<=)
    symh_min=-30.0,        # minimum SYM-H [nT] on the day and the previous day (>)
    ae_max=500.0,          # maximum AE [nT] on the day and the previous day (<)
    min_omni_coverage=0.8,  # fraction of a day's OMNI samples that must be finite
    bin_min=15.0,          # UT time-of-day bin width [min]
    min_days=5,            # minimum quiet days contributing to a bin
    by_year=False,         # False: pool all years per calendar month
)


def acfg() -> dict:
    """DEFAULTS < config.toml [quiet_curve] < env (ANALYSIS_BY_YEAR=1, ANALYSIS_MIN_DAYS=n)."""
    c = {**DEFAULTS, **utils.load_config().get("quiet_curve", {})}
    if os.environ.get("ANALYSIS_BY_YEAR"):
        c["by_year"] = os.environ["ANALYSIS_BY_YEAR"].strip().lower() in {"1", "true", "yes", "y", "on"}
    if os.environ.get("ANALYSIS_MIN_DAYS"):
        c["min_days"] = int(os.environ["ANALYSIS_MIN_DAYS"])
    return c


# ---------------------------------------------------------------------------
# Step 2: quiet days
# ---------------------------------------------------------------------------

def load_kp(path: Path | None = None) -> pd.DataFrame:
    """Daily max Kp and Ap from the GFZ file (downloaded once, then cached)."""
    path = path or OUT_DIR / "Kp_ap_Ap_SN_F107_since_1932.txt"
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        logger.info(f"Downloading Kp from {KP_URL}")
        with urllib.request.urlopen(urllib.request.Request(KP_URL, headers={"User-Agent": "tid-analysis"}), timeout=120) as r:
            path.write_bytes(r.read())
    cols = ["yy", "mm", "dd", "days", "days_m", "bsr", "dB"] + [f"Kp{i}" for i in range(1, 9)] \
        + [f"ap{i}" for i in range(1, 9)] + ["Ap", "SN", "F107obs", "F107adj", "D"]
    df = pd.read_csv(path, comment="#", sep=r"\s+", header=None, names=cols)
    df["date"] = pd.to_datetime(dict(year=df.yy, month=df.mm, day=df.dd)).dt.normalize()
    kp = df[[f"Kp{i}" for i in range(1, 9)]].where(lambda x: x >= 0)  # -1 = missing
    return pd.DataFrame({"date": df.date, "kp_max": kp.max(axis=1), "ap_daily": df.Ap.where(df.Ap >= 0)}).set_index("date")


def load_omni_daily(years) -> pd.DataFrame:
    """Per-UT-day min SYM-H, max AE and finite-sample coverage (5-min OMNI)."""
    import pyomnidata
    frames = []
    for year in sorted(set(years)):
        d = pyomnidata.GetOMNI(year, Res=5)
        df = pd.DataFrame({"date": pd.to_datetime(d["Date"].astype(str), format="%Y%m%d"),
                           "SymH": d["SymH"].astype(float), "AE": d["AE"].astype(float)})
        g = df.groupby("date")
        frames.append(pd.DataFrame({"symh_min": g.SymH.min(), "ae_max": g.AE.max(),
                                    "omni_coverage": g.apply(lambda x: float((x.SymH.notna() & x.AE.notna()).mean()),
                                                            include_groups=False)}))
    return pd.concat(frames)


def quiet_day_table(years, c: dict | None = None) -> pd.DataFrame:
    """One row per UT day; ``quiet`` is True only if the day *and* the previous day are quiet."""
    c = c or acfg()
    first, last = pd.Timestamp(min(years) - 1, 12, 31), pd.Timestamp(max(years), 12, 31)
    days = pd.date_range(first, last)
    t = pd.DataFrame(index=days).join(load_kp()).join(load_omni_daily(list(range(min(years) - 1, max(years) + 1))))
    t.index.name = "date"
    t["day_ok"] = (t.kp_max <= c["kp_max"]) & (t.symh_min > c["symh_min"]) & (t.ae_max < c["ae_max"]) \
        & (t.omni_coverage >= c["min_omni_coverage"])
    t["day_ok"] = t.day_ok.fillna(False)  # unknown Kp/OMNI -> not quiet
    t["quiet"] = t.day_ok & t.day_ok.shift(1, fill_value=False)
    why = []
    for r, prev in zip(t.itertuples(), t.shift(1).itertuples()):
        reasons = []
        for tag, row in (("", r), ("prev ", prev)):
            if not (row.kp_max <= c["kp_max"]):
                reasons.append(f"{tag}Kp>{c['kp_max']:g} or missing")
            if not (row.symh_min > c["symh_min"]):
                reasons.append(f"{tag}SYM-H<={c['symh_min']:g} or missing")
            if not (row.ae_max < c["ae_max"]):
                reasons.append(f"{tag}AE>={c['ae_max']:g} or missing")
            if not (row.omni_coverage >= c["min_omni_coverage"]):
                reasons.append(f"{tag}OMNI coverage low")
        why.append("" if not reasons else ";".join(dict.fromkeys(reasons)))
    t["reason"] = why
    return t.loc[pd.Timestamp(min(years), 1, 1):]


# ---------------------------------------------------------------------------
# Step 3: quiet-day curve
# ---------------------------------------------------------------------------

def bin_day(seg: pd.DataFrame, param: str, bin_min: float) -> pd.Series:
    """Mean of one day's samples per UT time-of-day bin (index = bin start, hours)."""
    tod = seg["datetime"].dt.hour * 60 + seg["datetime"].dt.minute + seg["datetime"].dt.second / 60.0
    b = (tod // bin_min) * bin_min / 60.0
    return seg[param].groupby(b.to_numpy()).mean()


def quiet_curve(day_series: dict, c: dict | None = None) -> pd.DataFrame:
    """Per-bin distribution over days: median (the curve), quartiles, n_days.

    ``day_series`` maps a label (date) -> Series from ``bin_day``.
    """
    c = c or acfg()
    edges = np.arange(0, 24, c["bin_min"] / 60.0)
    mat = pd.DataFrame(day_series).reindex(edges)
    n = mat.notna().sum(axis=1)
    out = pd.DataFrame({"ut_hour": edges + c["bin_min"] / 120.0,  # bin centre
                        "median": mat.median(axis=1), "q25": mat.quantile(0.25, axis=1),
                        "q75": mat.quantile(0.75, axis=1), "n_days": n.to_numpy()})
    out.loc[out.n_days < c["min_days"], ["median", "q25", "q75"]] = np.nan
    return out.reset_index(drop=True)


def build_curves(stations=None, c: dict | None = None):
    """Yield (station, param, group, curve_df, day_series) for every group with data."""
    c = c or acfg()
    qc = quality.load_quality()
    qc = qc[qc.passed]
    if stations:
        qc = qc[qc.station.isin(stations)]
    qc = qc.assign(date=pd.to_datetime(qc.date))
    quiet = quiet_day_table(sorted(qc.date.dt.year.unique()), c)
    qc = qc[qc.date.map(quiet.quiet).fillna(False).astype(bool)]
    logger.info(f"{len(qc)} quiet QC-passed station-day-params; {int(quiet.quiet.sum())} quiet UT days overall")
    qcfg = quality.cfg()
    for (station, param), g in qc.groupby(["station", "param"]):
        series = {}
        for day, rows in g.groupby("date"):
            path = quality.DATA_DIR / str(day.year) / station / f"{day.date().isoformat()}.csv"
            seg, _ = quality.clean_day(quality.load_day(path, [param]), param, qcfg)
            if not seg.empty:
                series[day] = bin_day(seg, param, c["bin_min"])
        keys = sorted({(d.year if c["by_year"] else 0, d.month) for d in series})
        for year, month in keys:
            sel = {d: s for d, s in series.items() if d.month == month and (not c["by_year"] or d.year == year)}
            curve = quiet_curve(sel, c)
            yield station, param, (year, month), curve, sel


def plot_curves(station: str, param: str, groups: dict, c: dict, out_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    months = sorted(groups)
    ncol = 3
    fig, axes = plt.subplots(-(-len(months) // ncol), ncol, figsize=(12, 3.2 * -(-len(months) // ncol)),
                             sharex=True, constrained_layout=True, squeeze=False)
    units = "km" if param.lower().startswith("h") else "MHz"
    for ax, key in zip(axes.flat, months):
        curve, sel = groups[key]
        for s in sel.values():
            ax.plot(s.index + c["bin_min"] / 120.0, s.values, color="grey", alpha=0.25, lw=0.6)
        ax.fill_between(curve.ut_hour, curve.q25, curve.q75, color="red", alpha=0.2, lw=0)
        ax.plot(curve.ut_hour, curve["median"], color="red", lw=2)
        label = f"{pd.Timestamp(2000, key[1], 1):%B}" + (f" {key[0]}" if c["by_year"] else "")
        ax.set_title(f"{station} {param} - {label} ({len(sel)} quiet days)", fontsize=9, loc="left")
        ax.grid(ls=":", alpha=0.5)
        ax.set_xlim(0, 24)
    for ax in axes[:, 0]:
        ax.set_ylabel(f"{param} [{units}]")
    for ax in axes[-1]:
        ax.set_xlabel("UT [h]")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


def main():
    c = acfg()
    stations = [s.strip() for s in os.environ.get("ANALYSIS_STATIONS", "").split(",") if s.strip()] or None
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    years = sorted(pd.to_datetime(quality.load_quality().date).dt.year.unique())
    quiet = quiet_day_table(years, c)
    quiet.to_csv(OUT_DIR / "quiet_days.csv")
    logger.info(f"quiet_days.csv: {int(quiet.quiet.sum())}/{len(quiet)} days quiet")

    grouped: dict = {}
    for station, param, key, curve, sel in build_curves(stations, c):
        grouped.setdefault((station, param), {})[key] = (curve, sel)
    for (station, param), groups in grouped.items():
        table = pd.concat([cv.assign(year=k[0] or np.nan, month=k[1]) for k, (cv, _) in groups.items()])
        table.to_csv(OUT_DIR / f"quiet_curve_{station}_{param}.csv", index=False)
        plot_curves(station, param, groups, c, FIG_DIR / f"quiet_curve_{station}_{param}.png")


if __name__ == "__main__":
    main()
