"""Skymap (DDM) tilt + Morse wavelet analysis for one event day, after Molina et al.

For one station/day, from the SKY files downloaded by ``didbase_pipeline.py --sky``:

1. Parse every SKY file (``pynasonde`` ``SkyExtractor``); cache the day as a pickle.
2. Filter echoes to the bulk F-region plasma: range window, zenith angle <= 25 deg
   and a per-sounding Doppler outlier filter (``[skymap]`` in config.toml).
3. Per sounding: echo count, zenith/azimuth histograms and the echo *centre*
   (mean x/y position) as the ionospheric tilt. ASSUMPTION: SKY ``x_coord`` is
   East-West (+East) and ``y_coord`` North-South (+North), in degrees from zenith.
4. Continuous Morse wavelet transform (gamma=3, beta=20) of the E-W and N-S tilt.
5. Overlay the ``mstid.py`` events for the day.

    python skymap.py AL945 2024-04-02
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger

import quality
import utils

ROOT = Path(__file__).resolve().parent
OUT_DIR = Path(__import__("os").environ.get("SKYMAP_OUT_DIR", str(ROOT / "Data" / "analysis")))
FIG_DIR = Path(__import__("os").environ.get("SKYMAP_FIG_DIR", str(ROOT / "Data" / "figures" / "skymap")))


def scfg() -> dict:
    return utils.load_config()["skymap"]


# ---------------------------------------------------------------------------
# 1. Load
# ---------------------------------------------------------------------------

def _parse_one(path: str):
    from pynasonde.digisonde.parsers.sky import SkyExtractor
    logger.disable("pynasonde")
    try:
        e = SkyExtractor(path, extract_time_from_name=True, extract_stn_from_name=True)
        e.extract()
        d = e.to_pandas()
        return d if not d.empty else None
    except Exception as exc:  # corrupt/empty file: skip, but say so
        logger.enable("pynasonde")
        logger.warning(f"{Path(path).name}: {exc}")
        return None


def load_day(station: str, day: pd.Timestamp, data_dir: Path | None = None, use_cache: bool = True) -> pd.DataFrame:
    """All echoes of one day (one row per source) with a ``file_time`` column."""
    from joblib import Parallel, delayed
    day = pd.Timestamp(day)
    data_dir = Path(data_dir or quality.DATA_DIR)
    cache = OUT_DIR / "sky" / f"{station}_{day:%Y%m%d}.pkl"
    if use_cache and cache.is_file():
        return pd.read_pickle(cache)
    files = sorted((data_dir / str(day.year) / station / "sky" / day.date().isoformat()).glob("*.SKY"))
    if not files:
        raise FileNotFoundError(f"No SKY files for {station} {day.date()} (run didbase_pipeline.py --sky)")
    parts = Parallel(n_jobs=8)(delayed(_parse_one)(str(f)) for f in files)
    df = pd.concat([p for p in parts if p is not None], ignore_index=True)
    df = df.rename(columns={"datetime": "file_time"})
    cache.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(cache)
    return df


# ---------------------------------------------------------------------------
# 2-3. Filter + tilt
# ---------------------------------------------------------------------------

def add_geometry(df: pd.DataFrame) -> pd.DataFrame:
    """zenith [deg] = radius of (x, y); azimuth [deg, clockwise from North] assuming x=East, y=North."""
    out = df.copy()
    out["zenith"] = np.hypot(out.x_coord, out.y_coord)
    out["azimuth"] = np.degrees(np.arctan2(out.x_coord, out.y_coord)) % 360.0
    return out


def filter_echoes(df: pd.DataFrame, c: dict | None = None) -> pd.DataFrame:
    """Range window, zenith limit and per-sounding Doppler outlier rejection."""
    c = c or scfg()
    d = add_geometry(df)
    d = d[d.group_range.between(c["range_min_km"], c["range_max_km"]) & (d.zenith <= c["max_zenith_deg"])]
    med = d.groupby("file_time").spect_dop.transform("median")
    mad = (d.spect_dop - med).abs().groupby(d.file_time).transform("median") * 1.4826
    return d[(d.spect_dop - med).abs() <= c["doppler_mad_k"] * mad.clip(lower=1.0)]


def tilt_series(echoes: pd.DataFrame, c: dict | None = None) -> pd.DataFrame:
    """Per-sounding echo count and centre (mean x = E-W, mean y = N-S) [deg]."""
    c = c or scfg()
    g = echoes.groupby("file_time")
    t = pd.DataFrame({"n_echoes": g.size(), "ew": g.x_coord.mean(), "ns": g.y_coord.mean(),
                      "mean_dop": g.spect_dop.mean(), "mean_range": g.group_range.mean()})
    t.loc[t.n_echoes < c["min_echoes"], ["ew", "ns", "mean_dop"]] = np.nan
    return t


# ---------------------------------------------------------------------------
# 4. Morse wavelet
# ---------------------------------------------------------------------------

def morse_cwt(x: np.ndarray, dt_min: float, periods: np.ndarray, gamma: float = 3.0, beta: float = 20.0):
    """Analytic generalised-Morse CWT in the frequency domain.

    Peak-normalised wavelet ``psi(w) = 2 (w/wp)^beta exp(beta/gamma (1 - (w/wp)^gamma))`` with
    ``wp = (beta/gamma)^(1/gamma)``. ``x`` must be uniformly sampled and NaN-free.
    Returns ``(W, coi_min)``: complex ``(n_periods, n)`` and the cone-of-influence half-width
    [min] per period (time-width of the wavelet envelope, ``sqrt(beta*gamma)/wp * s``).
    """
    n = len(x)
    x = np.asarray(x, float) - np.nanmean(x)
    xf = np.fft.fft(x)
    w = 2 * np.pi * np.fft.fftfreq(n, d=dt_min)  # rad/min
    wp = (beta / gamma) ** (1.0 / gamma)
    W = np.empty((len(periods), n), complex)
    coi = np.empty(len(periods))
    for i, P in enumerate(periods):
        s = P * wp / (2 * np.pi)  # scale [min] such that the peak frequency maps to period P
        ws = w * s
        psi = np.zeros(n)
        pos = ws > 0
        r = ws[pos] / wp
        psi[pos] = 2.0 * np.exp(beta * np.log(r) + (beta / gamma) * (1 - r ** gamma))
        W[i] = np.fft.ifft(xf * psi)
        coi[i] = np.sqrt(beta * gamma) / wp * s
    return W, coi


def wavelet_of(series: pd.Series, c: dict | None = None):
    """CWT of a tilt series on the uniform ``cadence_min`` grid (gaps linearly interpolated; short gaps only)."""
    c = c or scfg()
    s = series.dropna()
    grid = pd.date_range(s.index.min().ceil(f"{c['cadence_min']}min"), s.index.max(), freq=f"{c['cadence_min']}min")
    y = np.interp(grid.astype("int64"), s.index.astype("int64"), s.to_numpy())
    periods = np.geomspace(c["period_min"], c["period_max"], 60)
    W, coi = morse_cwt(y, c["cadence_min"], periods, c["morse_gamma"], c["morse_beta"])
    return grid, periods, np.abs(W) ** 2, coi, y


def ridge_in_window(grid, periods, power, coi, start, end):
    """Period of maximum mean power inside [start, end], ignoring points inside the cone of influence."""
    t = (grid - grid[0]).total_seconds().to_numpy() / 60.0
    sel = (grid >= start) & (grid <= end)
    if sel.sum() < 3:
        return np.nan, np.nan
    mask = np.ones_like(power, bool)
    for i, h in enumerate(coi):
        mask[i] = (t >= t[0] + h) & (t <= t[-1] - h)
    p = np.where(mask & sel[None, :], power, np.nan)
    with np.errstate(all="ignore"):
        mean_p = np.nanmean(p, axis=1)
    if np.all(np.isnan(mean_p)):
        return np.nan, np.nan
    i = int(np.nanargmax(mean_p))
    return float(periods[i]), float(mean_p[i])


# ---------------------------------------------------------------------------
# 5. Figure for one event day
# ---------------------------------------------------------------------------

def plot_event_day(station: str, day: pd.Timestamp, out_path: Path | None = None) -> dict:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    c = scfg()
    day = pd.Timestamp(day)
    raw = load_day(station, day)
    ech = filter_echoes(raw, c)
    tilt = tilt_series(ech, c)
    logger.info(f"{station} {day.date()}: {raw.file_time.nunique()} soundings, {len(raw)} echoes -> {len(ech)} kept; "
                f"{int(tilt.ew.notna().sum())} tilt samples")

    events = pd.DataFrame()
    ev_path = OUT_DIR / "mstid_events.csv"
    if ev_path.is_file():
        e = pd.read_csv(ev_path, parse_dates=["date", "start", "end"])
        events = e[(e.station == station) & (e.date == day)]

    fig, ax = plt.subplots(7, 1, figsize=(10, 17), sharex=True, constrained_layout=True)
    x0 = day
    # (0) foF2 against the quiet baseline, for context
    try:
        import mstid
        path = quality.DATA_DIR / str(day.year) / station / f"{day.date().isoformat()}.csv"
        seg, _ = quality.clean_day(quality.load_day(path, ["foF2"]), "foF2", quality.cfg())
        base = mstid.load_baseline(station, "foF2")
        b = mstid.baseline_at(base[base.month == day.month], seg["datetime"])
        ax[0].scatter(seg.datetime, seg.foF2, s=6, color="grey", alpha=0.6, lw=0)
        ax[0].fill_between(seg.datetime, b["q25"], b["q75"], color="red", alpha=0.2, lw=0)
        ax[0].plot(seg.datetime, b["median"], color="red", lw=1.4)
    except Exception as exc:
        logger.warning(f"foF2 context panel skipped: {exc}")
    ax[0].set_ylabel("foF2 [MHz]")
    # (1) echo counts
    ax[1].plot(tilt.index, tilt.n_echoes, color="tab:green")
    ax[1].plot(tilt.index, raw.groupby("file_time").size().reindex(tilt.index), color="grey", lw=0.8)
    ax[1].set_ylabel("echoes/sounding\n(green kept, grey all)")
    # (2,3) zenith / azimuth histograms vs time (Pradipta et al. 2015 style)
    for a, col, bins, lab in ((ax[2], "zenith", np.arange(0, c["max_zenith_deg"] + 0.5, 0.5), "zenith [deg]"),
                              (ax[3], "azimuth", np.arange(0, 365, 5), "azimuth [deg]")):
        times = np.sort(ech.file_time.unique())
        H = np.array([np.histogram(ech.loc[ech.file_time == t, col], bins=bins)[0] for t in times]).T
        a.pcolormesh(times, bins, H, shading="auto", cmap="viridis")
        a.set_ylabel(lab)
    # (4) tilt time series
    ax[4].plot(tilt.index, tilt.ew, color="tab:blue", lw=1.0, label="E-W (x)")
    ax[4].plot(tilt.index, tilt.ns, color="tab:orange", lw=1.0, label="N-S (y)")
    ax[4].axhline(0, color="k", lw=0.5)
    ax[4].set_ylabel("echo centre [deg]")
    ax[4].legend(fontsize=7, ncol=2)
    # (5,6) wavelets
    ridges = []
    for a, col in ((ax[5], "ew"), (ax[6], "ns")):
        grid, periods, power, coi, _ = wavelet_of(tilt[col], c)
        t = mdates.date2num(grid.to_pydatetime())
        a.contourf(t, periods, power, levels=20, cmap="turbo")
        a.set_yscale("log")
        a.set_ylabel(f"period [min]\n{'E-W' if col == 'ew' else 'N-S'} tilt")
        left = mdates.date2num((grid[0] + pd.to_timedelta(coi, "min")).to_pydatetime())
        right = mdates.date2num((grid[-1] - pd.to_timedelta(coi, "min")).to_pydatetime())
        a.plot(left, periods, "w--", lw=0.8)
        a.plot(right, periods, "w--", lw=0.8)
        for r in events.itertuples():
            P, pw = ridge_in_window(grid, periods, power, coi, r.start, r.end)
            ridges.append(dict(component=col, param=r.param, start=r.start, end=r.end,
                               mstid_period_min=r.period_min, wavelet_period_min=P, wavelet_power=pw))
    colours = {"foF2": "tab:blue", "hmF2": "tab:orange", "foE": "tab:green", "hE": "tab:red"}
    for a in ax:
        for r in events.itertuples():
            a.axvspan(r.start, r.end, color=colours.get(r.param, "k"), alpha=0.18, lw=0)
        a.grid(ls=":", alpha=0.4)
        a.tick_params(labelsize=8)
    ax[0].set_title(f"{station} {day.date()} skymap tilt + Morse wavelet  |  mstid events shaded: "
                    + ", ".join(f"{p}" for p in sorted(events.param.unique())) if len(events) else
                    f"{station} {day.date()} skymap tilt + Morse wavelet (no mstid events)", fontsize=9, loc="left")
    ax[-1].set_xlim(x0, x0 + pd.Timedelta(days=1))
    ax[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    ax[-1].xaxis.set_major_locator(mdates.HourLocator(interval=3))
    ax[-1].set_xlabel("UT   (white dashed: cone of influence)")
    out_path = out_path or FIG_DIR / station / f"skymap_{station}_{day:%Y%m%d}.png"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=110)
    plt.close(fig)
    logger.info(f"Saved: {out_path}")
    tilt.to_csv(OUT_DIR / f"skymap_tilt_{station}_{day:%Y%m%d}.csv")
    return dict(tilt=tilt, ridges=pd.DataFrame(ridges))


if __name__ == "__main__":
    res = plot_event_day(sys.argv[1], pd.Timestamp(sys.argv[2]))
    print(res["ridges"].to_string(index=False) if not res["ridges"].empty else "no events")
