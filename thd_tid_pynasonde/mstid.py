"""Per-day MSTID signatures in foF2/hmF2/foE/hE against the quiet-day baseline.

For every QC-passed station/day/parameter:

A. Residual: the day's cleaned samples minus the station's quiet-day median
   curve (``quiet_curve.py``) for that month, interpolated (24 h periodic) onto
   the sample times (NaN where the baseline is NaN). The residual is then
   de-spiked (Hampel, ``despike_*``) so isolated jumps/jitter cannot ring the filter.
B. MSTID band: ``quality.bandpass`` with the short-period edge raised to
   ``max(band_short_min, nyquist_margin * cadence)`` (so it stays well clear of the
   Nyquist period where scaler jitter lives) up to 60 min; output blanked inside
   gaps wider than ``bandpass.mask_gap_min``.
C. Metrics: band-pass RMS, peak-to-peak, ``amp90`` (90th percentile of |band-pass|,
   an amplitude that single spikes cannot dominate), dominant period
   (Lomb-Scargle), and wave-packet events: runs where the Hilbert envelope exceeds
   the null amplitude for >= ``min_event_min`` AND that are a coherent wave
   (Lomb-Scargle power >= ``min_ls_power`` at a period within the band, and
   duration >= ``min_cycles`` periods). Each event has start/end, amplitude,
   period, local solar time and a flag for being near the solar terminator.
D. Significance: the null is the ``null_percentile``-th percentile of ``amp90``
   over the *quiet* days of the same station x parameter x month. A day's
   ``amp_percentile`` is its rank within that distribution. ``enhanced`` =
   ``amp90`` above the null AND at least one non-terminator event.

Outputs (``Data/analysis/``): ``mstid_days.csv``, ``mstid_events.csv``, summary
figure per station, and per-day figures in ``Data/figures/mstid/<station>/``.

    ANALYSIS_STATIONS=AL945 python mstid.py
"""
from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import pandas as pd
from loguru import logger
from scipy.signal import hilbert, lombscargle

import quality
import quiet_curve
import utils

ROOT = Path(__file__).resolve().parent
OUT_DIR = quiet_curve.OUT_DIR
FIG_DIR = Path(os.environ.get("MSTID_FIG_DIR", str(ROOT / "Data" / "figures" / "mstid")))


def band_config(c: dict) -> dict:
    """Full config with the MSTID band's short-period limits swapped in for quality.bandpass."""
    full = dict(utils.load_config())
    full["bandpass"] = {**full["bandpass"], "period_min_short": c["band_short_min"],
                        "nyquist_margin": c["nyquist_margin"]}
    return full


def mcfg() -> dict:
    c = dict(utils.load_config()["mstid"])
    if os.environ.get("MSTID_PLOT"):
        c["plot"] = os.environ["MSTID_PLOT"].strip().lower()
    return c


# ---------------------------------------------------------------------------
# A. Baseline + residual
# ---------------------------------------------------------------------------

def load_baseline(station: str, param: str, out_dir: Path = OUT_DIR) -> pd.DataFrame | None:
    path = out_dir / f"quiet_curve_{station}_{param}.csv"
    return pd.read_csv(path) if path.is_file() else None


def baseline_at(curve: pd.DataFrame, times: pd.Series) -> dict:
    """Interpolate a monthly quiet curve to timestamps (24 h periodic).

    ``curve`` has ``ut_hour``, ``median``, ``q25``, ``q75`` for one month.
    A value is NaN when either bracketing bin is NaN. Returns arrays by name.
    """
    c = curve.sort_values("ut_hour").reset_index(drop=True)
    x = c.ut_hour.to_numpy()
    xe = np.concatenate([[x[-1] - 24.0], x, [x[0] + 24.0]])
    hours = (times.dt.hour + times.dt.minute / 60.0 + times.dt.second / 3600.0).to_numpy()
    right = np.searchsorted(xe, hours, side="right").clip(1, len(xe) - 1)
    left = right - 1
    out = {}
    for col in ("median", "q25", "q75"):
        ye = np.concatenate([[c[col].iloc[-1]], c[col].to_numpy(), [c[col].iloc[0]]])
        w = (hours - xe[left]) / (xe[right] - xe[left])
        out[col] = ye[left] * (1 - w) + ye[right] * w  # NaN propagates if either neighbour is NaN
    return out


# ---------------------------------------------------------------------------
# C. Events
# ---------------------------------------------------------------------------

def dominant_period(t_min: np.ndarray, y: np.ndarray, band: tuple[float, float]):
    """Lomb-Scargle peak period [min] and its normalised power within ``band``."""
    if len(y) < 8:
        return np.nan, np.nan
    periods = np.linspace(band[0], band[1], 120)
    p = lombscargle(t_min, y - y.mean(), 2 * np.pi / periods, normalize=True)
    i = int(np.argmax(p))
    return float(periods[i]), float(p[i])


def mean_period(t_min: np.ndarray, y: np.ndarray) -> float:
    """Period [min] from zero-crossing spacing (2 x mean half-period); NaN if < 2 crossings."""
    s = np.sign(y)
    idx = np.flatnonzero(s[:-1] * s[1:] < 0)
    if len(idx) < 2:
        return np.nan
    tz = t_min[idx] - y[idx] * (t_min[idx + 1] - t_min[idx]) / (y[idx + 1] - y[idx])
    return float(2 * np.mean(np.diff(tz)))


def detect_events(filtered: pd.Series, amp_thr: float, cadence: float, c: dict,
                  band: tuple[float, float] | None = None) -> list[dict]:
    """Coherent wave packets: runs of the Hilbert envelope above ``amp_thr``.

    ``filtered`` is the band-passed series on a uniform grid with NaN in masked
    gaps; each NaN-free chunk is analysed separately and its edges (``edge_min``)
    are ignored. A run is kept only if it lasts >= ``min_event_min``, its
    Lomb-Scargle peak (within ``band``) has normalised power >= ``min_ls_power``,
    and it spans >= ``min_cycles`` periods -- this rejects one-sample jitter and
    isolated jumps, which have no coherent period inside the band.
    """
    events = []
    valid = filtered.notna().to_numpy()
    if not valid.any() or not np.isfinite(amp_thr):
        return events
    band = band or (15.0, 60.0)
    cuts = np.flatnonzero(np.diff(valid.astype(int)) != 0) + 1
    for chunk in np.split(np.arange(len(filtered)), cuts):
        if not valid[chunk[0]] or len(chunk) < 8:
            continue
        y = filtered.iloc[chunk].to_numpy()
        env = np.abs(hilbert(y))
        trim = int(round(c["edge_min"] / cadence))
        if trim > 0:
            env[:trim] = np.nan
            env[-trim:] = np.nan
        above = np.nan_to_num(env, nan=0.0) > amp_thr
        runs = np.flatnonzero(np.diff(np.concatenate([[0], above.astype(int), [0]])) != 0).reshape(-1, 2)
        for a, b in runs:  # [a, b)
            dur = (b - a) * cadence
            if dur < c["min_event_min"]:
                continue
            t = filtered.index[chunk[a:b]]
            tmin = (t - t[0]).total_seconds().to_numpy() / 60.0
            period, power = dominant_period(tmin, y[a:b], band)
            if not (power >= c["min_ls_power"]) or dur < c["min_cycles"] * period:
                continue
            events.append(dict(start=t[0], end=t[-1], duration_min=dur,
                               amplitude=float(np.nanmax(env[a:b])), period_min=period,
                               ls_power=power, n_cycles=dur / period))
    return events


def near_terminator(start, end, lat: float, lon: float, margin_min: float) -> bool:
    """True if SZA crosses 90 deg within [start - margin, end + margin]."""
    t = pd.date_range(start - pd.Timedelta(minutes=margin_min), end + pd.Timedelta(minutes=margin_min), freq="10min")
    sza = np.array([utils._solar_zenith_angle(x.to_pydatetime(), lat, lon) for x in t])
    return bool(sza.min() < 90.0 <= sza.max())


def local_hour(ts: pd.Timestamp, lon: float) -> float:
    return float((ts.hour + ts.minute / 60.0 + lon / 15.0) % 24.0)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def station_coords(station: str) -> tuple[float, float]:
    path = quality.DATA_DIR / "conus_stations.csv"
    row = pd.read_csv(path).set_index("URSI").loc[station]
    return float(row.LAT), float(row.LONG)


def analyse_station_param(station: str, param: str, qc: pd.DataFrame, quiet: pd.Series, c: dict):
    """Returns (days_df, events_df, per-day arrays for plotting)."""
    base = load_baseline(station, param)
    if base is None:
        logger.warning(f"{station} {param}: no quiet curve (run quiet_curve.py); skipping")
        return pd.DataFrame(), pd.DataFrame(), {}
    qcfg = quality.cfg()
    band_cfg = band_config(c)
    lat, lon = station_coords(station)
    rows, store = [], {}
    for day in sorted(qc.date):
        month_curve = base[(base.month == day.month) & (base.year.isna() | (base.year == day.year))]
        if month_curve.empty:
            continue
        path = quality.DATA_DIR / str(day.year) / station / f"{day.date().isoformat()}.csv"
        seg, m = quality.clean_day(quality.load_day(path, [param]), param, qcfg)
        if seg.empty:
            continue
        b = baseline_at(month_curve, seg["datetime"])
        resid = seg.assign(**{param: seg[param].to_numpy() - b["median"]}).dropna(subset=[param])
        if len(resid) < 8:
            continue
        window = max(3, int(round(c["despike_window_min"] / m["cadence_min"])) | 1)
        spike = quality.hampel_mask(resid[param].to_numpy(), window, c["despike_k"], qcfg["min_scale"]["km" if param.lower().startswith("h") else "MHz"])
        n_despiked = int(spike.sum())
        resid = resid[~spike].reset_index(drop=True)
        filt, band = quality.bandpass(resid, param, m["cadence_min"], band_cfg)
        if filt is None:
            continue
        v = filt.dropna()
        if len(v) < 8:
            continue
        tmin = (v.index - v.index[0]).total_seconds().to_numpy() / 60.0
        period, power = dominant_period(tmin, v.to_numpy(), band)
        med_base = float(np.nanmedian(b["median"]))
        rms = float(np.sqrt(np.mean(v.to_numpy() ** 2)))
        amp90 = float(np.percentile(np.abs(v.to_numpy()), 90))
        rows.append(dict(station=station, param=param, date=day, month=day.month, quiet=bool(quiet.get(day, False)),
                         cadence_min=m["cadence_min"], band_short_min=band[0], band_long_min=band[1],
                         rms=rms, amp90=amp90, n_despiked=n_despiked, p2p=float(v.max() - v.min()), rms_pct_of_baseline=100 * rms / med_base, amp90_pct_of_baseline=100 * amp90 / med_base,
                         dominant_period_min=period, ls_power=power,
                         resid_mean=float(resid[param].mean()), n_baseline_days=int(month_curve.n_days.max())))
        store[day] = dict(seg=seg, base=b, resid=resid, filt=filt, band=band, cadence=m["cadence_min"])
    days = pd.DataFrame(rows)
    if days.empty:
        return days, pd.DataFrame(), store

    # D. significance against quiet days of the same month
    thr, pct, nnull = {}, [], {}
    for month, g in days.groupby("month"):
        null = g.loc[g.quiet, "amp90"].to_numpy()
        nnull[month] = len(null)
        thr[month] = np.percentile(null, c["null_percentile"]) if len(null) >= c["min_null_days"] else np.nan
    for r in days.itertuples():
        null = days.loc[(days.month == r.month) & days.quiet, "amp90"].to_numpy()
        pct.append(float((null <= r.amp90).mean()) if len(null) >= c["min_null_days"] else np.nan)
    days["amp_percentile"] = pct
    days["null_amp_thr"] = days.month.map(thr)
    days["n_null_days"] = days.month.map(nnull)

    events = []
    for r in days.itertuples():
        amp_thr = r.null_amp_thr  # amp90 ~ sinusoid amplitude, so compare the envelope to it directly
        for ev in detect_events(store[r.date]["filt"], amp_thr, r.cadence_min, c, band=store[r.date]["band"]):
            mid = ev["start"] + (ev["end"] - ev["start"]) / 2
            events.append(dict(station=station, param=param, date=r.date, quiet=r.quiet, **ev,
                               local_hour=local_hour(mid, lon), sza_deg=utils._solar_zenith_angle(mid.to_pydatetime(), lat, lon),
                               near_terminator=near_terminator(ev["start"], ev["end"], lat, lon, c["terminator_margin_min"])))
        store[r.date]["amp_thr"] = amp_thr
    events = pd.DataFrame(events)
    nev = events.groupby("date").size() if not events.empty else pd.Series(dtype=int)
    nev_clean = events[~events.near_terminator].groupby("date").size() if not events.empty else pd.Series(dtype=int)
    days["n_events"] = days.date.map(nev).fillna(0).astype(int)
    days["n_events_non_terminator"] = days.date.map(nev_clean).fillna(0).astype(int)
    days["enhanced"] = (days.amp90 > days.null_amp_thr) & (days.n_events_non_terminator > 0)
    return days, events, store


def plot_day(station: str, param: str, row, events: pd.DataFrame, st: dict, out_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt

    day = row.date
    fig, axes = plt.subplots(3, 1, figsize=(7, 7), sharex=True, constrained_layout=True)
    units = "km" if param.lower().startswith("h") else "MHz"
    seg, b = st["seg"], st["base"]
    t = seg["datetime"]
    ax = axes[0]
    ax.scatter(t, seg[param], s=6, color="grey", alpha=0.6, lw=0, label="data (QC segment)")
    ax.fill_between(t, b["q25"], b["q75"], color="red", alpha=0.2, lw=0, label="quiet IQR")
    ax.plot(t, b["median"], color="red", lw=1.6, label="quiet median")
    ax.set_ylabel(f"{param} [{units}]")
    ax.legend(fontsize=7, loc="best", framealpha=0.7, edgecolor="none")
    ax.set_title(f"{station} {param} {day.date()} | amp90 {row.amp90:.3g} ({row.amp_percentile:.0%} of quiet) | "
                 f"P={row.dominant_period_min:.0f} min | {'ENHANCED' if row.enhanced else 'not enhanced'}",
                 fontsize=8, loc="left")
    axes[1].plot(st["resid"]["datetime"], st["resid"][param], color="tab:blue", lw=0.8)
    axes[1].axhline(0, color="k", lw=0.5)
    axes[1].set_ylabel(f"data - quiet [{units}]")
    f = st["filt"]
    axes[2].plot(f.index, f.values, color="tab:red", lw=1.0)
    if np.isfinite(st["amp_thr"]):
        axes[2].axhline(st["amp_thr"], color="k", ls="--", lw=0.6)
        axes[2].axhline(-st["amp_thr"], color="k", ls="--", lw=0.6)
    axes[2].axhline(0, color="k", lw=0.5)
    for ev in events.itertuples():
        axes[2].axvspan(ev.start, ev.end, color="tab:orange" if ev.near_terminator else "tab:green", alpha=0.25, lw=0)
    axes[2].set_ylabel(f"band-pass {st['band'][0]:g}-{st['band'][1]:g} min [{units}]", fontsize=8)
    x0 = pd.Timestamp(day)
    axes[2].set_xlim(x0, x0 + pd.Timedelta(days=1))
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    axes[2].xaxis.set_major_locator(mdates.HourLocator(interval=3))
    axes[2].set_xlabel("UT   (green: event, orange: event near terminator)")
    for a in axes:
        a.grid(ls=":", alpha=0.5)
        a.tick_params(labelsize=8)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=120)
    plt.close(fig)


def plot_summary(station: str, days: pd.DataFrame, events: pd.DataFrame, quiet_tbl: pd.DataFrame, out_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(14, 3.8), constrained_layout=True)
    params = list(quality.cfg()["parameters"])
    for i, p in enumerate(params):
        d = days[days.param == p]
        frac = d.groupby("month").enhanced.mean()
        ax[0].bar(frac.index + (i - 1.5) * 0.2, frac.values, width=0.2, label=p)
    ax[0].set_xlabel("month")
    ax[0].set_ylabel("fraction of days enhanced")
    ax[0].legend(fontsize=7)
    if not events.empty:
        e = events[~events.near_terminator]
        for p in params:
            ax[1].hist(e.loc[e.param == p, "local_hour"], bins=np.arange(0, 25, 1), histtype="step", lw=1.5, label=p)
        ax[1].legend(fontsize=7)
    ax[1].set_xlabel("local solar time of event [h] (non-terminator)")
    ax[1].set_ylabel("events")
    for p in params:
        d = days[days.param == p].set_index("date")
        kp = quiet_tbl.kp_max.reindex(d.index)
        ax[2].scatter(kp, d.amp90_pct_of_baseline, s=8, alpha=0.5, label=p)
    ax[2].set_xlabel("max Kp of the day")
    ax[2].set_ylabel("band-pass amp90 [% of baseline]")
    ax[2].legend(fontsize=7)
    for a in ax:
        a.grid(ls=":", alpha=0.5)
    fig.suptitle(f"{station} MSTID summary", fontsize=10)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130)
    plt.close(fig)
    logger.info(f"Saved: {out_path}")


def main():
    c = mcfg()
    stations = [s.strip() for s in os.environ.get("ANALYSIS_STATIONS", "").split(",") if s.strip()] or None
    qc = quality.load_quality()
    qc = qc[qc.passed].assign(date=lambda d: pd.to_datetime(d.date))
    if stations:
        qc = qc[qc.station.isin(stations)]
    quiet_tbl = quiet_curve.quiet_day_table(sorted(qc.date.dt.year.unique()))
    all_days, all_events = [], []
    for station in sorted(qc.station.unique()):
        s_days, s_events = [], []
        for param in quality.cfg()["parameters"]:
            sub = qc[(qc.station == station) & (qc.param == param)]
            days, events, store = analyse_station_param(station, param, sub, quiet_tbl.quiet, c)
            if days.empty:
                continue
            logger.info(f"{station} {param}: {len(days)} days, {int(days.enhanced.sum())} enhanced, {len(events)} events")
            s_days.append(days)
            s_events.append(events)
            if c["plot"] != "none":
                for r in days.itertuples():
                    if c["plot"] == "all" or r.enhanced:
                        ev = events[events.date == r.date] if not events.empty else events
                        plot_day(station, param, r, ev, store[r.date],
                                 FIG_DIR / station / f"{param}_mstid_{station}_{r.date:%Y%m%d}.png")
        if s_days:
            d = pd.concat(s_days)
            e = pd.concat([x for x in s_events if not x.empty]) if any(not x.empty for x in s_events) else pd.DataFrame()
            plot_summary(station, d, e, quiet_tbl, FIG_DIR / f"mstid_summary_{station}.png")
            all_days.append(d)
            all_events.append(e)
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    if all_days:
        pd.concat(all_days).to_csv(OUT_DIR / "mstid_days.csv", index=False)
    ev = [e for e in all_events if not e.empty]
    pd.concat(ev).to_csv(OUT_DIR / "mstid_events.csv", index=False) if ev else None
    logger.info(f"Wrote mstid_days.csv / mstid_events.csv to {OUT_DIR}")


if __name__ == "__main__":
    main()
