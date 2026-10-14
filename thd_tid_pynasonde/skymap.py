"""Skymap (DDM) tilt + Morse wavelet analysis for one event day, after Molina et al.

For one station/day, from the SKY files downloaded by ``didbase_pipeline.py --sky``:

1. Parse every SKY file (``pynasonde`` ``SkyExtractor``); cache the day as a pickle.
2. Filter echoes to the bulk F-region plasma: range window, zenith angle <= 25 deg
   and a per-sounding Doppler outlier filter (``[skymap]`` in config.toml).
3. Per sounding: echo count, zenith/azimuth histograms and the echo *centre*
   (median position) as the ionospheric tilt. Axis convention (that of
   ``pynasonde.digisonde.digi_plots.convert_to_rt``): SKY ``x_coord`` is positive
   toward NORTH and ``y_coord`` positive toward WEST (right-handed with z up), in
   degrees from zenith; so North = x, East = -y, azimuth = atan2(-y, x) clockwise from N.
4. Continuous Morse wavelet transform (gamma=3, beta=20) of the E-W and N-S tilt,
   with a 95% red-noise (AR(1) Monte-Carlo) significance contour.
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
    """zenith [deg] = radius of (x, y); azimuth [deg, clockwise from North]; x = North, y = West."""
    out = df.copy()
    out["north"] = out.x_coord
    out["east"] = -out.y_coord
    out["zenith"] = np.hypot(out.x_coord, out.y_coord)
    out["azimuth"] = np.degrees(np.arctan2(out.east, out.north)) % 360.0
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
    """Per-sounding echo count and centre (E-W = -y positive East, N-S = x positive North) [deg], de-spiked."""
    c = c or scfg()
    g = echoes.groupby("file_time")
    t = pd.DataFrame({"n_echoes": g.size(), "ew": -g.y_coord.median(), "ns": g.x_coord.median(),
                      "mean_dop": g.spect_dop.mean(), "mean_range": g.group_range.mean()})
    t.loc[t.n_echoes < c["min_echoes"], ["ew", "ns", "mean_dop"]] = np.nan
    for col in ("ew", "ns"):  # drop isolated bad soundings (Hampel) so one outlier cannot dominate the wavelet
        v = t[col].to_numpy()
        ok = ~np.isnan(v)
        if ok.sum() > 7:
            sp = quality.hampel_mask(v[ok], c["despike_window"], c["despike_k"], c["despike_min_scale_deg"])
            idx = np.flatnonzero(ok)[sp]
            t.iloc[idx, t.columns.get_loc(col)] = np.nan
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


def red_noise_threshold(y: np.ndarray, c: dict | None = None, quantile: float = 95.0, n_surr: int = 200, seed: int = 0):
    """Per-period power threshold from AR(1) red-noise surrogates (Monte Carlo, cf. Torrence & Compo 1998).

    The lag-1 autocorrelation and variance of ``y`` set the surrogates; the threshold at each period is
    the ``quantile``-th percentile of the surrogate wavelet power pooled over time and realisations.
    """
    c = c or scfg()
    y = np.asarray(y, float) - np.mean(y)
    alpha = float(np.clip(np.corrcoef(y[:-1], y[1:])[0, 1], 0.0, 0.99)) if len(y) > 3 else 0.0
    sigma = float(np.std(y))
    rng = np.random.default_rng(seed)
    periods = np.geomspace(c["period_min"], c["period_max"], 60)
    pool = []
    for _ in range(n_surr):
        e = rng.normal(0.0, sigma * np.sqrt(1 - alpha ** 2), len(y))
        x = np.empty(len(y))
        x[0] = rng.normal(0.0, sigma)
        for i in range(1, len(y)):
            x[i] = alpha * x[i - 1] + e[i]
        W, _ = morse_cwt(x, c["cadence_min"], periods, c["morse_gamma"], c["morse_beta"])
        pool.append(np.abs(W) ** 2)
    return np.percentile(np.concatenate(pool, axis=1), quantile, axis=1), alpha


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
# 5. Sunrise/sunset, Lomb-Scargle
# ---------------------------------------------------------------------------

def solar_terminators(day: pd.Timestamp, lat: float, lon: float) -> list[tuple[str, pd.Timestamp]]:
    """Sunrise/sunset UT times on ``day`` (SZA crossing 90 deg, 1-min grid; formula good to ~1 deg)."""
    t = pd.date_range(day, day + pd.Timedelta(days=1), freq="1min")
    sza = np.array([utils._solar_zenith_angle(x.to_pydatetime(), lat, lon) for x in t])
    out = []
    for i in np.flatnonzero((sza[:-1] - 90.0) * (sza[1:] - 90.0) < 0):
        out.append(("sunrise" if sza[i + 1] < sza[i] else "sunset", t[i]))
    return out


def lomb_scargle(series: pd.Series, periods: np.ndarray):
    """Normalised Lomb-Scargle power of an irregularly sampled tilt series vs period [min]."""
    from scipy.signal import lombscargle
    s = series.dropna()
    if len(s) < 8:
        return np.full(len(periods), np.nan)
    t = (s.index - s.index[0]).total_seconds().to_numpy() / 60.0
    return lombscargle(t, s.to_numpy() - s.mean(), 2 * np.pi / periods, normalize=True)


# ---------------------------------------------------------------------------
# 6. Per-event folder: figures/skymap/<station>/<date>/
# ---------------------------------------------------------------------------

EVENT_COLOURS = {"foF2": "tab:blue", "hmF2": "tab:orange", "foE": "tab:green", "hE": "tab:red"}


def _decorate(axes, day, events, terms, xlabel=True):
    import matplotlib.dates as mdates
    for a in np.atleast_1d(axes):
        for r in events.itertuples():
            a.axvspan(r.start, r.end, color=EVENT_COLOURS.get(r.param, "k"), alpha=0.15, lw=0)
        for lab, t in terms:
            a.axvspan(t - pd.Timedelta(minutes=30), t + pd.Timedelta(minutes=30), color="gold", alpha=0.3, lw=0)
            a.axvline(t, color="goldenrod", lw=1.0, ls="--")
        a.grid(ls=":", alpha=0.4)
        a.tick_params(labelsize=8)
        a.set_xlim(day, day + pd.Timedelta(days=1))
    last = np.atleast_1d(axes)[-1]
    last.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    last.xaxis.set_major_locator(mdates.HourLocator(interval=3))
    if xlabel:
        last.set_xlabel("UT   (gold: sunrise/sunset +/- 30 min; shaded colours: mstid events)")


def _event_table(station, day):
    path = OUT_DIR / "mstid_events.csv"
    if not path.is_file():
        return pd.DataFrame(columns=["param", "start", "end", "period_min"])
    e = pd.read_csv(path, parse_dates=["date", "start", "end"])
    return e[(e.station == station) & (e.date == day)]


def plot_ionosonde(station, day, events, terms, out):
    """foF2 and hmF2 (grey) against the quiet-day median/IQR (red)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import mstid
    fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True, constrained_layout=True)
    path = quality.DATA_DIR / str(day.year) / station / f"{day.date().isoformat()}.csv"
    raw = quality.load_day(path, ["foF2", "hmF2"])
    for a, p, u in ((ax[0], "foF2", "MHz"), (ax[1], "hmF2", "km")):
        a.scatter(raw.datetime, raw[p], s=7, color="grey", alpha=0.6, lw=0, label="data")
        base = mstid.load_baseline(station, p)
        if base is not None and not base[base.month == day.month].empty:
            b = mstid.baseline_at(base[base.month == day.month], raw["datetime"])
            a.fill_between(raw.datetime, b["q25"], b["q75"], color="red", alpha=0.2, lw=0)
            a.plot(raw.datetime, b["median"], color="red", lw=1.5, label="quiet median / IQR")
        a.set_ylabel(f"{p} [{u}]")
    ax[0].legend(fontsize=7, loc="best")
    ax[0].set_title(f"{station} {day.date()}  foF2 / hmF2   events: "
                    + (", ".join(f"{r.param} {r.start:%H:%M}-{r.end:%H:%M}" for r in events.itertuples()) or "none"),
                    fontsize=8, loc="left")
    _decorate(ax, day, events, terms)
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_skymaps(ech, station, day, out_dir, per_page=9):
    """3x3 polar skymaps (North up, azimuth clockwise, radius = zenith deg); + Doppler cross, - circle."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    times = np.sort(ech.file_time.unique())
    groups = {t: g for t, g in ech.groupby("file_time")}
    vmax = float(np.nanpercentile(np.abs(ech.spect_dop), 98)) or 1.0
    n_pages = -(-len(times) // per_page)
    for page in range(n_pages):
        fig, axes = plt.subplots(3, 3, figsize=(11, 11), subplot_kw=dict(projection="polar"), constrained_layout=True)
        for a in axes.flat:
            a.set_visible(False)
        for a, t in zip(axes.flat, times[page * per_page:(page + 1) * per_page]):
            g = groups[t]
            a.set_visible(True)
            a.set_theta_zero_location("N")
            a.set_theta_direction(-1)
            a.set_rlim(0, 25)
            a.set_rticks([5, 10, 15, 20, 25])
            a.tick_params(labelsize=6)
            th = np.radians(g.azimuth.to_numpy())
            for mask, marker in ((g.spect_dop >= 0, "+"), (g.spect_dop < 0, "o")):
                sc = a.scatter(th[mask.to_numpy()], g.zenith[mask], c=g.spect_dop[mask], cmap="RdBu", vmin=-vmax,
                               vmax=vmax, s=9 if marker == "+" else 7, marker=marker, linewidths=0.6)
            a.set_title(f"{pd.Timestamp(t):%H:%M:%S}  n={len(g)}", fontsize=8)
        cb = fig.colorbar(sc, ax=axes, shrink=0.5, pad=0.02)
        cb.set_label("Doppler line number (blue + cross, red - circle)", fontsize=8)
        fig.suptitle(f"{station} {day.date()}  skymaps (filtered echoes)  page {page + 1}/{n_pages}", fontsize=10)
        fig.savefig(out_dir / f"skymaps_{page + 1:02d}.png", dpi=90)
        plt.close(fig)
    return n_pages


def plot_zenith_azimuth(ech, tilt, station, day, events, terms, out):
    """Echo-count histograms of zenith and azimuth vs time (all soundings) + echo-centre lines."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    times = np.sort(ech.file_time.unique())
    fig, ax = plt.subplots(2, 1, figsize=(11, 7), sharex=True, constrained_layout=True)
    c = scfg()
    for a, col, bins, lab in ((ax[0], "zenith", np.arange(0, c["max_zenith_deg"] + 0.5, 0.5), "zenith [deg]"),
                              (ax[1], "azimuth", np.arange(0, 365, 5), "azimuth [deg from N, cw]")):
        H = np.array([np.histogram(ech.loc[ech.file_time == t, col], bins=bins)[0] for t in times]).T
        m = a.pcolormesh(times, 0.5 * (bins[:-1] + bins[1:]), H, shading="nearest", cmap="viridis")
        a.set_ylabel(lab)
        fig.colorbar(m, ax=a, pad=0.01, label="echoes")
    cz = np.hypot(tilt.ew, tilt.ns)
    ca = np.degrees(np.arctan2(tilt.ew, tilt.ns)) % 360.0
    ax[0].plot(tilt.index, cz, "w.", ms=2.5)
    ax[1].plot(tilt.index, ca, "w.", ms=2.5)
    ax[0].set_title(f"{station} {day.date()}  echo zenith / azimuth distribution (white: echo centre)", fontsize=9, loc="left")
    _decorate(ax, day, events, terms)
    fig.savefig(out, dpi=110)
    plt.close(fig)


def plot_periods(tilt, station, day, events, terms, out, c):
    """Morse-wavelet scalograms of E-W / N-S tilt (+ Lomb-Scargle periodogram on the right)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.dates as mdates
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(12, 7), sharey="row", constrained_layout=True,
                             gridspec_kw=dict(width_ratios=[4, 1]))
    ridges = []
    for row, col, lab in ((0, "ew", "E-W tilt"), (1, "ns", "N-S tilt")):
        a, b = axes[row]
        n_ok = int(tilt[col].notna().sum())
        if n_ok < c["min_tilt_samples"]:  # sparse skymaps (few echoes/sounding): no meaningful wavelet
            for ax_ in (a, b):
                ax_.text(0.5, 0.5, f"{lab}: only {n_ok} tilt samples (< {c['min_tilt_samples']}; "
                         f"min_echoes={c['min_echoes']}) - wavelet skipped", ha="center", va="center",
                         transform=ax_.transAxes, fontsize=9)
            _decorate([a], day, events, terms, xlabel=(row == 1))
            ridges.append(dict(component=col, param="skipped: too few tilt samples", start=pd.NaT, end=pd.NaT,
                               mstid_period_min=np.nan, wavelet_period_min=np.nan, wavelet_power=np.nan, n_tilt_samples=n_ok))
            continue
        grid, periods, power, coi, yy = wavelet_of(tilt[col], c)
        thr, alpha = red_noise_threshold(yy, c)
        a.contourf(mdates.date2num(grid.to_pydatetime()), periods, power, levels=20, cmap="turbo")
        a.contour(mdates.date2num(grid.to_pydatetime()), periods, power / thr[:, None], levels=[1.0], colors="w", linewidths=1.6)
        a.set_yscale("log")
        a.set_ylabel(f"period [min]\n{lab}  (AR1 a={alpha:.2f})")
        a.plot(mdates.date2num((grid[0] + pd.to_timedelta(coi, "min")).to_pydatetime()), periods, "w--", lw=0.8)
        a.plot(mdates.date2num((grid[-1] - pd.to_timedelta(coi, "min")).to_pydatetime()), periods, "w--", lw=0.8)
        for lo, hi in ((15, 60),):
            a.axhline(lo, color="w", lw=0.5, ls=":")
            a.axhline(hi, color="w", lw=0.5, ls=":")
        _decorate([a], day, events, terms, xlabel=(row == 1))
        if row == 0:
            a.set_title(f"{station} {day.date()}  Morse wavelet (gamma={c['morse_gamma']:g}, beta={c['morse_beta']:g}) | "
                        f"Lomb-Scargle   (white contour: 95% red-noise significance; dotted: 15-60 min band)", fontsize=8, loc="left")
        ls = lomb_scargle(tilt[col], periods)
        b.plot(ls, periods, color="k")
        b.set_xlabel("LS power")
        b.grid(ls=":", alpha=0.4)
        for r in events.itertuples():
            P, pw = ridge_in_window(grid, periods, power, coi, r.start, r.end)
            ridges.append(dict(component=col, param=r.param, start=r.start, end=r.end,
                               mstid_period_min=r.period_min, wavelet_period_min=P, wavelet_power=pw,
                               signif_ratio=(pw / float(np.interp(P, periods, thr))) if np.isfinite(P) else np.nan))
        ridges.append(dict(component=col, param="whole-day LS peak", start=grid[0], end=grid[-1],
                           mstid_period_min=np.nan, wavelet_period_min=float(periods[int(np.nanargmax(ls))]) if np.isfinite(ls).any() else np.nan,
                           wavelet_power=float(np.nanmax(ls)) if np.isfinite(ls).any() else np.nan))
    fig.savefig(out, dpi=110)
    plt.close(fig)
    return pd.DataFrame(ridges)


def analyse_event_day(station: str, day, base_dir: Path | None = None) -> Path:
    """Write the full figure set for one event day into ``figures/skymap/<station>/<date>/``."""
    import mstid
    c = scfg()
    day = pd.Timestamp(day)
    folder = (base_dir or FIG_DIR) / station / day.date().isoformat()
    folder.mkdir(parents=True, exist_ok=True)
    raw = load_day(station, day)
    ech = filter_echoes(raw, c)
    tilt = tilt_series(ech, c)
    lat, lon = mstid.station_coords(station)
    terms = solar_terminators(day, lat, lon)
    events = _event_table(station, day)
    logger.info(f"{station} {day.date()}: {raw.file_time.nunique()} soundings, {len(raw)} echoes -> {len(ech)} kept; "
                f"{int(tilt.ew.notna().sum())} tilt samples; sunrise/sunset "
                + ", ".join(f"{k} {t:%H:%M}" for k, t in terms))
    plot_ionosonde(station, day, events, terms, folder / "01_ionosonde_foF2_hmF2.png")
    n_pages = plot_skymaps(ech, station, day, folder)
    plot_zenith_azimuth(ech, tilt, station, day, events, terms, folder / "02_zenith_azimuth.png")
    ridges = plot_periods(tilt, station, day, events, terms, folder / "03_periods_wavelet_LS.png", c)
    tilt.to_csv(folder / "tilt.csv")
    ridges.to_csv(folder / "event_periods.csv", index=False)
    logger.info(f"Saved {n_pages} skymap pages + 3 figures in {folder}")
    return folder


if __name__ == "__main__":
    folder = analyse_event_day(sys.argv[1], pd.Timestamp(sys.argv[2]))
    print(pd.read_csv(folder / "event_periods.csv").to_string(index=False))
