"""MUSIC (MUltiple SIgnal Classification) MSTID analysis with optional simulated waves.

A ``.py`` port of ``notebooks/MUSIC-Simulate.ipynb`` from pyDARNmusic
(https://github.com/w2naf-academia/pyDARNmusic): SuperDARN FITACF -> ``musicArray`` ->
beam/time interpolation -> (optional simulated MSTIDs) -> FIR band-pass -> FFT -> cross-spectral
matrix Dlm -> MUSIC k-array -> signal detection. Every step's figure is written to disk
(``--out-dir``, default ``Data/figures/sd_music/<radar>_<date>/``) and the detected signals to
``signals.csv``; nothing is shown interactively.

    python music_simulate.py --download-testdata                 # fetch the 2011-05-09 Wallops sample
    python music_simulate.py --simulate custom                   # two injected MSTIDs + noise (notebook)
    python music_simulate.py --simulate off --radar bks --data-dir /sd-data ...   # real data, same chain

References: Samson et al. (1990), Bristow et al. (1994), Schmidt (1986).
"""
from __future__ import annotations

import argparse
import datetime as dt
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA_ROOT = ROOT.parent / "Data"
TESTDATA_URL = "https://raw.githubusercontent.com/w2naf-academia/pyDARNmusic/main/TestData"

# Signals injected by --simulate custom: cos(kx*x + ky*y - 2*pi*f*t + phi) * amp + dc
#   (amp, kx [rad/km], ky [rad/km], f [Hz], phi, dc offset)   -- as in the notebook
CUSTOM_SIGS = [
    (5, -0.01, 0.010, 0.0007, 0, 5.0),
    (10, 0.022, -0.023, 0.0004, 0, 5.0),
]


def parse_time(s: str) -> dt.datetime:
    return dt.datetime.fromisoformat(s)


def parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--radar", default="wal", help="3-letter radar code (default wal)")
    p.add_argument("--start", type=parse_time, default=dt.datetime(2011, 5, 9, 8), help="data window start (ISO)")
    p.add_argument("--end", type=parse_time, default=dt.datetime(2011, 5, 9, 19), help="data window end (ISO)")
    p.add_argument("--interest", nargs=2, type=parse_time, default=[dt.datetime(2011, 5, 9, 14), dt.datetime(2011, 5, 9, 16)],
                   metavar=("START", "END"), help="time of interest; the filter edges are trimmed around it")
    p.add_argument("--plot-time", type=parse_time, default=dt.datetime(2011, 5, 9, 14), help="time of the fan plots")
    p.add_argument("--data-dir", type=Path, default=DATA_ROOT / "sd_music" / "TestData",
                   help="FITACF tree: <dir>/<year>/fitacf/<radar>/<YYYYMMDD>.HHMM.*.<radar>.fitacf.bz2")
    p.add_argument("--gates", nargs=2, type=int, default=[30, 45], metavar=("LO", "HI"), help="range-gate limits")
    p.add_argument("--timeres", type=int, default=120, help="time resolution of the interpolated data [s]")
    p.add_argument("--numtaps", type=int, default=101, help="FIR filter taps")
    p.add_argument("--cutoff-low", type=float, default=0.0003, help="band-pass low cut-off [Hz]")
    p.add_argument("--cutoff-high", type=float, default=0.0012, help="band-pass high cut-off [Hz]")
    p.add_argument("--simulate", choices=["off", "default", "custom"], default="custom",
                   help="off: real data only; default: pyDARNmusic's default wave; custom: CUSTOM_SIGS + noise")
    p.add_argument("--noise", type=float, default=4.0, help="noise factor for --simulate custom")
    p.add_argument("--kmax", type=float, default=0.05, help="MUSIC k-array extent kx,ky [rad/km]")
    p.add_argument("--dk", type=float, default=0.001, help="MUSIC k-array step [rad/km]")
    p.add_argument("--out-dir", type=Path, default=None, help="figure/CSV folder")
    p.add_argument("--download-testdata", action="store_true", help="fetch the 2011-05-09 Wallops FITACF sample if missing")
    return p.parse_args(argv)


def download_testdata(data_dir: Path, radar: str = "wal") -> None:
    """The notebook's sample: pyDARNmusic TestData/2011/fitacf/wal, 2-hourly files 06:01-18:01 UT."""
    folder = data_dir / "2011" / "fitacf" / radar
    folder.mkdir(parents=True, exist_ok=True)
    for hhmm in ("0601", "0801", "1001", "1201", "1401", "1601", "1801"):
        name = f"20110509.{hhmm}.00.{radar}.fitacf.bz2"
        if (folder / name).is_file() and (folder / name).stat().st_size > 0:
            continue
        print(f"downloading {name}")
        urllib.request.urlretrieve(f"{TESTDATA_URL}/2011/fitacf/{radar}/{name}", folder / name)


def load_fitacf_records(radar: str, start: dt.datetime, end: dt.datetime, data_dir: Path, fit_sfx: str = "fitacf") -> list[dict]:
    """FITACF records for ``[start, end)`` from bz2 files under ``data_dir/<year>/<fit_sfx>/<radar>/``.

    Same file selection as ``pyDARNmusic.load_fitacf``, but reads with the current ``pydarnio.read_fitacf``:
    pyDARNmusic 0.2.0 still calls ``pydarnio.SDarnRead`` (removed in pydarnio 2.x) and so fails on a fresh install.
    """
    import bz2
    import glob

    import numpy as np
    import pydarnio
    from pydarn import time2datetime

    days = pd.date_range(start.date(), end.date(), freq="D")
    files = []
    for d in days:
        files += glob.glob(str(Path(data_dir) / f"{d:%Y}" / fit_sfx / radar / f"{d:%Y%m%d}*{radar}*.{fit_sfx}.bz2"))
    files = sorted(files)
    t0 = [dt.datetime.strptime(Path(f).name[:13], "%Y%m%d.%H%M") for f in files]
    keep = [(f, t) for f, t in zip(files, t0) if t <= end]
    before = [(start - t).total_seconds() for _, t in keep]
    if any(b >= 0 for b in before):  # the file that already covers ``start`` begins before it
        i = int(np.argmin([b for b in before if b >= 0]))
        keep = keep[i:]
    records = []
    for f, _ in keep:
        try:
            recs, _bad = pydarnio.read_fitacf(bz2.open(f).read())
        except Exception as exc:  # empty/corrupt file: skip it like pyDARNmusic does
            print(f"skipping {f}: {exc}")
            continue
        records += recs
    return [r for r in records if start <= time2datetime(r) < end]


class Figures:
    """Numbered PNG saver: 01_name.png, 02_name.png ..."""

    def __init__(self, folder: Path):
        self.folder, self.n = folder, 0
        folder.mkdir(parents=True, exist_ok=True)

    def save(self, fig, name: str):
        # pyDARNmusic plot calls return a Figure, or an object holding it in ``.fig`` (musicRTP / musicFan), or None
        fig = fig if hasattr(fig, "savefig") else getattr(fig, "fig", None) or plt.gcf()
        self.n += 1
        path = self.folder / f"{self.n:02d}_{name}.png"
        fig.savefig(path, dpi=110, bbox_inches="tight")
        plt.close(fig)
        print(f"  figure {path.name}")
        return path


def run(args: argparse.Namespace) -> Path:
    # Imported here so --help works without the SuperDARN stack.
    import cartopy.crs as ccrs
    from pydarn import RadarID, SuperDARNRadars
    from pyDARNmusic import (beamInterpolation, calculateDlm, calculateFFT, calculateKarr, defineLimits,
                             detectSignals, determineRelativePosition, filterTimes, musicFan,
                             musicRTP, plotDlm, plotFullSpectrum, plotKarr, plotKarrDetected,
                             plotRelativeRanges, simulator, spectrumMultiPlot, timeInterpolation,
                             timeSeriesMultiPlot, music)

    if args.download_testdata:
        download_testdata(args.data_dir, args.radar)
    out = args.out_dir or DATA_ROOT / "figures" / "sd_music" / f"{args.radar}_{args.start:%Y%m%d}"
    figs = Figures(out)

    # --- load + musicArray -------------------------------------------------------------------
    fitacf = load_fitacf_records(args.radar, args.start, args.end, args.data_dir)
    if not fitacf:
        raise FileNotFoundError(f"No FITACF for {args.radar} {args.start}..{args.end} under {args.data_dir} "
                                "(try --download-testdata for the Wallops 2011-05-09 sample)")
    # fovModel='GS': ground-scatter mapping (Bristow et al. 1994); only ground-scatter echoes are used.
    data = music.musicArray(fitacf, sTime=args.start, eTime=args.end, fovModel="GS")
    print("data sets:", data.get_data_sets())

    figs.save(musicRTP(data), "rtp_original")
    figs.save(musicFan(data, time=args.plot_time), "fan_original")

    # --- selection: gate + time limits (the FIR filter eats numtaps/2 samples at each end) ----
    defineLimits(data, gateLimits=list(args.gates))
    defineLimits(data, timeLimits=filterTimes(args.interest[0], args.interest[1], args.timeres, args.numtaps))
    data.active.applyLimits()
    figs.save(musicRTP(data), "rtp_limited")

    # --- interpolation to a gap-free beam x gate x time cube ------------------------------------
    beamInterpolation(data)
    radar_id = fitacf[0]["stid"]
    hw = SuperDARNRadars.radars[RadarID(int(radar_id))].hardware_info.geographic
    proj = ccrs.Orthographic(hw.lon, hw.lat)
    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(121, projection=proj, aspect="auto")
    musicFan(data, plotZeros=True, axis=ax, dataSet="originalFit", time=args.plot_time, autoScale=True)
    ax = fig.add_subplot(122, projection=proj, aspect="auto")
    musicFan(data, plotZeros=True, axis=ax, time=args.plot_time)
    figs.save(fig, "fan_beam_interpolated")
    timeInterpolation(data, timeRes=args.timeres)
    figs.save(timeSeriesMultiPlot(data, dataSet="timeInterpolated", dataSet2="beamInterpolated"), "timeseries_interp")
    determineRelativePosition(data)  # local Cartesian grid (km) of every cell, used for the wave vector
    figs.save(plotRelativeRanges(data, time=args.plot_time), "relative_ranges")

    # --- optional simulated MSTIDs inserted into the processing chain ----------------------------
    if args.simulate == "default":
        simulator(data)
    elif args.simulate == "custom":
        simulator(data, sigs=CUSTOM_SIGS, noiseFactor=args.noise)
    if args.simulate != "off":
        figs.save(musicRTP(data), f"rtp_simulated_{args.simulate}")
        figs.save(musicFan(data, time=args.plot_time), f"fan_simulated_{args.simulate}")

    # --- FIR band-pass --------------------------------------------------------------------------
    filt = music.filter(data, numtaps=args.numtaps, cutoff_low=args.cutoff_low, cutoff_high=args.cutoff_high)
    figs.save(filt.plotTransferFunction(xmax=0.004), "filter_transfer")
    figs.save(filt.plotImpulseResponse(), "filter_impulse")
    data.active.applyLimits()  # drop the filter-edge samples
    figs.save(musicRTP(data, autoScale=True), "rtp_filtered")
    figs.save(musicFan(data, time=args.plot_time, autoScale=True), "fan_filtered")

    # --- spectra, cross-spectral matrix, MUSIC k-array -------------------------------------------
    data.active.printHistory()
    calculateFFT(data)
    figs.save(spectrumMultiPlot(data, xlim=(-0.0025, 0.0025)), "spectrum_real_imag")
    figs.save(spectrumMultiPlot(data, plotType="magnitude", xlim=(0, 0.0025)), "spectrum_magnitude")
    figs.save(spectrumMultiPlot(data, plotType="phase", xlim=(0, 0.0025)), "spectrum_phase")
    figs.save(plotFullSpectrum(data, xlim=(0, 0.00175)), "spectrum_full")
    calculateDlm(data)
    figs.save(plotDlm(data), "dlm")
    calculateKarr(data, kxMax=args.kmax, kyMax=args.kmax, dkx=args.dk, dky=args.dk)
    figs.save(plotKarr(data), "karr")

    # --- signal detection --------------------------------------------------------------------------
    detectSignals(data)
    figs.save(plotKarr(data), "karr_detected")
    figs.save(plotKarrDetected(data), "karr_detection_groups")
    data.active.printHistory()

    info = data.active.sigDetect.info
    sig = pd.DataFrame(info)
    sig.to_csv(out / "signals.csv", index=False)
    print(f"detected {len(sig)} signal(s):")
    print(sig.to_string(index=False))
    return out


if __name__ == "__main__":
    run(parse_args())
