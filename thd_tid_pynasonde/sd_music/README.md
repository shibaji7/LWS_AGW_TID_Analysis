# sd_music — SuperDARN MUSIC (MSTID) analysis

`music_simulate.py` is a script version of pyDARNmusic's `MUSIC-Simulate.ipynb`:
FITACF -> interpolation -> (optional simulated MSTIDs) -> FIR band-pass -> FFT ->
cross-spectral matrix -> MUSIC k-array -> signal detection. Figures and `signals.csv`
go to `Data/figures/sd_music/<radar>_<date>/` (gitignored).

```bash
pip install pyDARNmusic          # env: apep
python music_simulate.py --download-testdata                  # Wallops 2011-05-09 sample, then runs
python music_simulate.py --simulate default                   # pyDARNmusic's default wave
python music_simulate.py --simulate off --radar bks --data-dir /path/to/sd-data \
    --start 2023-04-24T14 --end 2023-04-24T20 --interest 2023-04-24T16 2023-04-24T18 \
    --plot-time 2023-04-24T17                                 # real data, same chain
```

Real data layout: `<data-dir>/<year>/fitacf/<radar>/<YYYYMMDD>.HHMM.<..>.<radar>.fitacf.bz2`.

Notes
- pyDARNmusic 0.2.0 does not pin `pydarnio`; a fresh install gets pydarnio 2.x, which removed
  `SDarnRead`, so `pyDARNmusic.load_fitacf` fails. `load_fitacf_records` here reads with
  `pydarnio.read_fitacf` instead (same file selection) — no environment pinning needed.
- pydarn 4.x keys its radar table by `RadarID`; the script converts the record's `stid`.
- Injected waves (`CUSTOM_SIGS`): `(amp, kx, ky, f, phi, dc)`, kx/ky in rad/km, f in Hz.
